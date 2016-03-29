#include "common.h"

int** init_table_with_delta(int M,int N)
{
  int **table=allocate_2d_rectangle(M,N);
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) { 
      table[i][j] = rand_int(3)-1;
      //table[i][j] = 1;//-1+2*((i+j)%2);
    }
  }
  return table;
}

#define INDEX(e) ((e)<(minenergy)) ? (0) : (((e) > (minenergy+estep*energies)) ? (energies+1) : ((int)ceil((e-minenergy)/estep)))

int wang_landau(int M, int N, double delta, double minenergy, double estep, int energies, double flatness, long Niter, double T0, double Tstep, long Tn) {
  int **table=init_table_with_delta(M,N);
  int estat[energies+2];
  double edensity[energies];
  long i=0;
  int x,y;
  int k;
  int newstate;
  long iter=0;
  long energy, newenergy;
  double maximum, minimum, lnf=1;
  double prob=0;
  //  minenergy=min_energy(M,N,table);
  //  estep=-((double)2*minenergy-0.0001)/energies;
  energy=state_energy(M,N,table);
  //    printf("%ld %ld %lf %d\n",min_energy(M,N,table),energy, estep, (int)floor((energy-minenergy)/estep));    

  for (k=0;k<energies+2;k++) {
    estat[k]=0;
    edensity[k]=0;
  }
    
  while (1) {
    i++;
    x=rand_int(M);
    y=rand_int(N);
    //    energy=state_energy(M,N,table);    
    newenergy=energy;
    newstate=(table[x][y]+1+(nrand()<0.5));
    if (newstate>1) newstate=newstate-3;
    for (k=0; k<lattice_neighbors; k++){
      newenergy+=(table[x][y]-newstate)*table[mod(x+lattice_dx[k],M)][mod(y+lattice_dy[k],N)];
    }
    newenergy+=delta*(newstate*newstate-table[x][y]*table[x][y]);
    
    prob=exp(edensity[INDEX(energy)]-edensity[INDEX(newenergy)]);
    if (nrand()<prob) {
      table[x][y]=newstate;
      energy=newenergy;
      //printf("%lf %ld %ld %lf\n",prob,newenergy,(int)floor((energy-minenergy)/estep),edensity[(int)floor((energy-minenergy)/estep)]);
    }
    estat[INDEX(energy)]+=1;
    edensity[INDEX(energy)]+=lnf;
    //    printf("%ld %lf %d\n",energy,estep,(int)floor((energy-minenergy)/estep));
    if ((i>100000) && (i%10000==0)) {
      maximum=0;
      minimum=1e100;
      for (k=0; k<energies+2; k++) {
	if ((estat[k]!=0) && (edensity[k]>lnf)) {
	  if (estat[k]>maximum) maximum=estat[k];
	  if (estat[k]<minimum) minimum=estat[k];
	}
      }
      if (i%1000000==0) printf("%lf %lf\n",maximum,minimum);
      
      if (2*minimum> (maximum+minimum)*flatness) {
	lnf *=0.5;
	iter++;
	printf("%ld %lf %ld\n",iter,lnf,i);	
	if (iter>=Niter) break;
	for (k=0;k<energies+2;k++)
	  estat[k]=0;
	i=0;
      }
    }
  }
  minimum=1e100;

  for (k=0;k<energies+2; k++) {
    if ((estat[k]!=0) && (edensity[k]<minimum)) minimum=edensity[k];
    printf("%lf %ld %lf\n", minenergy+estep*k, estat[k], edensity[k]);    
  }
  if (iter<Niter) printf("Not enough steps: %ld %lf %ld\n",iter,lnf,i);	
  //  printf("Thermodynamics:\n");
  //  printf("T Ev Cv F:\n");  
//  double w=0, Z=0, Ev=0, E2v=0, cv;
//  double T=T0;
//  for (i=0; i<Tn; i++) {
//    w=0;
//    Z=0;
//    Ev=0;
//    E2v=0;
//    for (k=0;k<energies+2; k++) {
//      w=exp(edensity[k]-minimum+log(2.0)-(minenergy+k*estep)/T);
//      Z+=w;
//      Ev+=w*(minenergy+k*estep);
//      E2v +=w*(minenergy+k*estep)*(minenergy+k*estep);
//    }
//    Ev = Ev/Z;
//    cv=(E2v/Z-Ev*Ev)/(T*T);
//    //    printf("%lf %lf %lf %lf\n", T, Ev/(M*N),cv/(M*N),T*log(Z)/(M*N));
//    T+=Tstep;
//  }
  free_2d_array(M,table);
  //  printf("%lf\n",lnf);
  return 0;
}


int main(int argc, char **argv)
{
  lattice_neighbors=square_lattice_neighbors;
  lattice_dx=square_lattice_dx;
  lattice_dy=square_lattice_dy;
  
  gmp_randinit_mt(state);
  gmp_randseed_ui(state, time(NULL));
  mpf_init2(rop,256);
  //  wang_landau(5,25,1.0,20,0.8,5,1.5, 0.05, 30);
  //  return 0;

  
  double T=1.5;
  double step=0.01;
  int M=50;
  int N=50;
  int steps = 100;
  long measures=10000;
  double delta=1.9655;
  char *method="";

  int  opt;
  double minenergy=(-2+delta)*M*N;
  double estep=0.02;
  long energy_steps=100*M*N;

  double flatness=0.8;
  int iter=5;
  
  while ((opt = getopt(argc, argv, "M:N:T:S:s:d:E:e:k:f:i:L:")) != -1) {
    switch (opt) {
    case 'M':
      M=atoi(optarg);
      break;
    case 'N':
      N=atoi(optarg);
      break;
    case 'T':
      T = atof(optarg);
      break;
    case 'S':
      step = atof(optarg);
      break;
    case 's':
      steps = atoi(optarg);
      break;
    case 'd':
      delta = atof(optarg);
      break;
    case 'E':
      minenergy=atof(optarg);
      break;
    case 'e':
      estep=atof(optarg);
      break;
    case 'k':
      energy_steps=atol(optarg);
      break;
    case 'f':
      flatness=atof(optarg);
      break;
    case 'i':
      iter=atoi(optarg);
      break;
    case 'L':
      if ((strlen(optarg)==6) && (strcmp("triang",optarg)==0)) {
	lattice_neighbors=triangular_lattice_neighbors;
	lattice_dx=triangular_lattice_dx;
	lattice_dy=triangular_lattice_dy;
      }
      break;
      
    default: /* '?' */
      fprintf(stderr, "Usage: %s [-M xsize] [-N ysize] [-T temperature] [-S temperature step] [-s number of steps] [-d delta] [-E min_energy] [-e estep] [-k energy_steps] [-f flatness] [-i iterations] [-L lattice:square|triang]  name\n",
	      argv[0]);
      exit(EXIT_FAILURE);
    }
  }

  wang_landau(M,N,delta,minenergy, estep, energy_steps, flatness, iter, T, step, steps);
  return 0;
}
