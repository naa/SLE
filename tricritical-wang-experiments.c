#include "common.h"

int** init_table_with_delta(int M,int N)
{
  int **table=allocate_2d_rectangle(M,N);
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) { 
      //      table[i][j] = rand_int(3)-1;
      table[i][j] = (rand_int(10)>1)? 1:0;
      //table[i][j] = 1;//-1+2*((i+j)%2);
    }
  }
  return table;
}

long calc_energy(int M,int N,int **table,long delta) {
  int i, j,k;
  long e=0;

  for (i=0; i<M; i++) 
    for (j=0;j<N; j++) {
        for (k=0; k<lattice_neighbors; k++)
	  e-=table[i][j]*table[mod(i+lattice_dx[k],M)][mod(j+lattice_dy[k],N)]*10000;
	e+=2*delta*table[i][j]*table[i][j];
    }
  return e/2;
}

#define INDEX(e) ((e)<(minenergy)) ? (0) : (((e) > (minenergy+estep*energies)) ? (energies+1) : ((e-minenergy)/estep))

int wang_landau(int M, int N, long delta, long long minenergy, long estep, long energies, double flatness, long Niter, double T0, double Tstep, long Tn) {
  int **table=init_table_with_delta(M,N);
  long estat[energies+2];
  double edensity[energies];
  long i=0;
  int x,y;
  int k;
  int newstate;
  long iter=0;
  long long energy, newenergy;
  long nume;
  double maximum, minimum, lnf=1;
  double total, mean;
  double prob=0;
  //  minenergy=min_energy(M,N,table);
  //  estep=-((double)2*minenergy-0.0001)/energies;
  energy=calc_energy(M,N,table,delta);
  while ((INDEX(energy)<0) || (INDEX(energy)>energies)) {
    fprintf(stderr,"energy: %ld, index: %ld\n",energy,INDEX(energy));
    free_2d_array(M,table);    
    table=init_table_with_delta(M,N);
    energy=calc_energy(M,N,table,delta);
  }
  //    printf("%ld %ld %lf %d\n",min_energy(M,N,table),energy, estep, (int)floor((energy-minenergy)/estep));    


  for (k=0;k<energies+2;k++) {
    estat[k]=0;
    edensity[k]=0;
  }

  //  estat[0]=1;
  //  edensity[0]=lnf;
  
  while (1) {

    x=rand_int(M);
    y=rand_int(N);
    newenergy=energy;
    //    energy=state_energy(M,N,table);    
    newstate=(table[x][y]+1+(nrand()<0.5));
    if (newstate>1) newstate=newstate-3;
    for (k=0; k<lattice_neighbors; k++){
      newenergy+=(table[x][y]-newstate)*table[mod(x+lattice_dx[k],M)][mod(y+lattice_dy[k],N)]*10000;
    }
    newenergy+=delta*(newstate*newstate-table[x][y]*table[x][y]);
    if (INDEX(newenergy)>energies) continue;

    i++;    
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
      //      for (k=0; k<energies+2; k++) {
      total=0;
      nume=0;
      long first=-1;
      long last=0;
      long minindex=0;
      for (k=0;k<energies+1 ; k++) {
	if (edensity[k]>0) {// || ((k==0) && (edensity[k]==0))) {
	  if (first==-1) first=k;
	  last=k;
	  total+=estat[k];
	  nume++;
	  //	  if (estat[k]>maximum) maximum=estat[k];
	  if (estat[k]<minimum) { minimum=estat[k]; minindex=k; }
	}
      }

      if (i%1000000==0)   fprintf(stderr,"i:%ld, %ld steps, %lf %lf, %ld of %ld, first: %ld, last: %ld, min energy: %lf, max energy: %lf, min index: %ld\n",iter, i,total/nume,minimum, nume, energies,first,last,minenergy/10000.0+estep*first/10000.0,minenergy/10000.0+estep*last/10000.0,minindex);
      
      if (minimum> total/nume*flatness) {
	lnf *=0.5;
	iter++;
	fprintf(stderr, "iteration %ld %lf %ld\n",iter,lnf,i);	
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
    if (edensity[k]>0)
      printf("%lf %ld %lf\n", (minenergy+estep*k)/10000.0, estat[k], edensity[k]);    
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
  long delta=19655;
  char *method="";

  int  opt;
  long long minenergy=0;
  long estep=0;
  long energy_steps=0;

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
      delta = atoi(optarg);
      break;
    case 'E':
      minenergy=atol(optarg);
      break;
    case 'e':
      estep=atol(optarg);
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
  if (minenergy==0) 
    minenergy=M*N*(-20000+delta);
  //  long long maxenergy=M*N*(20000+delta);
  long long maxenergy=-(minenergy*10);
  if ((estep==0) && (energy_steps!=0))
    estep=(maxenergy-minenergy)/energy_steps;
  else if ((estep!=0) && (energy_steps==0))
    energy_steps=(long)((maxenergy-minenergy)/estep);
  else if ((estep!=0) && (energy_steps!=0)) {
    maxenergy=minenergy+estep*energy_steps;
  } else {
    estep=10000;
    energy_steps=(long)((maxenergy-minenergy)/estep)+1;
  }
  fprintf(stderr,"mine:%ld, estep:%ld, energy_steps:%ld, maxe:%ld\n", minenergy, estep, energy_steps, maxenergy);
  wang_landau(M,N,delta,minenergy, estep, energy_steps, flatness, iter, T, step, steps);
  return 0;
}
