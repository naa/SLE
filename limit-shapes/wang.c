#include "common.h"

#define INDEX(e) ((e)<(minenergy)) ? (0) : (((e) > (minenergy+estep*energies)) ? (energies+1) : ((int)ceil((e-minenergy)/estep)))

int contsignal=0;

#include<signal.h>
#include<unistd.h>

void sig_handler(int signo)
{
  if (signo == SIGINT) {
    fprintf(stderr,"received SIGINT\n");
    contsignal=1;
  }
}


void init_with_dimers(int M, int N, int **table, int K){
  int i,j;
  int nv;
  table[0][0]=K;
  table[M-1][N-1]=0;
  for (i=0;i<M;i++) {
    for (j=0;j<N;j++) {
      if (((i==0) &&(j==0)) || ((i==M-1) && (j==N-1))) nv=table[i][j];
      else if (i==0)  nv=1+rand_int(table[i][j-1] < K-1 ? table[i][j-1]: K-1);
      else if (j==0)  nv=1+rand_int(table[i-1][j] < K-1 ? table[i-1][j]: K-1);
      else nv=rand_int(1+(table[i-1][j]<table[i][j-1] ? table[i-1][j]:table[i][j-1]));
      table[i][j]=nv;
    }
  }
}

void init_with_min_dimers(int M, int N, int **table, int K){
  int i,j;
  int nv;
  table[0][0]=K;
  table[M-1][N-1]=0;
  for (i=0;i<M;i++) {
    for (j=0;j<N;j++) {
      if (((i==0) &&(j==0)) || ((i==M-1) && (j==N-1))) nv=table[i][j];
      else if (i==0)  nv=1;
      else if (j==0)  nv=1;
      else nv=0;
      table[i][j]=nv;
    }
  }
}

long state_energy_dimers(int M, int N, int** table) {
  int i, j;
  long e=0;

  for (i=0; i<M; i++) 
    for (j=0;j<N; j++)
      if (table[i][j]>=0)
	e+=table[i][j];
  return e;
}

void print_table (FILE* F, int M, int N, int **table)
{
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
	fprintf(F, "%d ",table[i][j]); 
    }
    fprintf(F, "\n");
  }
    fprintf(F, "\n");  
}

void read_table (FILE* F, int M, int N, int **table)
{
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      fscanf(F, " %d",&(table[i][j])); 
    }
  }
}


int wang_landau(int M, int N, int **table, int K, long minenergy, long estep, long energies, double flatness, long Niter, double T0, double Tstep, long Tn) {

  print_table(stderr,M,N,table);
  long estat[energies+2];
  double edensity[energies];
  long i;
  int x,y;
  int k;
  long iter=0;
  long energy, newenergy;
  int newheight;
  double maximum, minimum, lnf=1;
  double prob=0;
  //  minenergy=min_energy(M,N,table);
  //  estep=-((double)2*minenergy-0.0001)/energies;
  energy=state_energy_dimers(M,N,table);

  //    printf("%ld %ld %lf %d\n",min_energy(M,N,table),energy, estep, (int)floor((energy-minenergy)/estep));    

  for (k=0;k<energies+2;k++) {
    estat[k]=0;
    edensity[k]=0;
  }

  //  estat[INDEX(energy)]+=1;
  //  edensity[INDEX(energy)]+=lnf;
  
  for (i=1; ; i++) {
    //    while(1) {
      x=rand_int(M);
      y=rand_int(N);

      newenergy=energy+2 * (nrand() > 0.5) - 1;
      newheight=table[x][y]+newenergy-energy;

      //      if ((newheight<0)|| (newheight>K)) continue;
      //      else break;

//      if (((x==0)&&(y==0)) || ((x==M-1) && (y==N-1))) continue; 
//      if ((x==0) && (y==N-1) && ((newheight<=0)|| (newheight>=K))) continue;
//      if ((x==M-1) && (y==0) && ((newheight<=0)|| (newheight>=K))) continue;
//
      
      //      if ((x==0) && (newheight<=table[x][y-1]) && (newheight>=table[x][y+1]) && (newheight>=table[x+1][y])) break;
      //      if ((y==0) && (x==M-1) && (newheight<=table[x-1][y]) && (newheight>=table[x][y+1])) break;
      //      if ((y==0) && (newheight<=table[x-1][y]) && (newheight>=table[x+1][y]) && (newheight>=table[x][y+1])) break;       

      //    }
      if ((newheight<0)|| (newheight>K))
	prob=0;
      else if (((x==0) || (table[x-1][y]<0) || (newheight<=table[x-1][y])) &&
	       ((y==0) || (table[x][y-1]<0) || (newheight<=table[x][y-1])) &&
	       ((x==M-1) || (table[x+1][y]<0) || (newheight>=table[x+1][y])) &&
	       ((y==N-1) || (table[x][y+1]<0) || (newheight>=table[x][y+1]))) 
	prob=exp(edensity[INDEX(energy)]-edensity[INDEX(newenergy)]);
      else
	prob=0;
    //fprintf(stderr, "%lf %ld %ld %lf\n",prob,newenergy,(int)floor((energy-minenergy)/estep),edensity[(int)floor((energy-minenergy)/estep)]);
    
    if (nrand()<prob) {
      table[x][y]=newheight;
      energy=newenergy;
      //      print_table(stderr,M,N,table);      
    }
    estat[INDEX(energy)]+=1;
    edensity[INDEX(energy)]+=lnf;
    //    printf("%ld %lf %d\n",energy,estep,(int)floor((energy-minenergy)/estep));


    if ((i>100*M*N) && (i%(10*M*N)==0)) {

    //    if ((i>10) && (i%10==0)) {      
      maximum=0;
      minimum=1e100;
      //      for (k=0; k<energies+2; k++) {
      double total=0;
      long nume=0;
      long first=-1;
      long last=0;
      long minindex=0;
      for (k=0;k<energies+1 ; k++) {
	//	if (edensity[k]>0) {// || ((k==0) && (edensity[k]==0))) {
	if (estat[k]>0) {
	  if (first==-1) first=k;
	  last=k;
	  total+=(double)(estat[k]);
	  nume++;
	  //	  if (estat[k]>maximum) maximum=estat[k];
	  if (estat[k]<minimum) { minimum=estat[k]; minindex=k; }
	}
      }

      if (i%(1000*M*N)==0)   fprintf(stderr,"i:%ld, %ld steps, %lf %lf, %ld of %ld, first: %ld, last: %ld, min energy: %ld, max energy: %ld, min index: %ld\n",iter, i,total/((double)nume),minimum, nume, energies,first,last,minenergy+estep*first,minenergy+estep*last,minindex);
      
      if ((minimum> total/nume*flatness) || (contsignal)) {
	contsignal=0;
	lnf *=0.5;
	iter++;
	fprintf(stderr, "iteration %ld %lf %ld\n",iter,lnf,i);	
	if (iter>=Niter) break;
	for (k=0;k<energies+2;k++)
	  estat[k]=0;
	i=0;
      }
    }
      
//    if ((i>1000000) && (i%10000==0)) {
//      maximum=0;
//      minimum=1e100;
//      for (k=0; k<energies+2; k++) {
//	if ((estat[k]!=0) || (edensity[k]>0.00001)) {
//	  if (estat[k]>maximum) maximum=estat[k];
//	  if (estat[k]<minimum) minimum=estat[k];
//	}
//      }
//      //      printf("%lf %lf\n",mean,minimum);
//      
//      if (2*minimum> (maximum+minimum)*flatness) {
//	lnf *=0.5;
//	iter++;
//	//	printf("%ld %lf %ld\n",iter,lnf,i);	
//	if (iter>=Niter) break;
//	for (k=0;k<energies+2;k++)
//	  estat[k]=0;
//      }
//    }
  }
  minimum=1e100;

  for (k=0;k<energies+2; k++) {
    if ((estat[k]!=0) && (edensity[k]<minimum)) minimum=edensity[k];
    printf("%ld %ld %lf\n", minenergy+estep*k, estat[k], edensity[k]);    
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
  int K=50;
  int steps = 100;
  long measures=10000;
  char *method="";

  int  opt;
  double minenergy=-M*N*2;
  double estep=4;
  long energy_steps=M*N;

  double flatness=0.8;
  int iter=5;
  char* fname="";
  
  while ((opt = getopt(argc, argv, "M:N:K:T:S:s:E:e:k:f:i:L:F:")) != -1) {
    switch (opt) {
    case 'M':
      M=atoi(optarg);
      break;
    case 'N':
      N=atoi(optarg);
      break;
    case 'K':
      K = atoi(optarg);
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
    case 'E':
      minenergy=atoi(optarg);
      break;
    case 'e':
      estep=atoi(optarg);
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
    case 'F':
      fname=optarg;
      break;
    default: /* '?' */
      fprintf(stderr, "Usage: %s [-M xsize] [-N ysize] [-K zsize] [-T temperature] [-S temperature step] [-s number of steps] [-E min_energy] [-e estep] [-k energy_steps] [-f flatness] [-i iterations] [-L lattice:square|triang] [-F initial_data_file_name] \n",
	      argv[0]);
      exit(EXIT_FAILURE);
    }
  }
  //  fprintf(stderr, "%d %d %d %lf %lf %ld %lf %d %lf %lf %d",M,N,K,minenergy, estep, energy_steps, flatness, iter, T, step, steps);
  if (signal(SIGINT, sig_handler) == SIG_ERR)
    fprintf(stderr,"\ncan't catch SIGINT\n");


  int **table=allocate_2d_rectangle(M,N);
  FILE* input=fopen(fname,"r");
  if (input!=NULL)
    read_table(input,M,N,table);
  else
    init_with_dimers(M,N,table,K);
  wang_landau(M,N,table,K,minenergy, estep, energy_steps, flatness, iter, T, step, steps);
  return 0;
}
