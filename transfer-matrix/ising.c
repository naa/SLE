#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <getopt.h>

int main(int argc, char **argv)
{
  
  double T=1.5;
  double step=0.01;
  int M=10;
  int N=50;
  int steps = 100;

  int flags, opt;
  int nsecs, tfnd;

  while ((opt = getopt(argc, argv, "M:N:T:S:s:")) != -1) {
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
        
    default: /* '?' */
      fprintf(stderr, "Usage: %s [-M xsize] [-N ysize] [-T temperature] [-S temperature step] [-s number of steps]\n",
	      argv[0]);
      exit(EXIT_FAILURE);
    }
  }

  double *lnZ;
  //  double *lnZbar;
  int twotoM;
  long double Z;
  twotoM=1<<M;
  lnZ=(double*)malloc(twotoM*sizeof(double));
  //  lnZbar=(double*)malloc(twotoM*sizeof(double));  
  int i,j,k,l;
  for (i=0;i<steps;i++) {
    for (j=0; j<twotoM;j++){
      lnZ[j]=0;
      for (l=0;l<M-1;l++) {
	lnZ[j]+=2/T*(((j>>l)%2)^(((~j)>>(l+1))%2));
      }
      lnZ[j]+=2/T*(((j>>(M-1))%2)^((~j)&1))-M/T;
      //      printf("%d %lf\n",j,lnZ[j]);
    }
    for (k=0;k<N; k++) {
      for (l=0;l<M;l++) {
	for (j=0;j<twotoM;j++) {
	  lnZ[j]=1/T+lnZ[j]+log(1+exp(-2/T+lnZ[j^(1<<l)]-lnZ[j]));
	}
      }
      for (j=0; j<twotoM;j++){      
	for (l=0;l<M-1;l++) {
	  lnZ[j]+=2/T*(((j>>l)%2)^(((~j)>>(l+1))%2));
	}
	lnZ[j]+=2/T*((((~j)>>(M-1))%2)^((~j)&1))-M/T;
      }
    }
    Z=0.0;
    for (j=0;j<twotoM;j++) {
      Z+=expl(lnZ[j]);
    }
    printf("%d %lf %lf\n",M,T,-(double)(logl(Z)*T/(M*N)));
    T=T+step;
  }
  return 0;
}
