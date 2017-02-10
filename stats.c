#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <getopt.h>

int main(int argc, char **argv)
{
  double from=0;
  double to=1000;
  double step=0.1;
  long steps;
  int flags, opt;
  
  while ((opt = getopt(argc, argv, "f:t:s:")) != -1) {
    switch (opt) {
    case 'f':
      from=atof(optarg);
      break;
    case 't':
      to=atof(optarg);
      break;
    case 's':
      step = atof(optarg);
      break;
      
    default: /* '?' */
      fprintf(stderr, "Usage: %s [-f from time] [-t to time] [-s time step]\n",
	      argv[0]);
      exit(EXIT_FAILURE);
    }
  }

  steps=(long)((to-from)/step)+1;
  long *allvisits=(long *)calloc(steps,sizeof(long));
  double *allxin=(double*)calloc(steps,sizeof(double));
  double *allyn=(double*)calloc(steps,sizeof(double));
  double *allxinsquare=(double*)calloc(steps,sizeof(double));
  double *allynsquare=(double*)calloc(steps,sizeof(double));
  double *allxinyn=(double*)calloc(steps,sizeof(double));

  long *visits=(long *)calloc(steps,sizeof(long));
  double *xin=(double*)calloc(steps,sizeof(double));
  double *yn=(double*)calloc(steps,sizeof(double));
  double *xinsquare=(double*)calloc(steps,sizeof(double));
  double *ynsquare=(double*)calloc(steps,sizeof(double));
  double *xinyn=(double*)calloc(steps,sizeof(double));

  
  int err;

  double t,xi,y;
  long index,i;

  err=scanf("{");
  while(err!=EOF){
    err=scanf("{");
    while((err=scanf("{%lf,%lf,%lf},",&t,&xi,&y))==3) {
      index=(long)(floor((t-from)/step));
      //      fprintf(stderr,"%ld %.10g %.10g %.10g\n",index,t,xi,y);      
      if ((index>=0) && (index<steps)) {
	visits[index]++;
	xin[index]+=xi;
	yn[index]+=y;
	xinsquare[index]+=xi*xi;
	ynsquare[index]+=y*y;
	xinyn[index]+=xi*y;
      }
    }
    err=scanf("}, ");
  }

  for (i=0;i<steps;i++) {
    if (visits[i]>0) {
      printf("%.10g %ld %.10g %.10g %.10g %.10g %.10g\n",from+step*i,visits[i],xin[i]/visits[i],xinsquare[i]/visits[i],yn[i]/visits[i],ynsquare[i]/visits[i],xinyn[i]/visits[i]);
    }
  }

  free(visits);
  free(xin);
  free(yn);
  free(xinsquare);
  free(ynsquare);
  free(xinyn);
  
  return 0;
}
