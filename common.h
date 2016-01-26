#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <gmp.h>
#include <getopt.h>


gmp_randstate_t state;
mpf_t rop;

int square_lattice_dx[]={-1,0,+1,0};
int square_lattice_dy[]={0,-1,0,+1};
int square_lattice_neighbors=4;


int triangular_lattice_dx[]={-1,0,+1,0,-1,1};
int triangular_lattice_dy[]={0,-1,0,+1,-1,1};
int triangular_lattice_neighbors=6;


int lattice_neighbors;
int* lattice_dx;
int* lattice_dy;


int **allocate_2d_array(int size)
{
  int **table;
  table = (int**)malloc(size * sizeof(int*));
  for (int i = 0; i < size; i++) {
    table[i] = (int*)malloc(size * sizeof(int));
  }
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      table[i][j] = 0;
    }
  }
  return table;
}

int **allocate_2d_rectangle(int M, int N)
{
  int **table;
  table = (int**)malloc(M * sizeof(int*));
  for (int i = 0; i < M; i++) {
    table[i] = (int*)malloc(N * sizeof(int));
  }
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      table[i][j] = 0;
    }
  }
  return table;
}

void free_2d_array(int size, int **table){
  for (int i = 0; i < size; i++) {
    free(table[i]);
  }
  free(table);
}

int mod(int x, int m) {
    return (x%m + m)%m;
}

int rand_int(int max) 
{  
  return (int)floor(((double) max * gmp_urandomm_ui(state,(long)RAND_MAX+1)) / (RAND_MAX+1.0));
}

double nrand()
{
  mpf_urandomb (rop,  state, 256);
  return mpf_get_d(rop);
  //  return (double)(gmp_urandomm_ui(state,(long)RAND_MAX+1)) / (double)RAND_MAX ;
}

void init_table_pm_boundary(int M, int N, int ** table)
{
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) { 
      table[i][j] = 2 * (nrand() > 0.5) - 1;
      if (( i == 0 ) || ( ( j == N-1 ) && ( i < M - 1)))
	table[i][j]=-1;
      if (( i == N-1 ) || (( j == 0 ) && ( i > 0 )))
	table[i][j]=1;
    }
  }
}

void init_table_no_boundary(int M,int N, int ** table)
{
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) { 
      table[i][j] = 2 * (nrand() > 0.5) - 1;
      //table[i][j] = 1;//-1+2*((i+j)%2);
    }
  }
}

long calc_magnetization(int M, int N, int ** table)
{
  long res=0;
  for (int i=0; i<M; i++)
    for (int j=0; j<N; j++)
      res+=table[i][j];
  return res;
}

long state_energy(int M, int N, int** table) {
  int i, j,k;
  long e=0;

  for (i=0; i<M; i++) 
    for (j=0;j<N; j++) 
        for (k=0; k<lattice_neighbors; k++)
	  e-=table[i][j]*table[mod(i+lattice_dx[k],M)][mod(j+lattice_dy[k],N)];
  return e/2;
}

int** init_table_with_p(int M, int N, double p){
  int **table=allocate_2d_rectangle(M,N);
  long num=(long)((1.0-p)*M*N);
  long i=0,x,y;
  init_table_no_boundary(M,N,table);
  while (i<num) {
    x=rand_int(M);
    y=rand_int(N);
    if (table[x][y]!=0) {
      table[x][y]=0;
      i++;
    }
  }
  //  printf("%ld %ld %lf %lf\n", (long)size*size,num, (double)num/(size*size),p);
  return table;
}

long min_energy(int M, int N, int** table) {
  int i, j,k;
  long e=0;

  for (i=0; i<M; i++) 
    for (j=0;j<N; j++) 
        for (k=0; k<lattice_neighbors; k++)
	  e-=abs(table[i][j]*table[mod(i+lattice_dx[k],M)][mod(j+lattice_dy[k],N)]);
  return e/2;
}


void init_with_p(int M, int N, int **table, double p){
  long num=(long)((1.0-p)*M*N);
  long i=0,x,y;
  int j;
  for ( i = 0; i < M; i++) {
    for ( j = 0; j < N; j++) { 
      //      table[i][j] =  2*((i+j)%2)-1;
      table[i][j]=2 * (nrand() > 0.5) - 1;
    }
  }
  i=0;
  while (i<num) {
    x=rand_int(M);
    y=rand_int(N);
    if (table[x][y]!=0) {
      table[x][y]=0;
      i++;
    }
  }
}
