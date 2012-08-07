#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

int **allocate_2d_array(int size) {
  int **table;
  table = (int**)malloc(size * sizeof(int*));
  for (int i = 0; i < size; i++) {
    table[i] = (int*)malloc(size * sizeof(int));
  }
  return table;
}

void free_2d_array(int size, int **table){
  for (int i = 0; i < size; i++) {
    free(table[i]);
  }
  free(table);
}


int rand_int(int max) 
{  
  return (int)floor(((double) max * rand()) / (RAND_MAX+1.0));
}

double nrand()
{
  return (double)rand() / (double)RAND_MAX ;
}


void init_table(int size, int **table)
{
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) { 
      table[i][j] = 2 * (nrand() > 0.5) - 1;
      if ((i == 0) || (i == size - 1) || (j == 0) || (j == size - 1)) {
	table[i][j] = 2 * (2 * i > size) - 1;
      }
    }
  }
}
 
int modify_cell (int size,int **table, int i, int j, double beta)
{
  double p1 = exp( beta * (table[i - 1][j] + table[i + 1][j] + table[i][j - 1] + table[i][j + 1]));
  double p2 = exp(-beta * (table[i - 1][j] + table[i + 1][j] + table[i][j - 1] + table[i][j + 1]));
  int nval = 2 * (((p1 + p2) * nrand()) < p1) - 1;
  int res = (nval != table[i][j]);
  table[i][j]=nval;
  return res;
}

void evolve_table ( int size, int** table, double beta, long long niter ) 
{
  int i,j;
  init_table(size,table);

  for (long long k = 0; k < niter; k++){
    i=1 + rand_int(size-2);
    j=1 + rand_int(size-2);
    modify_cell (size,table,i, j,beta);
  }
}

int ** calc_interface(int size, int **table)
{
  int **interface = allocate_2d_array(size);
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      interface[i][j] = 0;
    }
  }

  int x1, y1, x2, y2, dx, dy; // dx = x2 - x1; dy = y2 - y1;
  int c1, c2;

  for (x1 = 0; (x1 < size - 1) && (table[x1][0] == table[x1 + 1][0]); x1++) {}
  x2 = x1 + 1; 
  y1 = 0; 
  y2 = 0; 

  for (int k = 0; ((y2 < size - 1) || (y1 < size - 1) ) && (k < size * size); k++)    {
    interface[x1][y1] = 1;
    interface[x2][y2] = 2;
    dx = x2 - x1; dy = y2 - y1;
    c1 = table[x1 - dy][y1 + dx];
    c2 = table[x2 - dy][y2 + dx];
    if ((3 * c1 + c2) == (3 * (+1) + (+1))) {x2 = x1 - dy; y2 = y1 + dx;}
    if ((3 * c1 + c2) == (3 * (-1) + (-1))) {x1 = x2 - dy; y1 = y2 + dx;}
    if ((3 * c1 + c2) == (3 * (-1) + (+1))) {x1 = x1 - dy; y1 = y1 + dx; x2 = x2 - dy; y2 = y2 + dx;}
    //      if ((3 * c1 + c2) == (3 * (+1) + (-1))) {x1 = x2 - dy; y1 = y2 + dx;}
    if ((3 * c1 + c2) == (3 * (+1) + (-1))) {
      if (nrand() <0.5) {
	x1 = x2 - dy; y1 = y2 + dx;
      } else {
	x2 = x1 - dy; y2 = y1 + dx;
      }
    }
  }
  interface[x1][y1] = 1;
  interface[x2][y2] = 2;
  return interface;
}

void print_table (FILE* F, int size, int **table)
{
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      if (table[i][j] > 0) 
	fprintf(F, "+"); 
      else
	fprintf(F, "-");
    }
    fprintf(F, "\n");
  }
}


void print_table_xpm (FILE* F, int size, int **table, int **interface)
{
  fprintf(F, "/* XPM */\nstatic char *pict[] = {\n");
  fprintf(F, "\"%d %d %d %d\",\n", size, size, 4, 1);
  fprintf(F, "\"+ c #ffffff\",\n");
  fprintf(F, "\"- c #000000\",\n");
  //fprintf(F, "\"2 c #ffcfcf\",\n");
  //fprintf(F, "\"1 c #300000\",\n");
  fprintf(F, "\"2 c #ff0000\",\n");
  fprintf(F, "\"1 c #ff0000\",\n");

  for (int i = 0; i < size; i++) {
    fprintf(F, "\"");
    for (int j = 0; j < size; j++)
      if (interface[i][j]==0) {
	if (table[i][j] > 0) 
	  fprintf(F, "+"); 
	else 
	  fprintf(F, "-");
      } else {
	fprintf(F, "%d", interface[i][j]);
      }
    fprintf(F, "\"");
    if (i < size - 1) {
      fprintf(F, ",");
    }
    fprintf(F, "\n");
  }

  fprintf(F, "};");
}
 

int main(int argc, char **argv) 
{
  int size=512;
  double beta=2;
  long long niter=size*size;
  FILE* Fout=stdout;

  if (argc>1) size= atoi(argv[1]);
  if (argc>2) beta=atof(argv[2]);
  if (argc>3) niter=atoll(argv[3]);
  if (argc>4) Fout = fopen(argv[4],"w");

  srand(time(NULL));
  
  int **table=allocate_2d_array(size);

  evolve_table(size,table,beta, niter);
  int **interface=calc_interface(size, table);

  print_table_xpm(Fout, size, table, interface);
  fclose(Fout);

  free_2d_array(size,table);
  free_2d_array(size,interface);
  return 0;
}
