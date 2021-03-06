#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#define SWAP(x,y) {float* temp; temp = x; x = y; y = temp;}

/*** Skeleton for Lab 2 ***/

/***** Globals ******/
float *a; /* The coefficients */
float *x;  /* The unknowns */
float *b;  /* The constants */
float err; /* The absolute relative error */
int num = 0;  /* number of unknowns */


/****** Function declarations */
void check_matrix(); /* Check whether the matrix will converge */
void get_input();  /* Read input from file */
int serial_solver(); /* Serial version of the solving method */
int parallel_solver(); /* Parallel version of the solving method */

/********************************/



/* Function definitions: functions are ordered alphabetically ****/
/*****************************************************************/

/* 
   Conditions for convergence (diagonal dominance):
   1. diagonal element >= sum of all other elements of the row
   2. At least one diagonal element > sum of all other elements of the row
 */
void check_matrix()
{
  int bigger = 0; /* Set to 1 if at least one diag element > sum  */
  int i, j;
  float sum = 0;
  float aii = 0;
  
  for(i = 0; i < num; i++)
  {
    sum = 0;
    aii = fabs(a[i * num + i]);
    
    for(j = 0; j < num; j++)
       if( j != i)
	 sum += fabs(a[i * num + j]);
       
    if( aii < sum)
    {
      printf("The matrix will not converge.\n");
      exit(1);
    }
    
    if(aii > sum)
      bigger++;
    
  }
  
  if( !bigger )
  {
     printf("The matrix will not converge\n");
     exit(1);
  }
}


/******************************************************/
/* Read input from file */
/* After this function returns:
 * a[][] will be filled with coefficients and you can access them using a[i][j] for element (i,j)
 * x[] will contain the initial values of x
 * b[] will contain the constants (i.e. the right-hand-side of the equations
 * num will have number of variables
 * err will have the absolute error that you need to reach
 */
void get_input(char filename[], int rank)
{
  FILE * fp;
  int i,j;  

  fp = fopen(filename, "r");
  if(!fp)
  {
    printf("Cannot open file %s\n", filename);
    exit(1);
  }

  fscanf(fp,"%d ",&num);
  fscanf(fp,"%f ",&err);

  /* Now, time to allocate the matrices and vectors */
  a = (float*)malloc(num * num * sizeof(float));
  if(!a)
  {
    printf("Cannot allocate a!\n");
    exit(1);
  }

  x = (float *) malloc(num * sizeof(float));
  if(!x)
  {
	 printf("Cannot allocate x!\n");
	 exit(1);
  }

  b = (float *) malloc(num * sizeof(float));
  if(!b)
  {
	 printf("Cannot allocate b!\n");
	 exit(1);
  }

  /* Now .. Filling the blanks */ 
 
  /* The initial values of Xs */
  for(i = 0; i < num; i++)
    fscanf(fp,"%f ", &x[i]);

  for(i = 0; i < num; i++)
  {
   for(j = 0; j < num; j++)
     fscanf(fp,"%f ",&a[i * num + j]);
   
   /* reading the b element */
   fscanf(fp,"%f ",&b[i]);
  }

  fclose(fp); 
}



/******************************************************/
/* Implement the described method in serial computation.
 * Just for me to get a sense of the computation method. 
 */
int serial_solver()
{
  float *new_x = (float *) malloc(num * sizeof(float));
  float curr_err;
  int nit = 0;

  do 
  {
    nit++;
    curr_err = 0;

    for(size_t i=0; i<num; i++)
    {
      new_x[i] = b[i];
      for(size_t j=0; j<num; j++) 
      {
        if(j != i)
        {
          new_x[i] -= a[i * num + j] * x[j];
        }
      }
      new_x[i] /= a[i * num + i];
    }

    for(size_t i=0; i<num; i++) 
    {
      //printf("curr err: %f cal: %f  err: %f\n", curr_err, fabs((new_x[i] - x[i]) / new_x[i]), err); 
      curr_err = MAX(curr_err, fabs((new_x[i] - x[i]) / new_x[i]));

      x[i] = new_x[i];
    }
    
  } while(curr_err > err);

  free(new_x);
  return nit;
}

/******************************************************/
/* Calculate the error for each variable
 * And output the largest one 
 */
float curr_max_err(float* new_x, float* old_x)
{
  float max_err = 0;
  for(size_t i=0; i<num; i++) 
  {
    //printf("max err: %f new_x: %f  old_x:%f\n", max_err, new_x[i], old_x[i]);
    max_err = MAX(max_err, fabs((new_x[i] - old_x[i]) / new_x[i]));
  }
  return max_err;
}

/******************************************************/
/* Parallel by MPI
 * Idea: 
 * 1. We split the data a, b, x across the processes. So they keep a_local, b_local and x_local
 * 2. Each round, each process calculates the new variable value in their assigned x_local from the old_x. 
 * 3. Then, the newly calculated x (x_local) from each process have to gather up and share across all processes for next round of calculation.
 * 4. Repeat 2 & 3 until the solution converge (compare old_x and gathered new_x).
 * 5. Place the converged solution to x array and output at root process.
 */
int parallel_solver(int comm_size, int rank) 
{
  //broadcast the number of variables and the absolute relative error
  MPI_Bcast(&num, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&err, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  int *a_sendcounts;    // array describing how many elements to send to each process
  int *a_displs;
  int *bx_sendcounts;  
  int *bx_displs;

  a_sendcounts = malloc(sizeof(int)*comm_size);
  a_displs = malloc(sizeof(int)*comm_size);
  bx_sendcounts = malloc(sizeof(int)*comm_size);
  bx_displs = malloc(sizeof(int)*comm_size);

  // calculate send counts and displacements
  int rem = num % comm_size;
  int sum = 0;
  for (int i = 0; i < comm_size; i++) {
      a_sendcounts[i] = (num/comm_size) * num;
      bx_sendcounts[i] = num/comm_size;
      if (rem > 0) {
          a_sendcounts[i] += num;
          bx_sendcounts[i]++;
          rem--;
      }

      a_displs[i] = sum * num;
      bx_displs[i] = sum;
      sum += bx_sendcounts[i];
  }

  //for debug
  /*if (0 == rank) {
    for (int i = 0; i < comm_size; i++) {
      printf("bx_sendcounts[%d] = %d\tbx_displs[%d] = %d\n", i, bx_sendcounts[i], i, bx_displs[i]);
    }
  }*/

  float *a_local; 
  float *x_local;  
  float *b_local;  

  int num_pproc = bx_sendcounts[rank];
 
  a_local = (float *) malloc(num_pproc * num * sizeof(float));
  x_local = (float *) malloc(num_pproc * sizeof(float));
  b_local = (float *) malloc(num_pproc * sizeof(float));

  //scatter the matrix A
  MPI_Scatterv(a, a_sendcounts, a_displs, MPI_FLOAT, a_local, 
    num_pproc * num, MPI_FLOAT, 0, MPI_COMM_WORLD);
  //scatter the b
  MPI_Scatterv(b, bx_sendcounts, bx_displs, MPI_FLOAT, b_local, 
    num_pproc, MPI_FLOAT, 0, MPI_COMM_WORLD);
  //scatter the x
  MPI_Scatterv(x, bx_sendcounts, bx_displs, MPI_FLOAT, x_local, 
    num_pproc, MPI_FLOAT, 0, MPI_COMM_WORLD);

  int nit = 0;
  float *old_x = (float *) malloc(num * sizeof(float));
  float *new_x = (float *) malloc(num * sizeof(float));
  MPI_Allgatherv(x_local, num_pproc, MPI_FLOAT, new_x,
      bx_sendcounts, bx_displs, MPI_FLOAT, MPI_COMM_WORLD);

  do
  {
    nit++;
    SWAP(new_x, old_x);
    for(size_t i_local=0; i_local < num_pproc; i_local++)
    {
      int i_global = bx_displs[rank] + i_local;
      x_local[i_local] = b_local[i_local];
      for(size_t j=0; j<num; j++)
      {
        if(j != i_global)
        {
          x_local[i_local] -= a_local[i_local * num + j] * old_x[j];
        }
      }
      x_local[i_local] /= a_local[i_local * num + i_global];
    }

    MPI_Allgatherv(x_local, num_pproc, MPI_FLOAT, new_x,
      bx_sendcounts, bx_displs, MPI_FLOAT, MPI_COMM_WORLD);
    
  } while (curr_max_err(new_x, old_x) > err);

  //put the final solution to x at root process for output
  if(rank==0) 
    SWAP(new_x, x);

  free(a_local);
  free(b_local);
  free(x_local);

  return nit;
}

void clean_up() {
  free(a);
  free(b);
  free(x);
}

/************************************************************/

int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);

  int i;
  int nit = 0; /* number of iterations */
  FILE * fp;
  char output[100] ="";
  
  if( argc != 2)
  {
   printf("Usage: ./gsref filename\n");
   exit(1);
  }

  int comm_size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(rank == 0) 
  {
    double MPI_Wtime(void);
    double start, finish;

    start = MPI_Wtime();

    /* Read the input file and fill the global data structure above */ 
    get_input(argv[1], rank);

    finish = MPI_Wtime();
    //printf("Time taken for reading input = %e seconds\n", finish-start);

    /* Check for convergence condition */
    /* This function will exit the program if the coffeicient will never converge to 
    * the needed absolute error. 
    * This is not expected to happen for this programming assignment.
    */
    check_matrix();
  }

  //nit = serial_solver();
  nit = parallel_solver(comm_size, rank);

  if(rank == 0) 
  {
    /* Writing results to file */
    sprintf(output,"%d.sol",num);
    fp = fopen(output,"w");
    if(!fp)
    {
     printf("Cannot create the file %s\n", output);
     exit(1);
    }
      
    for( i = 0; i < num; i++)
     fprintf(fp,"%f\n",x[i]);

    printf("total number of iterations: %d\n", nit);

    fclose(fp); 
  }

  MPI_Finalize();

  clean_up();

  exit(0);

}
