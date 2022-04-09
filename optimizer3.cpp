#include "mpi.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#define BLOCK_LOW(id,p,n)  ((1LL*id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1)

int main(int argc, char *argv[])
{
   double elapsed_time; /* Parallel execution time */
   int    count;        /* Local prime count */
   int    global_count; /* Global prime count */
   int    high_value;   /* Highest value on this proc */
   int    low_value;    /* Lowest value on this proc */
   int    first;        /* Index of first multiple */
   int    id;           /* Process ID number */
   int    index;        /* Index of current prime */
   bool  *marked;       /* Portion of 2,...,'n' */
   int    n;            /* Sieving from 2, ..., 'n' */
   int    p;            /* Number of processes */
   int    prime;        /* Current prime */
   int    size;         /* Elements in 'marked' */
   int    no;

   MPI_Init (&argc, &argv);
   MPI_Comm_rank (MPI_COMM_WORLD, &id);
   MPI_Comm_size (MPI_COMM_WORLD, &p);
   MPI_Barrier(MPI_COMM_WORLD);
   /* Start the timer */
   elapsed_time = -MPI_Wtime();

   n = atoi(argv[1]);
   int chunk = 131072;
   if(argc==3)chunk=atoi(argv[2]);
   no = (n-1)/2; // nubmer of odds in 3 to n

	/* Figure out this process's share of the array, as
      well as the integers represented by the first and
      last array elements */
   low_value = 2 * BLOCK_LOW(id, p, no) + 3;
   high_value = 2 * BLOCK_HIGH(id, p, no) + 3;
   size = (high_value-low_value)/2+1;

	/* Allocate this process' share of the array */
	marked = (bool *) malloc(size);
	for (int i=0; i<size; i++) marked[i] = 0;
	
	int sqr=sqrt(n),sqrr=sqrt(sqr);
	int primes_size = (sqr - 3)/2 + 1; // number of adds in 3 .. sqr 
	bool* primes = (bool *) malloc(primes_size);
	for (int i=0; i<primes_size; i++) primes[i] = 0;

	index = 0;
	prime = 3;
	do {
		for (int i = (prime*prime-3)/2; i < primes_size; i += prime) primes[i] = 1;
		while (primes[++index]);
		prime = (index<<1)+3;
	} while (prime <= sqrr);
	
	for (int sec = 0; sec < size; sec += chunk) {
		index = 0;
		prime = 3;
		int lv = 2*((low_value-3)/2+sec)+3;
		do {
			first = (prime*prime-lv)/2;
			if(first<0) first=(first%prime+prime)%prime;
			for (int i = first+sec; i < first+sec+chunk && i < size; i += prime)
				marked[i] = 1;
			while (primes[++index]);
			prime = 2*index + 3;
		} while (prime*prime <= n);
	}
	
	count = 0;
	for (int i=0; i<size; i++)if (!marked[i]) count++;
	MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	/* Stop the timer */
	elapsed_time += MPI_Wtime();

	/* Print the results */
   if (!id) {
      printf ("There are %d primes less than or equal to %d\n",global_count+1, n);
      printf ("SIEVE (%d) %10.6f\n", p, elapsed_time);
   }
   free(marked);
   free(primes);
   MPI_Finalize ();
   return 0;
}