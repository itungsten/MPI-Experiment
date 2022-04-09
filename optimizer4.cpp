#pragma GCC optimize(3,"Ofast","inline")
#include "mpi.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#define BLOCK_LOW(id,p,n)  ((1LL*id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1)
int aux[16]={0,0,1,2,2,3,3,4,4,4,4,5,5,6,6,6};
int main(int argc, char *argv[])
{
   double elapsed_time; /* Parallel execution time */
   int    count;        /* Local prime count */
   int    global_count; /* Global prime count */
   int    high_block;   /* Highest value on this proc */
   int    low_block;    /* Lowest value on this proc */
   int    first;        /* Index of first multiple */
   int    id;           /* Process ID number */
   int    index;        /* Index of current prime */
   bool  *marked;       /* Portion of 2,...,'n' */
   bool  *marked2;       /* Portion of 2,...,'n' */
   int    n;            /* Sieving from 2, ..., 'n' */
   int    p;            /* Number of processes */
   int    prime;        /* Current prime */
   int    num_block;         /* Elements in 'marked' */
   int    no;

   MPI_Init (&argc, &argv);
   MPI_Comm_rank (MPI_COMM_WORLD, &id);
   MPI_Comm_size (MPI_COMM_WORLD, &p);
   MPI_Barrier(MPI_COMM_WORLD);
   /* Start the timer */
   elapsed_time = -MPI_Wtime();

   n = atoi(argv[1]);
   int chunk = 128*1024;
   if(argc==3)chunk=atoi(argv[2]);

   if(n<=15){
      if(!id){
         elapsed_time += MPI_Wtime();
         printf ("There are %d primes less than or equal to %d\n",aux[n], n);
         printf ("SIEVE (%d) %10.6f\n", p, elapsed_time);   
      }
      MPI_Finalize ();
      return 0;  
   }
	
	int sqr=sqrt(n),sqrr=sqrt(sqr);
	int primes_size = (sqr - 3)/2 + 1; // number of adds in 3 .. sqr 
	bool* primes = (bool *) malloc(primes_size);
	memset(primes,0,primes_size);

	index = 0; prime = 3;
	do {
		for (int i = (prime*prime-3)/2; i < primes_size; i += prime) primes[i] = 1;
		while (primes[++index]);
		prime = (index<<1)+3;
	} while (prime <= sqrr);
   std::vector<int> primeVec; primeVec.reserve(primes_size);
   index = 1; prime = 5;
   do{
      primeVec.push_back(prime);
      while (primes[++index]);
      prime = (index<<1)+3;
   }while(prime <= sqr);

   no = (n+1)/6; // nubmer of blocks from 5,7... 

	/* Figure out this process's share of the array, as
      well as the integers represented by the first and
      last array elements */
   low_block = BLOCK_LOW(id, p, no);
   high_block = BLOCK_HIGH(id, p, no);
   num_block = high_block-low_block+1;

	/* Allocate this process' share of the array */
	marked = (bool *) malloc(num_block<<1);
   marked2 = marked+num_block;
   memset(marked,0,num_block<<1);

	for (int sec = 0; sec < num_block; sec += chunk) {
      for(auto prime : primeVec){
         first = (prime*prime-5)/6-low_block-sec;
         if(first<0) first=(first%prime+prime)%prime;
         for (int i = first+sec; i < first+sec+chunk && i < num_block; i += prime)marked[i] = 1;
      }
	}
   for (int sec = 0; sec < num_block; sec += chunk) {
      for(auto prime : primeVec){
         int offset=1<<(prime%6==1)+1;
         first = (prime*(prime+offset)-5)/6-low_block-sec;
         if(first<0) first=(first%prime+prime)%prime;
         for (int i = first+sec; i < first+sec+chunk && i < num_block; i += prime)marked2[i] = 1;
      }
	}
	count = 0;
	for (int i=0; i<num_block<<1; i++)if (!marked[i]) count++;
	MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	/* Stop the timer */
	elapsed_time += MPI_Wtime();

	/* Print the results */
   if (!id) {
      printf ("There are %d primes less than or equal to %d\n",global_count+2, n);
      printf ("SIEVE (%d) %10.6f\n", p, elapsed_time);
   }
   free(marked);
   free(primes);
   MPI_Finalize ();
   return 0;
}