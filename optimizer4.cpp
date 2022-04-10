#pragma GCC optimize(3,"Ofast","inline")
#include "mpi.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#define BLOCK_LOW(id,p,n)  ((1LL*id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1)
#define MIN(a,b) ((a)<(b)?(a):(b))
int aux[8]={0,0,1,2,2,3,3,4};
int begArr[8]={31,7,11,13,17,19,23,29};
int posArr[30]={0,0,0,0,0,0,0,1,0,0,0,2,0,3,0,0,0,4,0,5,0,0,0,6,0,0,0,0,0,7};
int offsetArr[8][8]={
{0,6,0,24,6,0,24,0},
{6,24,6,6,24,24,6,24},
{10,16,20,4,26,10,14,20},
{12,12,12,18,12,18,18,18},
{16,4,26,16,14,4,26,14},
{18,0,18,0,0,12,0,12},
{22,22,2,28,2,28,8,8},
{28,10,8,10,20,22,20,2}
};
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

   if(n<=7){
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
   index = 2; prime = 7;
   do{
      primeVec.push_back(prime);
      while (primes[++index]);
      prime = (index<<1)+3;
   }while(prime <= sqr);

   no = (n+23)/30; // nubmer of blocks from 7... 
   low_block = BLOCK_LOW(id, p, no);
   high_block = BLOCK_HIGH(id, p, no);
   num_block = high_block-low_block+1;
	marked = (bool *) malloc(MIN(num_block,chunk));

   count=0;
   for (int sec = 0; sec < num_block; sec += chunk) {
      for(int k=0;k<8;++k){
         int block=MIN(chunk,num_block-sec);
         if(n<begArr[k]+(low_block+sec+block-1)*30)block--;
         memset(marked,0,block);
         for(auto prime : primeVec){
            int pos=posArr[prime%30];
            int offset=offsetArr[k][pos];
            first = (prime*(prime+offset)-7)/30-low_block-sec;
            if(first<0) first=(first%prime+prime)%prime;
            for (int i = first; i < block; i += prime)marked[i] = 1;
         }
         for (int i=0; i<block; i++)if (!marked[i]) count++;
      }
   }
	
	MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	/* Stop the timer */
	elapsed_time += MPI_Wtime();

	/* Print the results */
   if (!id) {
      printf ("There are %d primes less than or equal to %d\n",global_count+3, n);
      printf ("SIEVE (%d) %10.6f\n", p, elapsed_time);
   }
   free(marked);
   free(primes);
   MPI_Finalize ();
   return 0;
}