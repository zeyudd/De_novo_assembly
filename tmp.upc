#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>
#include <upc.h>
#include <upc_io.h>

#include "upc_packingDNAseq.h"
#include "upc_kmer_hash.h"


int main(int argc, char *argv[]){

	static shared [1] int *gptr[THREADS];
	shared int *ptr = upc_alloc(sizeof(int));

	*ptr = MYTHREAD;

	gptr[MYTHREAD] = ptr;

	if(MYTHREAD ==0){
		int i = 0;
		for(; i < THREADS; i++){
			printf("ptr[%d]=%d\n", i, *ptr[i]);
		}
	}

	return 0;
}
