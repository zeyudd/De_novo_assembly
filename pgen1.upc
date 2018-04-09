#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>
#include <upc.h>
#include <upc_io.h>

//#include "upc_packingDNAseq.h"
//#include "upc_kmer_hash.h"

//upc_lock_t heap_lock;
typedef struct data_t data_t;
struct data_t{
	char mer[2];
	int pos;
};


int main(int argc, char *argv[]){

	shared [1] data_t * data;
	upc_barrier;

	data = (shared [1] data_t * )upc_all_alloc(THREADS, sizeof(data_t));


	data[MYTHREAD].pos = MYTHREAD;
	data[MYTHREAD].mer[0] = 'A' + MYTHREAD;
	data[MYTHREAD].mer[1] = '\0';

	printf("data[%d] = (%d, %s)\n", MYTHREAD, data[MYTHREAD].pos, data[MYTHREAD].mer);

	upc_barrier;

#if 0
	if(MYTHREAD == 0){
		int i;
		for(i = 0; i < THREADS; i++){
			printf("data[%d] = (%d, %s)\n", i, data[i].pos, data[i].mer);
		}
	}
#endif
	return 0;
}
