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
	char *mer;
	int pos;
};


int main(int argc, char *argv[]){

	shared [1] data_t * data;

	shared [2] char * mers;

	upc_barrier;

	data = (shared [1] data_t * )upc_all_alloc(THREADS, sizeof(data_t));

	mers = (shared [2] char *)upc_all_alloc(THREADS, 2*sizeof(char)); 

	if(data == NULL){
		printf("upc_all_alloc failed!\n");
		return 0;

	}
	upc_barrier;

	data[MYTHREAD].pos = MYTHREAD;
	data[MYTHREAD].mer = &mers[2*MYTHREAD];
	data[MYTHREAD].mer[0] = 'A' + MYTHREAD;
	data[MYTHREAD].mer[1] = 'D' + MYTHREAD;

	printf("Thread%d: data[%d] = (%d, %s)\n",MYTHREAD, MYTHREAD, data[MYTHREAD].pos, data[MYTHREAD].mer);

	upc_barrier;
	sleep(1);

#if 1
	if(MYTHREAD == 0){
		int i;
		for(i = 0; i < THREADS; i++){
			//printf("data[%d] = (%d, %s)\n", i, data[i].pos, data[i].mer);
			//printf("data[%d] = (%d, %s)\n", i, data[i].pos, &mers[i*2]);
			printf("Thread0: data[%d] = (%d, %c, %c)\n", i, data[i].pos, mers[2*i], mers[2*i+1]);
		}
	}
#endif
	return 0;
}
