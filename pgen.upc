#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>
#include <upc.h>
#include <upc_io.h>

#include "packingDNAseq.h"
#include "kmer_hash.h"


int main(int argc, char *argv[]){

	/** Declarations **/
	double inputTime=0.0, constrTime=0.0, traversalTime=0.0;

	/** Read input **/
	upc_barrier;
	inputTime -= gettime();
	///////////////////////////////////////////
	//   code for input file reading : begin //
	///////////////////////////////////////////
	char * input_UFX_name = argv[1];
	int64_t nKmers = getNumKmersInUFX(input_UFX_name);
	
	int64_t avg_nKmers = nKmers / THREADS;
	int64_t rem_nKmers = nKmers % THREADS;
	
	int64_t my_lines_to_read = (MYTHREAD < rem_nKmers)? (avg_nKmers+1) : avg_nKmers;
	int64_t my_lines_to_skip = (MYTHREAD <= rem_nKmers)? (avg_nKmers * MYTHREAD + MYTHREAD) : (avg_nKmers * MYTHREAD + rem_nKmers);

	int64_t my_read_size = my_lines_to_read * LINE_SIZE;
	int64_t my_read_offset = my_lines_to_skip * LINE_SIZE;

	unsigned char* my_buffer = 
	  (unsigned char*) malloc((my_read_size) * sizeof(unsigned char)); // local buffer
	
	upc_file_t *input_file;
	input_file = upc_all_fopen(input_UFX_name, UPC_RDONLY | UPC_INDIVIDUAL_FP, 0, NULL);
	upc_all_fseek(input_file, my_read_offset*sizeof(unsigned char), UPC_SEEK_SET);
	int64_t cur_chars_read = upc_all_fread_local(input_file, my_buffer, sizeof(unsigned char), my_read_size,
	                                             UPC_IN_ALLSYNC | UPC_OUT_ALLSYNC);
	upc_all_fclose(input_file);
	
	printf("Reading Finished on Thread#%d of %d threads.\n", MYTHREAD, THREADS);
	

	///////////////////////////////////////////
	//   code for input file reading : end   //
	///////////////////////////////////////////
	upc_barrier;
	inputTime += gettime();

	/** Graph construction **/
	constrTime -= gettime();
	///////////////////////////////////////////
	// Your code for graph construction here //
	///////////////////////////////////////////
	upc_barrier;
	constrTime += gettime();

	/** Graph traversal **/
	traversalTime -= gettime();
	////////////////////////////////////////////////////////////
	// Your code for graph traversal and output printing here //
	// Save your output to "pgen.out"                         //
	////////////////////////////////////////////////////////////
	upc_barrier;
	traversalTime += gettime();

	/** Print timing and output info **/
	/***** DO NOT CHANGE THIS PART ****/
	if(MYTHREAD==0){
		printf("%s: Input set: %s\n", argv[0], argv[1]);
		printf("Number of UPC threads: %d\n", THREADS);
		printf("Input reading time: %f seconds\n", inputTime);
		printf("Graph construction time: %f seconds\n", constrTime);
		printf("Graph traversal time: %f seconds\n", traversalTime);
	}
	return 0;
}
