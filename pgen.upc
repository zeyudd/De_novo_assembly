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

	/** Declarations **/
	double inputTime=0.0, constrTime=0.0, traversalTime=0.0;

	/** Read input **/
	upc_barrier;
	inputTime -= gettime();

	///////////////////////////////////////////
	//   code for input file reading : begin //
	
	char * input_UFX_name = argv[1];
	int64_t nKmers = getNumKmersInUFX(input_UFX_name);
	
	int64_t avg_nKmers = nKmers / THREADS;
	int64_t rem_nKmers = nKmers % THREADS;
	
	int64_t my_lines_to_read = (MYTHREAD < rem_nKmers)? (avg_nKmers+1) : avg_nKmers;
	int64_t my_lines_to_skip = (MYTHREAD <= rem_nKmers)? (avg_nKmers * MYTHREAD + MYTHREAD) : (avg_nKmers * MYTHREAD + rem_nKmers);

	int64_t my_read_size = my_lines_to_read * LINE_SIZE;
	int64_t my_read_offset = my_lines_to_skip * LINE_SIZE;

	unsigned char* my_buffer = 
	  (unsigned char*) malloc((my_read_size) * sizeof(unsigned char));
	
	upc_file_t *input_file;
	input_file = upc_all_fopen(input_UFX_name, UPC_RDONLY | UPC_INDIVIDUAL_FP, 0, NULL);
	upc_all_fseek(input_file, my_read_offset*sizeof(unsigned char), UPC_SEEK_SET);
	int64_t cur_chars_read = upc_all_fread_local(input_file, my_buffer, sizeof(unsigned char), my_read_size,
	                                             UPC_IN_ALLSYNC | UPC_OUT_ALLSYNC);
	upc_all_fclose(input_file);
	
	printf("Thread#%d of %d: read %d kMers, skip %d kMers.\n", 
	  MYTHREAD, THREADS, my_lines_to_read, my_lines_to_skip);
	
	
	//   code for input file reading : end   //
	///////////////////////////////////////////
	upc_barrier;
	inputTime += gettime();

	/** Graph construction **/
	constrTime -= gettime();
	///////////////////////////////////////////
	//  code for graph construction: begin   //
	
    /* Initialize lookup table that will be used for the DNA packing routines */
    init_LookupTable();


	/* Creates a hash table and (pre)allocates memory for the memory heap */
/*
shared hash_table_t* upc_create_hash_table(int64_t nEntries, shared memory_heap_t *memory_heap)
{
   shared hash_table_t *result;
   int64_t n_buckets = nEntries * LOAD_FACTOR;

   result = (shared hash_table_t*) upc_all_alloc(sizeof(hash_table_t), 1);
   result->size = n_buckets;
   result->table = (shared bucket_t*) upc_all_alloc(n_buckets , sizeof(bucket_t));
   
   if (result->table == NULL) {
      fprintf(stderr, "ERROR: Could not allocate memory for the hash table: %lld buckets of %lu bytes\n", n_buckets, sizeof(bucket_t));
      exit(1);
   }
   
   memory_heap->heap = (shared kmer_t *) upc_all_alloc(nEntries, sizeof(kmer_t));
   if (memory_heap->heap == NULL) {
      fprintf(stderr, "ERROR: Could not allocate memory for the heap!\n");
      exit(1);
   }
   memory_heap->posInHeap = 0;
   
   return result;
}
*/	
	shared hash_table_t *hashtable;
    shared memory_heap_t memory_heap;

	hashtable = upc_create_hash_table(nKmers, &memory_heap);










#if 0
	int64_t n = lines_to_read; // total kmers processed by the current thread
	// generate hash table and heap
	shared int64_t* next_index = (shared int64_t*)upc_all_alloc(nKmers, sizeof(int64_t));
    shared kmer_t* memory_heap = (shared kmer_t*)upc_all_alloc(nKmers, sizeof(kmer_t)); 
  
    // initialize hash_table
   int64_t tablesize = nKmers * LOAD_FACTOR;
   shared int64_t* hash_table = (shared int64_t*)upc_all_alloc(tablesize, sizeof(int64_t));
   int64_t i;
   upc_forall(i=0; i<tablesize; i++; &hash_table[i]) // initial hash table
    hash_table[i] = -1;
   upc_forall(i=0; i<nKmers; i++; &next_index[i]) // initial next index
    next_index[i] = -1;
  
   upc_barrier; // need to synchronize before we actually start!
  
   int64_t k = lines_to_ignore; // global kmer index
   int64_t ptr = 0;
  
  
  
	while (ptr < cur_chars_read) {
	
    char left_ext = (char) buffer[ptr+KMER_LENGTH+1];
    char right_ext = (char) buffer[ptr+KMER_LENGTH+2];

    /* Add k-mer to hash table */
    add_kmer(next_index, k, hash_table, tablesize, memory_heap, &buffer[ptr], left_ext, right_ext);
    
    /*
    //// test whether it is inserted to the hash_table
    {
      char packedKmer[KMER_PACKED_LENGTH];
      packSequence(&buffer[ptr], (unsigned char*) packedKmer, KMER_LENGTH);
      int64_t hashval = hashkmer(tablesize, (char*) packedKmer);
      int64_t p = hash_table[hashval];
      while(p != k && p != -1) {
        p = next_index[p];
      }
      if(p != k) {
        printf("<ADD Failure!> Cannot find the kmer index [%d] in the hash_table on Thread#%d!\n", k, MYTHREAD);
      }
      
      if(ok==1) {
      
      kmer_t tmp_kmer;
      upc_memget(& tmp_kmer, &memory_heap[k], sizeof(kmer_t));
      
      
      char unpackedKmer[KMER_LENGTH+1];
      unpackSequence((unsigned char*) tmp_kmer.kmer,  (unsigned char*) unpackedKmer, KMER_LENGTH);
      if(memcmp(unpackedKmer, &buffer[ptr], KMER_LENGTH * sizeof(char)) != 0) {
        printf("<Fetch ERROR> the fetched kmer cannot be correctly recovered on Thread#%d!\n", MYTHREAD);
      }

      }
    }
    
    */
    /* Move to the next k-mer in the input buffer */
    ptr += LINE_SIZE;
    k ++; 
    
    //i++;
  }
#endif


	
	//  code for graph construction: end     //
	///////////////////////////////////////////
	upc_barrier;
	constrTime += gettime();

	/** Graph traversal **/
	traversalTime -= gettime();
	////////////////////////////////////////////////////////////
	// code for graph traversal and output printing: begin    //
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
