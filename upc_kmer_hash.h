#ifndef UPC_KMER_HASH_H
#define UPC_KMER_HASH_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h> 
#include <string.h>
#include "upc_contig_generation.h"
#include <upc.h>

/* Auxiliary function for computing hash values */
int64_t hashseq(int64_t  hashtable_size, char *seq, int size)
{
   unsigned long hashval;
   hashval = 5381;
   for(int i = 0; i < size; i++) {
      hashval = seq[i] +  (hashval << 5) + hashval;
   }
   
   return hashval % hashtable_size;
}

/* Returns the hash value of a kmer */
int64_t hashkmer(int64_t  hashtable_size, char *seq)
{
   return hashseq(hashtable_size, seq, KMER_PACKED_LENGTH);
}

/* Looks up a kmer in the hash table and returns a pointer to that entry */
shared kmer_t* lookup_kmer(shared bucket_t *hashtable, int64_t hashtable_size, const unsigned char *kmer)
{
    char packedKmer[KMER_PACKED_LENGTH];
    packSequence(kmer, (unsigned char*) packedKmer, KMER_LENGTH);
    int64_t hashval = hashkmer(hashtable_size, (char*) packedKmer);
    bucket_t cur_bucket;
    shared kmer_t *result;
   
    cur_bucket = hashtable[hashval];
    result = cur_bucket.head;
   
    for (; result!=NULL; ) {
        if ( memcmp(packedKmer, (unsigned char*)result->kmer, KMER_PACKED_LENGTH * sizeof(char)) == 0 ) {
            return result;
        }
        result = result->next;
   }
   return NULL;
}

/* Adds a kmer and its extensions in the hash table (note that a memory heap should be preallocated. ) */
int add_kmer(shared bucket_t *hashtable, int64_t hashtable_size, shared kmer_t *heap, int64_t *posInHeap, 
            shared [KMER_PACKED_LENGTH] char * kmers, const unsigned char *kmer, char left_ext, char right_ext)
{
   /* Pack a k-mer sequence appropriately */
   char packedKmer[KMER_PACKED_LENGTH];
   packSequence(kmer, (unsigned char*) packedKmer, KMER_LENGTH);
   int64_t hashval = hashkmer(hashtable_size, (char*) packedKmer);

   int64_t pos = *posInHeap;
      
   /* Add the contents to the appropriate kmer struct in the heap */
   upc_memput((kmers[pos * KMER_PACKED_SIZE], packedKmer, KMER_PACKED_LENGTH * sizeof(char));
   heap[pos].kmer = &kmers[pos * KMER_PACKED_SIZE];
   heap[pos].l_ext = left_ext;
   heap[pos].r_ext = right_ext;
   
   /* Fix the next pointer to point to the appropriate kmer struct */
   heap[pos].next = hashtable[hashval].head;
   /* Fix the head pointer of the appropriate bucket to point to the current kmer */
   hashtable[hashval].head = &(heap[pos]);
   
   *posInHeap += THREADS; 
   return 0;
}

/* Adds a k-mer in the start list by using the memory heap (the k-mer was "just added" in the memory heap at position posInHeap - 1) */
void addKmerToStartList(shared kmer_t *heap, int64_t posInHeap, start_kmer_t **startKmersList)
{
    start_kmer_t *new_entry;
    shared kmer_t *ptrToKmer;
   
    int64_t prevPosInHeap = posInHeap - THREADS;
    ptrToKmer = &(heap[prevPosInHeap]);
    new_entry = (start_kmer_t*)malloc(sizeof(start_kmer_t));
    new_entry->next = (*startKmersList);
    new_entry->kmerPtr = ptrToKmer;
    (*startKmersList) = new_entry;
}

/* Deallocation functions */
int dealloc_heap(memory_heap_t *memory_heap)
{
   upc_free(memory_heap->heap);
   return 0;
}

int dealloc_hashtable(hash_table_t *hashtable)
{
   upc_free(hashtable->table);
   return 0;
}


#endif // UPC_KMER_HASH_H
