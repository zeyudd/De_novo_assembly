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
size_t hashseq(size_t  hashtable_size, char *seq, int size)
{
   unsigned long hashval;
   hashval = 5381;
   for(int i = 0; i < size; i++) {
      hashval = seq[i] +  (hashval << 5) + hashval;
   }
   
   return hashval % hashtable_size;
}

/* Returns the hash value of a kmer */
size_t hashkmer(size_t hashtable_size, char *seq)
{
    return hashseq(hashtable_size, seq, KMER_PACKED_LENGTH);
}
   

/* Looks up a kmer in the hash table and returns a pointer to that entry */
size_t lookup_kmer(shared [KMER_PACKED_LENGTH] char *kmer_c, shared [1] kmer_t *kmer_i, 
                    shared [LOAD_FACTOR] bucket_t *hashtable, size_t hashlen, const unsigned char *kmer)
{
    char packedKmer[KMER_PACKED_LENGTH];
    packSequence(kmer, (unsigned char*) packedKmer, KMER_LENGTH);
    size_t hashval = hashkmer(hashlen, (char*) packedKmer);

    bucket_t cur_bucket;
    size_t result;
    
    upc_memget(&cur_bucket, hashtable + hashval, sizeof(bucket_t));
    result = cur_bucket.head;
    char packed_kmer_buf[KMER_PACKED_LENGTH];
    for (; result!=hashlen; ) {
        upc_memget(packed_kmer_buf, kmer_c + result * KMER_PACKED_LENGTH, KMER_PACKED_LENGTH);
        if ( memcmp(packedKmer, packed_kmer_buf, KMER_PACKED_LENGTH * sizeof(char)) == 0 ) {
            return result;
        }
        result = kmer_i[result].next;
   }
   return hashlen;
}

/* Adds a kmer and its extensions in the hash table (note that a memory heap should be preallocated. ) */
int add_kmer(shared kmer_t *kmer_i, shared [KMER_PACKED_LENGTH] char *kmer_c, size_t *posInHeap, 
             shared bucket_t *hashtable, size_t hashlen, const unsigned char *kmer_to_add, char left_ext, char right_ext)
{
   /* Pack a k-mer sequence appropriately */
    char packedKmer[KMER_PACKED_LENGTH+1];
    packedKmer[KMER_PACKED_LENGTH] = '\0';
    packSequence(kmer_to_add, (unsigned char*) packedKmer, KMER_LENGTH);

    size_t hashval = hashkmer(hashlen, (char*) packedKmer);
    size_t pos = *posInHeap;
      
    /* Add the contents to the appropriate kmer struct in the heap */ 
    upc_memput(kmer_c + pos * KMER_PACKED_LENGTH, packedKmer, KMER_PACKED_LENGTH * sizeof(char));

    kmer_i[pos].l_ext = left_ext;
    kmer_i[pos].r_ext = right_ext;
    /* Fix the next pointer to point to the appropriate kmer struct */
    kmer_i[pos].next = hashtable[hashval].head;
    /* Fix the head pointer of the appropriate bucket to point to the current kmer */
    hashtable[hashval].head = pos;

    #if 0
    char unpackedKmer[KMER_LENGTH+1];
    unpackedKmer[KMER_LENGTH] = '\0';
    unpackSequence((unsigned char *)&kmer_c[pos * KMER_PACKED_LENGTH], unpackedKmer, KMER_LENGTH);
    printf("THREAD%d: packed kmer %s, pos=%d\n", MYTHREAD, unpackedKmer, pos);
    #endif
    *posInHeap += THREADS; 
    return 0;
}

/* Adds a k-mer in the start list by using the memory heap (the k-mer was "just added" in the memory heap at position posInHeap - 1) */
void addKmerToStartList(shared kmer_t *kmer_i, size_t posInHeap, start_kmer_t **startKmersList)
{
    start_kmer_t *new_entry;   
    size_t prevPosInHeap = posInHeap - THREADS;  
    new_entry = (start_kmer_t*)malloc(sizeof(start_kmer_t));
    new_entry->next = (*startKmersList);
    new_entry->kmerPtr = prevPosInHeap;
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
