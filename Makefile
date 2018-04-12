CC = gcc
UPCC = upcc

KMER_LENGTH 		= 19
KMER_PACKED_LENGTH 	= $(shell echo $$((($(KMER_LENGTH)+3)/4)))

# Add -std=gnu99 to CFLAGS if use gnu compiler
CFLAGS 	= -O3 -std=gnu99 
DEFINE 	= -DKMER_LENGTH=$(KMER_LENGTH) -DKMER_PACKED_LENGTH=$(KMER_PACKED_LENGTH)
HEADERS	= contig_generation.h kmer_hash.h packingDNAseq.h
UPC_HEADERS = upc_contig_generation.h upc_kmer_hash.h upc_packingDNAseq.h
LIBS	=
SCRIPTS = job-serial job-upc job-scale-single-node job-scale-multi-node init.sh
SUBMIT  = $(HEADERS) $(UPC_HEADERS) serial.c pgen.upc sort.cpp $(SCRIPTS) Makefile members.txt report.pdf

TARGETS	= serial pgen

all: 	$(TARGETS)

serial: serial.c $(HEADERS)
		$(CC) $(CFLAGS) -o $@ $< -DKMER_LENGTH=$(KMER_LENGTH) -DKMER_PACKED_LENGTH=$(KMER_PACKED_LENGTH) $(LIBS)

pgen:	pgen.upc $(UPC_HEADERS)
		$(UPCC) $(UPCFLAGS) -Wc,"$(CFLAGS)" -o $@ $< $(DEFINE) $(LIBS)

clean :
	rm -f *.o
	rm -rf $(TARGETS)

submit:
	tar -czvf ding_yuxin_guanhong_hw3.tar.gz $(SUBMIT)
