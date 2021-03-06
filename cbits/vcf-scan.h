#include <stdint.h>
#include <stdlib.h>

typedef struct scanner {
    // working buffer
    char *working_buffer ;
    char *next_work ;
    char *last_work ;

    // output structure
    char *refseq ;          // points to chrom
    char *erefseq ;         // points to chrom
    char *alleles ;         // points to ref, must be copied out immediately
    char *ealleles ;
    uint32_t pos ;
    uint8_t gts[2] ;        // genotypes (two bytes per sample)
} scanner ;

char *scan_hdr( struct scanner *sc ) ;
int   skip_hdr( struct scanner *sc ) ;
int   scan_vcf1( struct scanner *sc ) ;

