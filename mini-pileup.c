/* Pileup using samtools' code.  Example of minimum ugliness, I think. */

#include "mini-pileup.h"

#include <stdlib.h>
#include <stdio.h>

struct plp_aux_t {
    htsFile *fp ;
    bam_hdr_t *h ;
    bam_plp_t iter ;
    int min_mapq ;
    int ignore_indels ;
} ;

static int has_indels( bam1_t *b ) 
{
    for( int i = 0 ; i != b->core.n_cigar ; ++i ) {
        switch( bam_cigar_op( bam_get_cigar(b)[i] ) ) {
            case BAM_CINS:
            case BAM_CDEL:
            case BAM_CSOFT_CLIP:
                return 1 ;
        }
    }
    return 0 ;
}


// Fetches stuff, puts it into *b, returns 1 if something was fetched.
// Gets called by bam_plp_auto to fetch fresh input.  We just read from
// a single bam file.
static int plp_func( void *data, bam1_t *b )
{
    struct plp_aux_t *ma = (struct plp_aux_t*)data ;
    int ret ;
    do {
        ret = sam_read1(ma->fp, ma->h, b);
    } while( ret >= 0 && (
                b->core.tid < 0 || 
                (b->core.flag&BAM_FUNMAP) ||
                b->core.qual < ma->min_mapq ||
                (ma->ignore_indels && has_indels(b))) ) ;
    return ret;
}

struct plp_aux_t *pileup_init( int min_mapq_, int ignore_indels_, const char *fn )
{
    struct plp_aux_t *data = malloc( sizeof( struct plp_aux_t ) ) ;
    if( data ) {
        data->fp = hts_open( fn, "rb" ) ;
        if( data->fp ) {
            // XXX  This header is never freed
            data->h = sam_hdr_read( data->fp ) ;
            data->iter = bam_plp_init( plp_func, data ) ;
            data->min_mapq = min_mapq_ ;
            data->ignore_indels = ignore_indels_ ;
            return data ;
        }
        free( data ) ;
    }
    return 0 ;
}

// All right, fuck this---there doesn't seem to be a clean way to get
// the reference sequence; we get it from a FastA file instead.  ≈:-/

// An array of 32 bit words is supplied.  For each base, we store one
// integer:  bits 0..6:  quality
//           bits 7..7:  reversed?
//           bits 8..11: base
//           bits 12..21: position in the read
//           bits 22..31: length of the read

int pileup_step( struct plp_aux_t *data, int *tid, int *pos, int vsize, uint32_t *vec )
{
    int n_plp ;
    const bam_pileup1_t *plp = bam_plp_auto(data->iter, tid, pos, &n_plp) ;
    if( plp ) {
        for (int j = 0; j < n_plp && j < vsize; ++j) {
            int bqual = plp[j].is_del || plp[j].qpos >= plp[j].b->core.l_qseq
                      ? ( 15 << 8 )
                      : ( bam_seqi(bam_get_seq(plp[j].b), plp[j].qpos) << 8
                        | ( bam_get_qual(plp[j].b)[ plp[j].qpos ] & 0x7f )
                        | ( bam_is_rev(plp[j].b) ? 0x80 : 0x00 ) ) ;

            vec[j] = ( plp[j].b->core.l_qseq < 1024 ? plp[j].b->core.l_qseq : 1023 ) << 22
                   | ( plp[j].qpos < 1024 ? plp[j].qpos : 1023 ) << 12
                   | bqual ;
        }
        return n_plp < vsize ? n_plp : vsize ;
    }
    return -1 ;
}

void pileup_done( struct plp_aux_t *data )
{
    if( data ) {
        if( data->iter ) bam_plp_destroy(data->iter);
        if( data->fp ) hts_close(data->fp) ;
    }
    free( data ) ;
}

