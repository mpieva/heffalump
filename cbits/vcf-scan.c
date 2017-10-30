#include "vcf-scan.h"

// Scanner for (our) compressed VCF.  We ignore most of the file,
// hand-code the rest.
//
// The intention is to call from Haskell; that way we get
// multithreading and non-blocking I/O.  A scanner has a large input
// buffer, a small working buffer and a structure to put the output in.
// We scan until the working buffer is empty, then unzip a BGZF block,
// keep scanning.  We return when either the input buffer is
// sufficiently empty or we managed to produce a record.

#define inc() { ++p ; if( p == pe ) return 0 ; }
#define skip() { while( *p != '\t' ) inc() ; inc() ; }

// Scan the header.  Skips over meta info lines and the header.  If it
// completes, it NUL-terminates the header line and returns a pointer to
// the first sample.  Otherwise it may consume input and will return 0.
char *scan_hdr( struct scanner *sc )
{
    char *p = sc->next_work, *pe = sc->last_work ;
    // keep skipping over meta lines.  inc() bails out if we're out of
    // data, else we make progress by (permanently) skipping lines.
    if( p == pe ) return 0 ;
    while(1) {
        inc() ;
        if( p[-1] != '#' || p[0] != '#' ) break ;
        while( *p != '\n' ) inc() ;
        inc() ;
    }

    // next line should be header
    for( int i = 0 ; i != 9 ; ++i ) skip() ;   // skip CHROM..FORMAT
    char* sp = p ;
    while( *p != '\n' ) inc() ;
    inc() ;
    // if we get here, the whole header is available.  We NUL-terminate
    // and return the pointer to the samples.
    p[-1] = 0 ;
    sc->next_work = p;
    return sp ;
}

// Skip over header.  This is special dispensation for ONE SINGLE file
// that has a crooked header :(  Returns 0 if it is out of input, 1 when
// done.
int skip_hdr( struct scanner *sc )
{
    char *p = sc->next_work, *pe = sc->last_work ;
    // keep skipping over header lines.  inc() bails out if we're out of
    // data, else we make progress by (permanently) skipping lines.
    if( p != pe ) {
        while(1) {
            if( p[0] != '#' ) break ;
            while( *p != '\n' ) inc() ;
            inc() ;
        }
    }
    sc->next_work = p;
    return 1 ;
}

// Scan one entry, store in the appropriate fields in the scanner
// structure itself.  If it hits the end of the working buffer, it
// returns 0, pretends it hasn't scanned anything, and must be called
// again.  If it scans an invariant site, it loops without returning
// anything.
int scan_vcf1( struct scanner *sc )
{
    char *p = sc->next_work, *pe = sc->last_work ;
    while( p != pe ) {
        sc->refseq = p ;
        while( *p != '\t' ) inc() ;
        sc->erefseq = p ;
        inc() ;

        // parse position (positive int)
        unsigned pos = 0 ;
        while( *p != '\t' ) {
            pos = 10 * pos + *p - '0' ;
            inc() ;
        }
        sc->pos = pos ;
        inc() ;

        // skip id
        skip() ;
        sc->alleles = p ;

        // skip ref
        skip() ;
        char *thetab = p-1 ; // the tab between ref and alt

        // if alt=='.', we store no alt alleles 
        if( *p == '.' ) {
            skip() ;
            sc->ealleles = p-3 ;
        } else {
            skip() ;
            sc->ealleles = p-1 ;
        }

        // skip qual, filter, info, format, but not the last tab
        skip() ;
        skip() ;
        skip() ;
        while( *p != '\t' ) inc() ;

        // to end of line, read genotypes
        uint8_t *pgts = sc->gts ;
        while( *p != '\n' ) {
            inc() ;
            // scan one diploid GT.  have to allow for more than one
            // digit :(
            if( *p == '.' ) {
                *pgts++ = 0 ;
                inc() ;
            } else {
                uint8_t g = 0 ;
                while (*p >= '0' && *p <= '9') {
                    g = 10*g + *p - '0' ;
                    inc() ; 
                }
                *pgts++ = (g+1) << 1 ;
            }

            if( *p != '|' && *p != '/' ) {
                *pgts++ = -1 ;
            } else {
                uint8_t ph = *p == '|' ? 1 : 0 ;
                inc() ;

                if( *p == '.' ) {
                    *pgts++ = 0 ;
                    inc() ;
                } else {
                    uint8_t g = 0 ;
                    while (*p >= '0' && *p <= '9') {
                        g = 10*g + *p - '0' ;
                        inc() ;
                    }
                    *pgts++ = ((g+1) << 1) | ph ;
                }
            } 

            while( *p != '\t' && *p != '\n') inc() ;
        }
        // don't use the inc() macro here, we don't need more input
        sc->next_work = ++p ;
        *thetab = ',' ; // all alleles are now comma-separated
        return 1 ;
    }
    return 0 ;
}

void free_scanner( struct scanner* sc )
{
    if( sc ) {
        free( sc->working_buffer ) ;
        free( sc->input_buffer ) ;
        free( sc ) ;
    }
}


