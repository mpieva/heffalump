#include <htslib/hts.h>
#include <htslib/sam.h>

struct plp_aux_t ;

struct plp_aux_t *pileup_init( int mmapq, const char *fn ) ;
int pileup_step( struct plp_aux_t *data, int *tid, int *pos, int vsize, uint32_t *vec ) ;
void pileup_done( struct plp_aux_t *data ) ;
