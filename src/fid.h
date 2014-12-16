#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

/* constants */

#define PROG_NAME      "libfid"
#define PROG_VERSION   "v0.0.0"

#define FID_ALIGNMENT_SSE 16
#define FID_ALIGNMENT_AVX 32


#ifdef __cplusplus
extern "C" {
#endif

/* functions in fid1_matrix.c */

void strale_fid1f_column(int m, int n,
                         float ph,
                         float pb,
                         float pe,
                         float phh,
                         float phb,
                         float phe,
                         float pbh,
                         float pbb,
                         float pbe,
                         float peh,
                         float peb,
                         float pee);

void strale_fid1f_printcol(void);


/* functions in utils.c */

void fatal(const char * format, ...);
void * xmalloc(size_t size, const size_t alignment);


#ifdef __cplusplus
}
#endif
