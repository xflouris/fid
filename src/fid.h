#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>

/* constants */

#define PROG_NAME      "libfid"
#define PROG_VERSION   "v0.0.0"

#define FID_ALIGNMENT_SSE 16
#define FID_ALIGNMENT_AVX 32


#ifdef __cplusplus
extern "C" {
#endif

/* functions in fid1_matrix.c */

double strale_fid1d_column(int m, int n,
                           double ph,
                           double pb,
                           double pe);

void strale_fid1d_printcol(void);

void strale_fid1d_dump_tpm(void);

/* functions in utils.c */

void fatal(const char * format, ...);
void * xmalloc(size_t size, const size_t alignment);


#ifdef __cplusplus
}
#endif
