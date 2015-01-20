#include "fid.h"
#include <x86intrin.h>
#include <stdlib.h>

#define VEC_SIZE 2

static double * dd = NULL;
static long ddlen = 0;

static double phh; 
static double phb; 
static double phe;

static double pbh;
static double pbb;
static double pbe;

static double peh;
static double peb;
static double pee;

/* 
     
     FID 1 model implementation with doubles and no scaling

     Computation is performed on a single column parallel to the vertical
     sequence of the matrix, overwriting the column at every letter of the
     horizontal sequence.  In the end, the result is a single column (the last
     column of the matrix).  

     input:
     m: number of rows 
     n: number of columns
     ph: probability of homology
     pb: probability of birth
     pe: probability of extinction
     phh: probability of homology given homology
     phb: probability of homology given birth
     phe: probability of homology given extinction
     pbh: probability of birth given homology
     pbb: probability of birth given birth
     pbe: probability of birth given extinction
     peh: probability of extinction given homology
     peb: probability of extinction given birth
     pee: probability of extinction given extinction

     output:
     The column (dd) is a linear array of size 3*m, where m is the number of
     rows (length of vertical sequence). Elements dd[0..m-1] ents are the H
     components, dd[m..2*m-1] are the B componenets, and dd[2*m .. 3*m-1] the E
     components.

*/

static void pprint_sse(__m128d x)
{
  double * p = (double *) &x;

  printf("%f ", *p++);
  printf("%f ", *p++);
}

double strale_fid1d_matrix_sse(int m, int n,
                               double ph, double pb, double pe)
{

  double * ee = NULL;    /* pointer to the E section of the previous column */
  double * hh = NULL;    /* pointer to the H section of the previous column */
  double * bb = NULL;    /* pointer to the B section of the previous column */

  double * chh;          /* pointer to the E section of the current column */ 
  double * cbb;          /* pointer to the H section of the current column */
  double * cee;          /* pointer to the B section of the current column */

  __m128d PHH = _mm_set_pd(phh,phh);
  __m128d PHB = _mm_set_pd(phb,phb);
  __m128d PHE = _mm_set_pd(phe,phe);

  __m128d PBH = _mm_set_pd(pbh,pbh);
  __m128d PBB = _mm_set_pd(pbb,pbb);
  __m128d PBE = _mm_set_pd(pbe,pbe);

  __m128d PEH = _mm_set_pd(peh,peh);
  __m128d PEB = _mm_set_pd(peb,peb);
  __m128d PEE = _mm_set_pd(pee,pee);

  __m128d HH;
  __m128d BB;
  __m128d EE;
  __m128d XH;
  __m128d XB;
  __m128d XE;
  __m128d T;
  __m128d T2;
  __m128d T3;
  __m128d T4;

  long i, j;
  unsigned int y = (m + VEC_SIZE-1) & ~(VEC_SIZE-1);
  y += VEC_SIZE;


  /* reallocate matrix if necessary */
  if (3*n*y > ddlen)
  {
    free(dd);
    dd = (double *)xmalloc(3*n*y*sizeof(double), FID_ALIGNMENT_SSE);
    ddlen = 3*n*y;
  }

  /* init pointers to elements */
  chh = hh = dd;
  cbb = bb = hh + y;
  cee = ee = bb + y;

  double b0 = ph*phb + pb*pbb + pe*peb;
  double e0 = ph*phe + pb*pbe + pe*pee;
  double h0 = ph*phh + pb*pbh + pe*peh;

  /* precompute cell (0,1) */
  chh[1] = 0;
  cee[1] = 0;
  cbb[1] = b0;

  /* precompute cell (1,1) */
  chh[2] = h0;
  cbb[2] = peb*e0;
  cee[2] = pbe*b0;

  /* precompute the rest of column 1, i.e. cells (i,1) */
  for (i = 3; i < y; ++i)
  {
    cee[i] = chh[i-1]*phe + cbb[i-1]*pbe + cee[i-1]*pee;
    chh[i] = peh*e0;
    e0 *= pee;
    cbb[i] = peb*e0;
  }

  /* iterate through columns starting from the second */
  for (j = 1; j < n; ++j)
  {
    /* point each entry to the beginning (first element) */
    chh =  ee + y;
    cbb = chh + y;
    cee = cbb + y;

    /* compute element (j,0) */
    chh[1] = 0;
    cee[1] = 0;
    cbb[1] = bb[1]*pbb;

    XH  = _mm_set_pd ( 0x0000000000000000, 0x0000000000000000 );
    XB  = _mm_set_pd (              bb[1], 0x0000000000000000 );
    XE  = _mm_set_pd ( 0x0000000000000000, 0x0000000000000000 );

    /* iterate cells of a column */
    for (i = 2; i < y; i += VEC_SIZE)
    {
      /* the loop vectorizes the following three lines */
      //cbb[i] = hh[i]*phb + bb[i]*pbb + ee[i]*peb;
      //chh[i] = hh[i-1]*phh + bb[i-1]*pbh + ee[i-1]*peh;
      //cee[i] = chh[i-1]*phe + cbb[i-1]*pbe + cee[i-1]*pee;

      HH = _mm_load_pd(hh+i);
      BB = _mm_load_pd(bb+i);
      EE = _mm_load_pd(ee+i);

      /* compute cbb */
      T = _mm_mul_pd(HH,PHB);
      T = _mm_add_pd(T,_mm_mul_pd(BB, PBB));
      T = _mm_add_pd(T,_mm_mul_pd(EE,PEB));
      _mm_store_pd(cbb+i,T);

      /* compute chh */
      XH = _mm_shuffle_pd(XH, HH, 0x01);
      XB = _mm_shuffle_pd(XB, BB, 0x01);
      XE = _mm_shuffle_pd(XE, EE, 0x01);

      T = _mm_mul_pd(XH,PHH);
      T = _mm_add_pd(T, _mm_mul_pd(XB,PBH));
      T = _mm_add_pd(T, _mm_mul_pd(XE,PEH));
      _mm_store_pd(chh+i,T);

      XH = HH;
      XB = BB;
      XE = EE;

      /* attempt to parallelize these two lines */
      //cee[i] = chh[i-1]*phe + cbb[i-1]*pbe + cee[i-1]*pee;
      //cee[i+1] = chh[i]*phe + cbb[i]*pbe + cee[i]*pee;

      T  = _mm_loadu_pd (chh+i-1);
      T  = _mm_mul_pd(T,PHE);
      T2 = _mm_loadu_pd(cbb+i-1);
      T2 = _mm_mul_pd(T2,PBE);
      T3 = _mm_loadu_pd(cee+i-1);
      T3 = _mm_mul_pd(T3,PEE);
      
      T4 = _mm_add_pd(T,T2);
      T = _mm_add_pd(T4,T3);

      T2 = _mm_shuffle_pd(T,T,0x00);
      T2 = _mm_mul_pd(T2,PEE);
      T2 = _mm_add_pd(T2,T4);

      T = _mm_shuffle_pd(T,T2,0x02);
      _mm_store_pd(cee+i,T);
    }

    hh = chh; bb = cbb; ee = cee;
  }
  
  return (chh[m+1]*phh + cbb[m+1]*pbh + cee[m+1]*peh);
}

void strale_fid1d_init_tpm(double lambda, double gamma)
{
  double el = exp(-lambda);
  
  phh = el/gamma/(1+lambda)+1 - 1/gamma;
  phb = lambda/(gamma*(1+lambda));
  phe = (1-el)/gamma/(1+lambda);

  pbh = el/gamma/(1+lambda);
  pbb = 1-1/gamma+lambda/gamma/(1+lambda);
  pbe = phe;

  peh = lambda*el/gamma/(1-el)/(1+lambda);
  peb = (1-el*(1+lambda))/gamma/(1-el)/(1+lambda);
  pee = 1-1/gamma+lambda/(lambda+1)/gamma;
}

void strale_fid1d_dump_tpm(void)
{
  printf ("Transition Probabilities\n");
  printf ("  |     H         B         E\n");        
  printf ("--+------------------------------\n");
  printf ("H | %.7f %.7f %.7f\n", phh, phb, phe); 
  printf ("B | %.7f %.7f %.7f\n", pbh, pbb, pbe); 
  printf ("E | %.7f %.7f %.7f\n\n", peh, peb, pee); 
}

int main(int argc, char * argv[])
{
  double h;

  strale_fid1d_init_tpm(0.208025*0.02, 2);
  strale_fid1d_dump_tpm();

  h = strale_fid1d_matrix_sse(9000,
                              9000,
                              1,
                              0,
                              0);
  
  printf ("-> Result-sse: %f\n", h);

  return (EXIT_SUCCESS);
}

