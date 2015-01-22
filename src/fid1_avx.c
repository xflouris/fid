#include "fid.h"
#include <x86intrin.h>
#include <stdlib.h>

#define VEC_SIZE 4

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

static void pprint_avx(const char * s, __m256d x)
{
  double * p = (double *) &x;

  printf ("%s",s);

  printf("%f ", *p++);
  printf("%f ", *p++);
  printf("%f ", *p++);
  printf("%f ", *p++);
  printf("\n");
}

double strale_fid1d_matrix_avx(int m, int n,
                               double ph, double pb, double pe)
{

  double * ee = NULL;    /* pointer to the E section of the previous column */
  double * hh = NULL;    /* pointer to the H section of the previous column */
  double * bb = NULL;    /* pointer to the B section of the previous column */

  double * chh;          /* pointer to the E section of the current column */ 
  double * cbb;          /* pointer to the H section of the current column */
  double * cee;          /* pointer to the B section of the current column */

  /* vectors of probabilities */
  __m256d PHH = _mm256_set_pd(phh, phh, phh, phh);
  __m256d PHB = _mm256_set_pd(phb, phb, phb, phb);
  __m256d PHE = _mm256_set_pd(phe, phe, phe, phe);

  __m256d PBH = _mm256_set_pd(pbh, pbh, pbh, pbh);
  __m256d PBB = _mm256_set_pd(pbb, pbb, pbb, pbb);
  __m256d PBE = _mm256_set_pd(pbe, pbe, pbe, pbe);

  __m256d PEH = _mm256_set_pd(peh, peh, peh, peh);
  __m256d PEB = _mm256_set_pd(peb, peb, peb, peb);
  __m256d PEE = _mm256_set_pd(pee, pee, pee, pee);

  __m256d HH;
  __m256d BB;
  __m256d EE;
  __m256d XH;
  __m256d XB;
  __m256d XE;
  __m256d T;
  __m256d T1;
  __m256d T2;
  __m256d T3;
  __m256d T4;

  __m256d XBB;
  __m256d XHH;
  __m256d XEE;

  long i, j;
  unsigned int y = (m + VEC_SIZE-1) & ~(VEC_SIZE-1);
  y += VEC_SIZE;


  /* reallocate matrix if necessary */
  if (3*n*y > ddlen)
  {
    free(dd);
    dd = (double *)xmalloc(3*n*y*sizeof(double), FID_ALIGNMENT_AVX);
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
  chh[3] = 0;
  cee[3] = 0;
  cbb[3] = b0;

  /* precompute cell (1,1) */
  chh[4] = h0;
  cbb[4] = peb*e0;
  cee[4] = pbe*b0;

  /* precompute the rest of column 1, i.e. cells (i,1) */
  for (i = 5; i < y; ++i)
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
    chh[3] = 0;
    cee[3] = 0;
    cbb[3] = bb[3]*pbb;

    /* set top cell of previous column */
    XH  = _mm256_set_pd ( 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000 );
    XB  = _mm256_set_pd (              bb[3], 0x0000000000000000, 0x0000000000000000, 0x0000000000000000 );
    XE  = _mm256_set_pd ( 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000 );

    /* set top cell of current column */
    XHH = _mm256_set_pd ( 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000 );
    XBB = _mm256_set_pd (             cbb[3], 0x0000000000000000, 0x0000000000000000, 0x0000000000000000 );
    XEE = _mm256_set_pd ( 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000 );

    /* iterate cells of a column */
    for (i = 4; i < y; i += VEC_SIZE)
    {
      /* the loop vectorizes the following three lines */
      //cbb[i] = hh[i]*phb + bb[i]*pbb + ee[i]*peb;
      //chh[i] = hh[i-1]*phh + bb[i-1]*pbh + ee[i-1]*peh;
      //cee[i] = chh[i-1]*phe + cbb[i-1]*pbe + cee[i-1]*pee;

      HH = _mm256_load_pd(hh+i);
      BB = _mm256_load_pd(bb+i);
      EE = _mm256_load_pd(ee+i);

      /* compute cbb */
      T = _mm256_mul_pd(HH,PHB);
#ifdef HAVE_FMA
      T = _mm256_fmadd_pd(BB,PBB,T);
      T = _mm256_fmadd_pd(EE,PEB,T);
#else
      T = _mm256_add_pd(T,_mm256_mul_pd(BB, PBB));
      T = _mm256_add_pd(T,_mm256_mul_pd(EE,PEB));
#endif
      _mm256_store_pd(cbb+i,T);

      /* store T2 to compute cee later (B component) */
      T2 = _mm256_permute2f128_pd(T, XBB, 0x03);
      T3 = _mm256_permute_pd(T, 0x00);
      T2 = _mm256_unpackhi_pd(T2, T3);
      XBB = T;
     
      /**** compute chh ****/

      /* concatenate XH and HH and rotate right by one */
      XH = _mm256_permute2f128_pd(HH, XH, 0x03);
      T  = _mm256_permute_pd(HH, 0x00);
      XH = _mm256_unpackhi_pd(XH, T);

      /* concatenate XB and BB and rotate right by one */
      XB = _mm256_permute2f128_pd(BB, XB, 0x03);
      T  = _mm256_permute_pd(BB, 0x00);
      XB = _mm256_unpackhi_pd(XB, T);

      /* concatenate XE and EE and rotate right by one */
      XE = _mm256_permute2f128_pd(EE, XE, 0x03);
      T  = _mm256_permute_pd(EE, 0x00);
      XE = _mm256_unpackhi_pd(XE, T);

      T = _mm256_mul_pd(XH,PHH);
#ifdef HAVE_FMA
      T = _mm256_fmadd_pd(XB,PBH,T);
      T = _mm256_fmadd_pd(XE,PEH,T);
#else
      T = _mm256_add_pd(T, _mm256_mul_pd(XB,PBH));
      T = _mm256_add_pd(T, _mm256_mul_pd(XE,PEH));
#endif
      _mm256_store_pd(chh+i,T);

      /* store T1 to compute cee later (H component) */
      T1 = _mm256_permute2f128_pd(T, XHH, 0x03);
      T3 = _mm256_permute_pd(T, 0x00);
      T1 = _mm256_unpackhi_pd(T1, T3);
      XHH = T;

      XH = HH;
      XB = BB;
      XE = EE;

      /* attempt to parallelize these two lines */
      //cee[i] = chh[i-1]*phe + cbb[i-1]*pbe + cee[i-1]*pee;
      //cee[i+1] = chh[i]*phe + cbb[i]*pbe + cee[i]*pee;
      //cee[i+2] = chh[i+1]*phe + cbb[i+1]*pbe + cee[i+1]*pee;
      //cee[i+3] = chh[i+2]*phe + cbb[i+2]*pbe + cee[i+2]*pee;

      /* first (and valid) of the four E component used to compute cee */
#ifdef HAVE_AVX2
      T3 = _mm256_permute4x64_pd(XEE, 0x03);
#else
      T3 = _mm256_permute2f128_pd(XEE,XEE, 0x01);
      T3 = _mm256_permute_pd(T3, 0x03);
#endif

//      //T1  = _mm_loadu_pd (chh+i-1);
      T1 = _mm256_mul_pd(T1,PHE);
//      //T2 = _mm_loadu_pd(cbb+i-1);
#ifdef HAVE_FMA
      T4 = _mm256_fmadd_pd(T2,PBE,T1);
#else
      T2 = _mm256_mul_pd(T2,PBE);
      T4 = _mm256_add_pd(T1,T2);
#endif
//      //T3 = _mm_loadu_pd(cee+i-1);

#ifdef HAVE_FMA
      T1 = _mm256_fmadd_pd(T3,PEE,T4);
#else
      T3 = _mm256_mul_pd(T3,PEE);
      
      /* compute correct cee[i] value */
      T1 = _mm256_add_pd(T4,T3);
#endif

      /* move the correct value to the slot [127:64] */
      T2 = _mm256_permute_pd(T1, 0x00);
#ifdef HAVE_FMA
      T2 = _mm256_fmadd_pd(T2,PEE,T4);
#else
      T2 = _mm256_mul_pd(T2, PEE);
      /* compute correct cee[i+1] */
      T2 = _mm256_add_pd(T2,T4);
#endif
      /* blend T1[63:0] with T2[127:64] */
      T1 = _mm256_blend_pd(T1,T2, 0x02);
        
      /* move the correct value to slot [191:128] */
#ifdef HAVE_AVX2
      T3 = _mm256_permute4x64_pd(T2,0x55);
#else
      T3 = _mm256_permute2f128_pd(T2,T2,0x01);
      T3 = _mm256_permute_pd(T3, 0x0C);
#endif
#ifdef HAVE_FMA
      T3 = _mm256_fmadd_pd(T3,PEE,T4);
#else
      T3 = _mm256_mul_pd(T3,PEE);
      /* compute correct cee[i+2] */
      T3 = _mm256_add_pd(T3,T4);
#endif

      /* move the correct value to slot [255:192] */
      T2 = _mm256_permute_pd(T3, 0x00);
#ifdef HAVE_FMA
      T2 = _mm256_fmadd_pd(T2,PEE,T4);
#else
      T2 = _mm256_mul_pd(T2,PEE);
      /* compute correct cee[i+3] */
      T2 = _mm256_add_pd(T2,T4);
#endif
      /* blend T3[191:128] with T2[255:192] */
      T2 = _mm256_blend_pd(T3,T2,0x08);

      /* blend T1[127:0] with T2[255:128] */
      T1 = _mm256_blend_pd(T1, T2, 0x0C);

      _mm256_store_pd(cee+i,T1);

      XEE = T1;
    }

    hh = chh; bb = cbb; ee = cee;
  }
  
  return (chh[m+VEC_SIZE-1]*phh + cbb[m+VEC_SIZE-1]*pbh + cee[m+VEC_SIZE-1]*peh);
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

  h = strale_fid1d_matrix_avx(9000,
                              9000,
                              1,
                              0,
                              0);
  
  printf ("-> Result-avx: %f\n", h);

  return (EXIT_SUCCESS);
}

