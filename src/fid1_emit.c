#include "fid.h"
#include <string.h>

static double * dd = NULL;
static long ddlen = 0;
#ifdef SCALING
static int * scale = NULL;
#endif

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
     
     FID 1 model implementation (double) with emissions, homologies, no scaling

     Computation of the dynamic programming matrix is performed on a column by
     column basis. In the end, the log probability is returned based on the
     content of rightn.

     input:
     s1: partials of sequence 1. Eeach element has four entries (order A,C,G,T)
     s2: partials of sequence 2. Eeach element has four entries (order A,C,G,T)
     m: number of elements in s1 
     n: number of elements in s2
     pstart: starting probability for homology (for cell 0,0)
     sm: 4x4 substitution matrix - linear array of concatenated rows (order
         A,C,G,T in both dimensions)
     leftn: starting event (0 = Homology, 1 = Birth, 2 = Extinction)
     rightn: ending event (0 = Homology, 1 = Birth, 2 = Extinction)

     matrix format:
     The matrix (dd) is a linear array of size 3*n*(m+1) doubles. It is
     arranged in columns. Each component H,B,E of the cells of a column are
     stored consequently.  Therefore, elements dd[3*(m+1)*i+0..3*(m+1)*i+m] are
     the H components, dd[3*(m+1)*i+m+1..3*(m+1)*i+2*m+1] the B components, and
     dd[3*(m+1)*i+2*m+2..3*m*i+3*m+2] the E components of the i-th column, for
     0 <= i < n. Note that the initialization row is included in the matrix,
     but not the initialization column.

     output:
     Depending on the content of rightn, the function returns the calculated
     log probability of the cell diagonally bottom-right, right or below the
     bottom-right cell of the DP matrix, depending if rightn was set to H, B or
     E, repsectively.

*/


double strale_fid1d_emit(const double * s1, const double * s2,
                         const int m, const int n,
                         double pstart,
                         const double * sm,
                         unsigned int leftn, unsigned int rightn,
                         const double * freqs)
{
#ifdef SCALING
  double SCALE_THRESHOLD = sqrt(__DBL_MIN__);
  double SCALE_FACTOR    = __DBL_MAX__;
#endif

  double * ee = NULL;    /* pointer to the E section of the previous column */
  double * hh = NULL;    /* pointer to the H section of the previous column */
  double * bb = NULL;    /* pointer to the B section of the previous column */

  double * chh;          /* pointer to the E section of the current column */ 
  double * cbb;          /* pointer to the H section of the current column */
  double * cee;          /* pointer to the B section of the current column */
  double hemit;          /* emission probability (homology) */
  double bemit;          /* emission probability (birth) */
  double eemit;          /* emission probability (extinction) */

#ifdef SCALING
  int * hscale;
  int * bscale;
  int * escale;
#endif

  double prod0;          /* temp storage for computing products in H emissions */
  double prod1;
  double prod2;
  double prod3;

  double sum0;           /* temp storage for computing sums in H emissions */
  double sum1;
  double sum2;
  double sum3;

  long i, j;

  /* reallocate matrix if necessary */
  if (3*n*(m+1) > ddlen)
  {
    free(dd);
    dd = xmalloc(3*n*(m+1)*sizeof(double), FID_ALIGNMENT_SSE);
#ifdef SCALING
    free(scale);
    scale   = xmalloc(3*n*(m+1)*sizeof(int), FID_ALIGNMENT_SSE);
#endif
  }

#ifdef SCALING
  memset(scale,0,3*n*(m+1)*sizeof(int));
#endif

  /* init pointers to elements */
  chh = hh = dd;
  cbb = bb = hh + m+1;
  cee = ee = bb + m+1;

#ifdef SCALING
  hscale = scale;
  bscale = hscale + m+1;
  escale = bscale + m+1;
#endif
  
  bemit = freqs[0]*s2[0*4+0] + freqs[1]*s2[0*4+1] + freqs[2]*s2[0*4+2] + freqs[3]*s2[0*4+3];
  eemit = freqs[0]*s1[0*4+0] + freqs[1]*s1[0*4+1] + freqs[2]*s1[0*4+2] + freqs[3]*s1[0*4+3];
  //double b0 = ph*phb + pb*pbb + pe*peb;
  double b0 = pstart*phb*bemit;
  double e0 = pstart*phe*eemit;

  /* precompute cell (0,1) */
  chh[0] = 0;
  cee[0] = 0;
  cbb[0] = b0;

  /* precompute cell (1,1) */
  sum0  = sm[0*4+0]*s2[0*4+0] + sm[0*4+1]*s2[0*4+1] + sm[0*4+2]*s2[0*4+2] + sm[0*4+3]*s2[0*4+3];
  sum1  = sm[1*4+0]*s2[0*4+0] + sm[1*4+1]*s2[0*4+1] + sm[1*4+2]*s2[0*4+2] + sm[1*4+3]*s2[0*4+3];
  sum2  = sm[2*4+0]*s2[0*4+0] + sm[2*4+1]*s2[0*4+1] + sm[2*4+2]*s2[0*4+2] + sm[2*4+3]*s2[0*4+3];
  sum3  = sm[3*4+0]*s2[0*4+0] + sm[3*4+1]*s2[0*4+1] + sm[3*4+2]*s2[0*4+2] + sm[3*4+3]*s2[0*4+3];
  
  prod0 = freqs[0]*s1[0*4+0];
  prod1 = freqs[0]*s1[0*4+1];
  prod2 = freqs[0]*s1[0*4+2];
  prod3 = freqs[0]*s1[0*4+3];

  hemit = prod0*sum0 + prod1*sum1 + prod2*sum2 + prod3*sum3;

  chh[1] = pstart*phh*hemit;
  cbb[1] = peb*e0*bemit;
  cee[1] = pbe*b0*eemit;

#ifdef SCALING
  if (chh[1] < SCALE_THRESHOLD)
   {
     printf ("Scaling H (1,1)\n");
     hscale[1] = 1;
     chh[1] *= SCALE_FACTOR;

   }
  if (cbb[1] < SCALE_THRESHOLD)
  {
     printf ("Scaling B (1,1)\n");
     bscale[1] = 1;
     cbb[1] *= SCALE_FACTOR;
  }
  if (cee[1] < SCALE_THRESHOLD)
  {
     printf ("Scaling E (1,1)\n");
     escale[1] = 1;
     cee[1] *= SCALE_FACTOR;
  }
#endif

  /* precompute the rest of column 1, i.e. cells (i,1) */
  for (i = 2; i <= m; ++i)
  {
    prod0 = freqs[0]*s1[(i-1)*4+0];
    prod1 = freqs[0]*s1[(i-1)*4+1];
    prod2 = freqs[0]*s1[(i-1)*4+2];
    prod3 = freqs[0]*s1[(i-1)*4+3];

    hemit = prod0*sum0 + prod1*sum1 + prod2*sum2 + prod3*sum3;
    eemit  = freqs[0]*s1[(i-1)*4+0] + freqs[1]*s1[(i-1)*4+1] + freqs[2]*s1[(i-1)*4+2] + freqs[3]*s2[(i-1)*4+3];

    cee[i] = (chh[i-1]*phe + cbb[i-1]*pbe + cee[i-1]*pee)*eemit;
    chh[i] = peh*e0*hemit;
    e0    *= pee*eemit;
    cbb[i] = peb*e0*bemit;


#ifdef SCALING
    if (chh[i] < SCALE_THRESHOLD)
     {
       printf ("Scaling H (%d,1)\n", i);
       hscale[i] = 1;
       chh[i] *= SCALE_FACTOR;
     }
    if (cbb[i] < SCALE_THRESHOLD)
    {
       printf ("Scaling B (%d,1)\n", i);
       bscale[i] = 1;
       cbb[i] *= SCALE_FACTOR;
    }
    if (cee[i] < SCALE_THRESHOLD)
    {
       printf ("Scaling E (%d,1)\n", i);
       escale[i] = 1;
       cee[i] *= SCALE_FACTOR;
    }
#endif
  }

  /* iterate through columns starting from the second */
  for (j = 1; j < n; ++j)
  {
    /* compute b emission for this column */
    bemit    = freqs[0]*s2[j*4+0] + freqs[1]*s2[j*4+1] + freqs[2]*s2[j*4+2] + freqs[3]*s2[j*4+3];
    
    /* compute part of the H emission that does not change for the column */
    sum0  = sm[0*4+0]*s2[j*4+0] + sm[0*4+1]*s2[j*4+1] + sm[0*4+2]*s2[j*4+2] + sm[0*4+3]*s2[j*4+3];
    sum1  = sm[1*4+0]*s2[j*4+0] + sm[1*4+1]*s2[j*4+1] + sm[1*4+2]*s2[j*4+2] + sm[1*4+3]*s2[j*4+3];
    sum2  = sm[2*4+0]*s2[j*4+0] + sm[2*4+1]*s2[j*4+1] + sm[2*4+2]*s2[j*4+2] + sm[2*4+3]*s2[j*4+3];
    sum3  = sm[3*4+0]*s2[j*4+0] + sm[3*4+1]*s2[j*4+1] + sm[3*4+2]*s2[j*4+2] + sm[3*4+3]*s2[j*4+3];

    /* point each entry to the beginning (first element) */
    chh =  ee + m+1;
    cbb = chh + m+1;
    cee = cbb + m+1;

#ifdef SCALING
    hscale = escale + m+1;
    bscale = hscale + m+1;
    escale = bscale + m+1;
#endif

    /* compute element (0,j) */
    chh[0] = 0;
    cee[0] = 0;
    cbb[0] = bb[0]*pbb*bemit;

    /* iterate cells of a column */
    for (i = 1; i <= m; ++i)
    {
      prod0 = freqs[0]*s1[(i-1)*4+0];
      prod1 = freqs[0]*s1[(i-1)*4+1];
      prod2 = freqs[0]*s1[(i-1)*4+2];
      prod3 = freqs[0]*s1[(i-1)*4+3];

      hemit = prod0*sum0 + prod1*sum1 + prod2*sum2 + prod3*sum3;
      eemit  = freqs[0]*s1[(i-1)*4+0] + freqs[1]*s1[(i-1)*4+1] + freqs[2]*s1[(i-1)*4+2] + freqs[3]*s2[(i-1)*4+3];
      
      cbb[i] = (hh[i]*phb    + bb[i]*pbb    +    ee[i]*peb) * bemit;
      chh[i] = (hh[i-1]*phh  + bb[i-1]*pbh  +  ee[i-1]*peh) * hemit;
      cee[i] = (chh[i-1]*phe + cbb[i-1]*pbe + cee[i-1]*pee) * eemit;

#ifdef SCALING
      if (chh[i] < SCALE_THRESHOLD)
       {
         printf ("Scaling H (%d,%d)\n", i, j);
         hscale[i] = 1;
         chh[i] *= SCALE_FACTOR;
       }
      if (cbb[i] < SCALE_THRESHOLD)
      {
         printf ("Scaling B (%d,%d)\n", i, j);
         bscale[i] = 1;
         cbb[i] *= SCALE_FACTOR;
      }
      if (cee[i] < SCALE_THRESHOLD)
      {
         printf ("Scaling E (%d,%d)\n", i, j);
         escale[i] = 1;
         cee[i] *= SCALE_FACTOR;
      }
#endif
    }

    hh = chh; bb = cbb; ee = cee;
  }

 
  if (rightn == 0) return (chh[m]*phh + cbb[m]*pbh + cee[m]*peh);
  if (rightn == 1) return (chh[m]*phb + cbb[m]*pbb + cee[m]*peb);
  
  /* rightn == 2 */
  return (chh[m]*phe + cbb[m]*pbe + cee[m]*pee);
}

/* code taken from STrAlE for setting up transition probability among events */
static void compute_probs(double lambda, double gamma1)
{
  if ( (lambda) >= 1e-7) 
   {
     double el = exp(-(lambda));
     double sl = (1+lambda);
  
     phh = 1-( (1+lambda-el)/(gamma1*sl) );
     phb = lambda/(gamma1*sl);
     phe = (1-el)/(gamma1*sl);
  
     pbh = el/(gamma1*sl);
     pbb = ((gamma1*sl)-1)/(gamma1*sl);
     pbe = phe;
  
     peh = (lambda*el)/(gamma1*sl*(1-el));
     peb = (1-el*sl)/(gamma1*sl*(1-el));
     pee = (gamma1*sl-1)/(gamma1*sl);
   }
  else 
   { //approximations to avoid nan and inf values
     double el = 1-lambda-((lambda*lambda)/2)-((lambda*lambda*lambda)/6);
     double sl = (1+lambda);
  
     phh = 1-( (1+lambda-el)/(gamma1*sl) );
     phb = lambda/(gamma1*sl);
     phe = (1-el)/(gamma1*sl);
  
     pbh = el/(gamma1*sl);
     pbb = ((gamma1*sl)-1)/(gamma1*sl);
     pbe = phe;
  
     peh = (el)/(gamma1*sl);
     peb = lambda/(gamma1*sl);
     pee = (gamma1*sl-1)/(gamma1*sl);
  }
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


double strale_fid1de_run(const double * s1, 
                       const double * s2,
                       unsigned int s1start, 
                       unsigned int s1end,
                       unsigned int s2start, 
                       unsigned int s2end,
                       double t, 
                       double lambda, 
                       double gamma,
                       unsigned int leftn, 
                       unsigned int rightn,
                       const double * substmatrix,
                       const double * basefreqs,
                       int allow_homologies,
                       double startprob)
{
  double prob;

  /* compute the 9 probability entries and output them (order H,B,E) */
  compute_probs(t*lambda, gamma);
  strale_fid1d_dump_tpm();

  if (rightn != 0 && rightn != 1 && rightn != 2) return (0);

  if (allow_homologies)
    prob = strale_fid1d_emit(s1, s2,
                             s1end - s1start, s2end - s2start,
                             startprob,
                             substmatrix,
                             leftn, rightn,
                             basefreqs);
  else
  {
    /* not yet implemented. Here will be the call to the function not allowing
       homologies. */
    prob = 0;
  }

  return prob;
}
                 

int main(int argc, char * argv[])
{
  double h;
  int i;

  /* set number of elements (not elements*4) for the two sequences */
  int s1_size = 50;
  int s2_size = 50;

  /* set up base frequencies */
  double basefreqs[4] = {0.25, 0.25, 0.25, 0.25};

  /* set up substitution matrix */
  double substmatrix[16] = {  0.1,  0.2,  0.3,  0.4,
                             0.23, 0.48, 0.19,  0.1,
                              0.3,  0.3,  0.3,  0.1,
                             0.25, 0.25, 0.25, 0.25};

                                

  /* allocate memory for all elements*entries (4) of the sequences */
  double * s1 = xmalloc(s1_size*4*sizeof(double), FID_ALIGNMENT_SSE);
  double * s2 = xmalloc(s2_size*4*sizeof(double), FID_ALIGNMENT_SSE);

  /* set some random values as the partials for sequence 1, and make sure that
     the four entries for each element sum up to 1 */
  for(i = 0; i < 4*s1_size; i+=4)
  {
    s1[i+0] = drand48();
    s1[i+1] = drand48()*(1 - s1[i+0]);
    s1[i+2] = drand48()*(1 - s1[i+0] - s1[i+1]);
    s1[i+3] = 1 - s1[i+0] - s1[i+1] - s1[i+2];
  }
  
  /* set some random values as the partials for sequence 2, and make sure that
     the four entries for each element sum up to 1 */
  for(i = 0; i < 4*s2_size; i+=4)
  {
    s2[i+0] = drand48();
    s2[i+1] = drand48()*(1 - s2[i+0]);
    s2[i+2] = drand48()*(1 - s2[i+0] - s2[i+1]);
    s2[i+3] = 1 - s2[i+0] - s2[i+1] - s2[i+2];
  }

  h = strale_fid1de_run(s1,             /* sequence 1 */
                        s2,             /* sequence 2 */
                        0,              /* sequence 1 start */
                        s1_size,        /* sequence 1 end */
                        0,              /* sequence 1 start */
                        s2_size,        /* sequence 2 end */
                        0.208025,       /* tau */
                        0.02,           /* lambda */
                        2,              /* gamma */
                        0,              /* left neighbor */
                        0,              /* right neighbor */
                        substmatrix,    /* substitution matrix */
                        basefreqs,      /* base frequencies */
                        1,              /* allow homologies? */
                        1);             /* left neighbor start probability */

  printf ("-> Prob: %.200f\n", h);
  printf ("Logprob: %f\n", log(h));
  
  return (EXIT_SUCCESS);
}
