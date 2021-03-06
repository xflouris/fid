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

static void set_sequence_1 (double * s);
static void set_sequence_2 (double * s);

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

  long i, j, k1, k2;

  /* set the starting probability of cell (0,0) given leftn and pstart */
  double startprob[3] = {0,0,0};

  assert(leftn >=0 && leftn <= 3);
  startprob[leftn] = pstart;  


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
  double b0 = (startprob[0]*phb + startprob[1]*pbb + startprob[2]*peb)*bemit;
  double e0 = (startprob[0]*phe + startprob[1]*pbe + startprob[2]*pee)*eemit;

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
  for (k1 = 4, i = 2; i <= m; ++i, k1 += 4)
  {
    prod0 = freqs[0]*s1[k1+0];
    prod1 = freqs[0]*s1[k1+1];
    prod2 = freqs[0]*s1[k1+2];
    prod3 = freqs[0]*s1[k1+3];

    hemit = prod0*sum0 + prod1*sum1 + prod2*sum2 + prod3*sum3;
    eemit  = freqs[0]*s1[k1+0] + freqs[1]*s1[k1+1] + freqs[2]*s1[k1+2] + freqs[3]*s1[k1+3];

    cee[i] = (chh[i-1]*phe + cbb[i-1]*pbe + cee[i-1]*pee)*eemit;
    chh[i] = peh*e0*hemit;
    e0    *= pee*eemit;
    cbb[i] = peb*e0*bemit;


#ifdef SCALING
    if (chh[i] < SCALE_THRESHOLD)
     {
       printf ("Scaling H (%ld,1)\n", i);
       hscale[i] = 1;
       chh[i] *= SCALE_FACTOR;
     }
    if (cbb[i] < SCALE_THRESHOLD)
    {
       printf ("Scaling B (%ld,1)\n", i);
       bscale[i] = 1;
       cbb[i] *= SCALE_FACTOR;
    }
    if (cee[i] < SCALE_THRESHOLD)
    {
       printf ("Scaling E (%ld,1)\n", i);
       escale[i] = 1;
       cee[i] *= SCALE_FACTOR;
    }
#endif
  }

  /* iterate through columns starting from the second */
  for (k2 = 4, j = 1; j < n; ++j, k2 += 4)
  {
    /* compute b emission for this column */
    bemit    = freqs[0]*s2[k2+0] + freqs[1]*s2[k2+1] + freqs[2]*s2[k2+2] + freqs[3]*s2[k2+3];
    
    /* compute part of the H emission that does not change for the column */
    sum0  = sm[0*4+0]*s2[k2+0] + sm[0*4+1]*s2[k2+1] + sm[0*4+2]*s2[k2+2] + sm[0*4+3]*s2[k2+3];
    sum1  = sm[1*4+0]*s2[k2+0] + sm[1*4+1]*s2[k2+1] + sm[1*4+2]*s2[k2+2] + sm[1*4+3]*s2[k2+3];
    sum2  = sm[2*4+0]*s2[k2+0] + sm[2*4+1]*s2[k2+1] + sm[2*4+2]*s2[k2+2] + sm[2*4+3]*s2[k2+3];
    sum3  = sm[3*4+0]*s2[k2+0] + sm[3*4+1]*s2[k2+1] + sm[3*4+2]*s2[k2+2] + sm[3*4+3]*s2[k2+3];

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
    for (k1 = 0, i = 1; i <= m; ++i, k1 += 4)
    {
      prod0 = freqs[0]*s1[k1+0];
      prod1 = freqs[0]*s1[k1+1];
      prod2 = freqs[0]*s1[k1+2];
      prod3 = freqs[0]*s1[k1+3];

      hemit = prod0*sum0 + prod1*sum1 + prod2*sum2 + prod3*sum3;
      eemit  = freqs[0]*s1[k1+0] + freqs[1]*s1[k1+1] + freqs[2]*s1[k1+2] + freqs[3]*s2[k1+3];
      
      cbb[i] = (hh[i]*phb    + bb[i]*pbb    +    ee[i]*peb) * bemit;
      chh[i] = (hh[i-1]*phh  + bb[i-1]*pbh  +  ee[i-1]*peh) * hemit;
      cee[i] = (chh[i-1]*phe + cbb[i-1]*pbe + cee[i-1]*pee) * eemit;

#ifdef SCALING
      if (chh[i] < SCALE_THRESHOLD)
       {
         printf ("Scaling H (%ld,%ld)\n", i, j);
         hscale[i] = 1;
         chh[i] *= SCALE_FACTOR;
       }
      if (cbb[i] < SCALE_THRESHOLD)
      {
         printf ("Scaling B (%ld,%ld)\n", i, j);
         bscale[i] = 1;
         cbb[i] *= SCALE_FACTOR;
      }
      if (cee[i] < SCALE_THRESHOLD)
      {
         printf ("Scaling E (%ld,%ld)\n", i, j);
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
    prob = strale_fid1d_emit(s1+s1start*4, s2+s2start*4,
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

  /* set number of elements (not elements*4) for the two sequences */
  int s1_size = 75;
  int s2_size = 75;

  /* set up base frequencies */
  double basefreqs[4] = {0.25, 0.25, 0.25, 0.25};

  /* set up substitution matrix */
  double substmatrix[16] = { 0.989897,   0.00247589, 0.00493067, 0.0026961,
                             0.00247589, 0.986566,   0.00314022, 0.00781766,
                             0.00493067, 0.00314022, 0.990568,   0.00136104,
                             0.0026961,  0.00781766, 0.00136104, 0.988125};

                                

  /* allocate memory for all elements*entries (4) of the sequences */
  double * s1 = xmalloc(s1_size*4*sizeof(double), FID_ALIGNMENT_SSE);
  double * s2 = xmalloc(s2_size*4*sizeof(double), FID_ALIGNMENT_SSE);

  set_sequence_1(s1);
  set_sequence_2(s2);

  h = strale_fid1de_run(s1,           /* sequence 1 */
                        s2,           /* sequence 2 */
                        1,              /* sequence 1 start */
                        10,             /* sequence 1 end */
                        1,              /* sequence 1 start */
                        10,             /* sequence 2 end */
                        0.0452143,      /* tau */
                        0.0026875,      /* lambda */
                        2.5,            /* gamma */
                        0,              /* left neighbor */
                        2,              /* right neighbor */
                        substmatrix,    /* substitution matrix */
                        basefreqs,      /* base frequencies */
                        1,              /* allow homologies? */
                        1);             /* left neighbor start probability */

  printf ("-> Prob: %.200f\n", h);
  printf ("Logprob: %f\n", log(h));
  
  return (EXIT_SUCCESS);
}

static void set_sequence_1 (double * s)
{
  s[0] = 1;
  s[1] = 1;
  s[2] = 1;
  s[3] = 1;
  s[4] = 1.58485e-07;
  s[5] = 6.48549e-06;
  s[6] = 0.00482394;
  s[7] = 3.19091e-08;
  s[8] = 1.62398e-05;
  s[9] = 8.99627e-08;
  s[10] = 0.00757243;
  s[11] = 2.07098e-08;
  s[12] = 1.60808e-09;
  s[13] = 6.55967e-07;
  s[14] = 5.49195e-08;
  s[15] = 1.67577e-05;
  s[16] = 1.17792e-07;
  s[17] = 2.99707e-05;
  s[18] = 1.58618e-09;
  s[19] = 5.66209e-07;
  s[20] = 6.41107e-08;
  s[21] = 0.00472963;
  s[22] = 6.47532e-06;
  s[23] = 3.03119e-07;
  s[24] = 1.60094e-05;
  s[25] = 6.47402e-06;
  s[26] = 0.94029;
  s[27] = 1.21187e-06;
  s[28] = 1.60094e-05;
  s[29] = 6.47402e-06;
  s[30] = 0.94029;
  s[31] = 1.21187e-06;
  s[32] = 4.71537e-06;
  s[33] = 3.99984e-05;
  s[34] = 1.19841e-06;
  s[35] = 0.92535;
  s[36] = 5.11893e-05;
  s[37] = 3.81444e-08;
  s[38] = 5.36668e-05;
  s[39] = 1.79623e-08;
  s[40] = 1.60094e-05;
  s[41] = 6.47402e-06;
  s[42] = 0.94029;
  s[43] = 1.21187e-06;
  s[44] = 3.85419e-08;
  s[45] = 1.93998e-06;
  s[46] = 1.00489e-08;
  s[47] = 0.00736094;
  s[48] = 1.84881e-08;
  s[49] = 3.53345e-08;
  s[50] = 4.09574e-06;
  s[51] = 3.9165e-06;
  s[52] = 3.85419e-08;
  s[53] = 1.93998e-06;
  s[54] = 1.00489e-08;
  s[55] = 0.00736094;
  s[56] = 0.936161;
  s[57] = 4.00707e-06;
  s[58] = 1.59603e-05;
  s[59] = 4.75371e-06;
  s[60] = 1.62398e-05;
  s[61] = 8.99627e-08;
  s[62] = 0.00757243;
  s[63] = 2.07098e-08;
  s[64] = 1.84881e-08;
  s[65] = 3.53345e-08;
  s[66] = 4.09574e-06;
  s[67] = 3.9165e-06;
  s[68] = 1.84881e-08;
  s[69] = 3.53345e-08;
  s[70] = 4.09574e-06;
  s[71] = 3.9165e-06;
  s[72] = 1.12219e-07;
  s[73] = 1.07254e-07;
  s[74] = 0.00210469;
  s[75] = 1.21133e-06;
  s[76] = 1.00912e-07;
  s[77] = 2.83633e-09;
  s[78] = 3.19196e-08;
  s[79] = 5.68716e-06;
  s[80] = 5.43768e-08;
  s[81] = 3.07158e-07;
  s[82] = 1.2102e-06;
  s[83] = 0.00207952;
  s[84] = 0.936161;
  s[85] = 4.00707e-06;
  s[86] = 1.59603e-05;
  s[87] = 4.75371e-06;
  s[88] = 0.004679;
  s[89] = 2.0647e-08;
  s[90] = 4.92073e-07;
  s[91] = 2.40713e-08;
  s[92] = 5.74638e-08;
  s[93] = 2.96842e-05;
  s[94] = 8.46627e-10;
  s[95] = 4.98247e-07;
  s[96] = 1.90906e-05;
  s[97] = 9.80175e-08;
  s[98] = 2.63735e-07;
  s[99] = 1.36769e-09;
  s[100] = 3.85419e-08;
  s[101] = 1.93998e-06;
  s[102] = 1.00489e-08;
  s[103] = 0.00736094;
  s[104] = 2.26188e-07;
  s[105] = 2.06757e-09;
  s[106] = 1.05909e-05;
  s[107] = 3.45121e-08;
  s[108] = 1.29504e-05;
  s[109] = 1.3366e-05;
  s[110] = 3.88584e-08;
  s[111] = 5.35171e-08;
  s[112] = 6.41107e-08;
  s[113] = 0.00472963;
  s[114] = 6.47532e-06;
  s[115] = 3.03119e-07;
  s[116] = 7.77283e-10;
  s[117] = 4.90635e-07;
  s[118] = 2.04064e-08;
  s[119] = 1.64701e-05;
  s[120] = 1.12219e-07;
  s[121] = 1.07254e-07;
  s[122] = 0.00210469;
  s[123] = 1.21133e-06;
  s[124] = 5.74638e-08;
  s[125] = 2.96842e-05;
  s[126] = 8.46627e-10;
  s[127] = 4.98247e-07;
  s[128] = 8.40372e-08;
  s[129] = 0.0117516;
  s[130] = 1.05708e-07;
  s[131] = 4.08362e-05;
  s[132] = 5.43768e-08;
  s[133] = 3.07158e-07;
  s[134] = 1.2102e-06;
  s[135] = 0.00207952;
  s[136] = 1.60094e-05;
  s[137] = 6.47402e-06;
  s[138] = 0.94029;
  s[139] = 1.21187e-06;
  s[140] = 1.23731e-07;
  s[141] = 3.44557e-09;
  s[142] = 4.09708e-08;
  s[143] = 5.73484e-06;
  s[144] = 1.29504e-05;
  s[145] = 1.3366e-05;
  s[146] = 3.88584e-08;
  s[147] = 5.35171e-08;
  s[148] = 0.00754744;
  s[149] = 6.46142e-08;
  s[150] = 1.62356e-05;
  s[151] = 5.44315e-08;
  s[152] = 5.74638e-08;
  s[153] = 2.96842e-05;
  s[154] = 8.46627e-10;
  s[155] = 4.98247e-07;
  s[156] = 1.84881e-08;
  s[157] = 3.53345e-08;
  s[158] = 4.09574e-06;
  s[159] = 3.9165e-06;
  s[160] = 4.00649e-06;
  s[161] = 0.00373667;
  s[162] = 8.92455e-08;
  s[163] = 3.54104e-07;
  s[164] = 0.00234557;
  s[165] = 6.21943e-08;
  s[166] = 4.17441e-08;
  s[167] = 1.2875e-08;
  s[168] = 1.60094e-05;
  s[169] = 6.47402e-06;
  s[170] = 0.94029;
  s[171] = 1.21187e-06;
  s[172] = 2.3324e-08;
  s[173] = 9.98528e-09;
  s[174] = 0.00129326;
  s[175] = 1.04647e-08;
  s[176] = 5.35849e-05;
  s[177] = 3.82693e-08;
  s[178] = 5.13118e-05;
  s[179] = 1.79576e-08;
  s[180] = 1.84881e-08;
  s[181] = 3.53345e-08;
  s[182] = 4.09574e-06;
  s[183] = 3.9165e-06;
  s[184] = 5.35849e-05;
  s[185] = 3.82693e-08;
  s[186] = 5.13118e-05;
  s[187] = 1.79576e-08;
  s[188] = 5.43768e-08;
  s[189] = 3.07158e-07;
  s[190] = 1.2102e-06;
  s[191] = 0.00207952;
  s[192] = 6.41107e-08;
  s[193] = 0.00472963;
  s[194] = 6.47532e-06;
  s[195] = 3.03119e-07;
  s[196] = 1.62398e-05;
  s[197] = 8.99627e-08;
  s[198] = 0.00757243;
  s[199] = 2.07098e-08;
  s[200] = 3.86025e-08;
  s[201] = 2.149e-05;
  s[202] = 2.08171e-05;
  s[203] = 3.38488e-08;
  s[204] = 1.84151e-09;
  s[205] = 2.63635e-07;
  s[206] = 6.78782e-06;
  s[207] = 6.78866e-08;
  s[208] = 1.12219e-07;
  s[209] = 1.07254e-07;
  s[210] = 0.00210469;
  s[211] = 1.21133e-06;
  s[212] = 1.62398e-05;
  s[213] = 8.99627e-08;
  s[214] = 0.00757243;
  s[215] = 2.07098e-08;
  s[216] = 4.77463e-06;
  s[217] = 3.5512e-07;
  s[218] = 2.05668e-08;
  s[219] = 0.00409609;
  s[220] = 4.77463e-06;
  s[221] = 3.5512e-07;
  s[222] = 2.05668e-08;
  s[223] = 0.00409609;
  s[224] = 4.00649e-06;
  s[225] = 0.00373667;
  s[226] = 8.92455e-08;
  s[227] = 3.54104e-07;
  s[228] = 2.06991e-05;
  s[229] = 1.63263e-09;
  s[230] = 1.15253e-07;
  s[231] = 6.22702e-08;
  s[232] = 0.00379859;
  s[233] = 4.01171e-06;
  s[234] = 1.58473e-07;
  s[235] = 9.58636e-08;
  s[236] = 1.59824e-05;
  s[237] = 5.30803e-08;
  s[238] = 1.79713e-08;
  s[239] = 1.52345e-05;
  s[240] = 9.53491e-08;
  s[241] = 4.08598e-05;
  s[242] = 3.18263e-08;
  s[243] = 0.0118425;
  s[244] = 3.94633e-06;
  s[245] = 0.915908;
  s[246] = 6.35625e-06;
  s[247] = 3.97139e-05;
  s[248] = 8.40372e-08;
  s[249] = 0.0117516;
  s[250] = 1.05708e-07;
  s[251] = 4.08362e-05;
  s[252] = 1.31959e-08;
  s[253] = 0.00292177;
  s[254] = 1.26925e-07;
  s[255] = 1.29466e-07;
  s[256] = 1.84881e-08;
  s[257] = 3.53345e-08;
  s[258] = 4.09574e-06;
  s[259] = 3.9165e-06;
  s[260] = 5.13731e-08;
  s[261] = 0.000126863;
  s[262] = 3.28692e-08;
  s[263] = 0.000133485;
  s[264] = 8.40372e-08;
  s[265] = 0.0117516;
  s[266] = 1.05708e-07;
  s[267] = 4.08362e-05;
  s[268] = 3.94633e-06;
  s[269] = 0.915908;
  s[270] = 6.35625e-06;
  s[271] = 3.97139e-05;
  s[272] =  1.12219e-07;
  s[273] =  1.07254e-07;
  s[274] =  0.00210469;
  s[275] =  1.21133e-06;
  s[276] = 4.91727e-07;
  s[277] = 3.30911e-08;
  s[278] = 0.00469642;
  s[279] = 6.26005e-09;
  s[280] = 0.004679;
  s[281] = 2.0647e-08;
  s[282] = 4.92073e-07;
  s[283] = 2.40713e-08;
  s[284] = 8.10097e-06;
  s[285] = 0.936125;
  s[286] = 1.30265e-05;
  s[287] = 8.12073e-05;
  s[288] = 8.10097e-06;
  s[289] = 0.936125;
  s[290] = 1.30265e-05;
  s[291] = 8.12073e-05;
  s[292] = 5.21154e-05;
  s[293] = 1.94029e-07;
  s[294] = 0.00478487;
  s[295] = 5.94685e-08;
  s[296] = 1;
  s[297] = 1;
  s[298] = 1;
  s[299] = 1;
}

static void set_sequence_2(double * s)
{
  s[0] = 1;
  s[1] = 1;
  s[2] = 1;
  s[3] = 1;
  s[4] = 0;
  s[5] = 0;
  s[6] = 1;
  s[7] = 0;
  s[8] = 0;
  s[9] = 0;
  s[10] = 1;
  s[11] = 0;
  s[12] = 0;
  s[13] = 1;
  s[14] = 0;
  s[15] = 0;
  s[16] = 0;
  s[17] = 0;
  s[18] = 0;
  s[19] = 1;
  s[20] = 0;
  s[21] = 1;
  s[22] = 0;
  s[23] = 0;
  s[24] = 0;
  s[25] = 0;
  s[26] = 1;
  s[27] = 0;
  s[28] = 0;
  s[29] = 0;
  s[30] = 1;
  s[31] = 0;
  s[32] = 0;
  s[33] = 0;
  s[34] = 0;
  s[35] = 1;
  s[36] = 1;
  s[37] = 0;
  s[38] = 0;
  s[39] = 0;
  s[40] = 0;
  s[41] = 0;
  s[42] = 1;
  s[43] = 0;
  s[44] = 0;
  s[45] = 1;
  s[46] = 0;
  s[47] = 0;
  s[48] = 0;
  s[49] = 0;
  s[50] = 0;
  s[51] = 1;
  s[52] = 0;
  s[53] = 1;
  s[54] = 0;
  s[55] = 0;
  s[56] = 1;
  s[57] = 0;
  s[58] = 0;
  s[59] = 0;
  s[60] = 0;
  s[61] = 0;
  s[62] = 1;
  s[63] = 0;
  s[64] = 0;
  s[65] = 0;
  s[66] = 0;
  s[67] = 1;
  s[68] = 0;
  s[69] = 0;
  s[70] = 0;
  s[71] = 1;
  s[72] = 0;
  s[73] = 0;
  s[74] = 1;
  s[75] = 0;
  s[76] = 0;
  s[77] = 0;
  s[78] = 1;
  s[79] = 0;
  s[80] = 0;
  s[81] = 0;
  s[82] = 0;
  s[83] = 1;
  s[84] = 1;
  s[85] = 0;
  s[86] = 0;
  s[87] = 0;
  s[88] = 0;
  s[89] = 0;
  s[90] = 1;
  s[91] = 0;
  s[92] = 1;
  s[93] = 0;
  s[94] = 0;
  s[95] = 0;
  s[96] = 0;
  s[97] = 0;
  s[98] = 1;
  s[99] = 0;
  s[100] = 0;
  s[101] = 1;
  s[102] = 0;
  s[103] = 0;
  s[104] = 1;
  s[105] = 0;
  s[106] = 0;
  s[107] = 0;
  s[108] = 1;
  s[109] = 0;
  s[110] = 0;
  s[111] = 0;
  s[112] = 0;
  s[113] = 1;
  s[114] = 0;
  s[115] = 0;
  s[116] = 0;
  s[117] = 0;
  s[118] = 1;
  s[119] = 0;
  s[120] = 0;
  s[121] = 0;
  s[122] = 1;
  s[123] = 0;
  s[124] = 1;
  s[125] = 0;
  s[126] = 0;
  s[127] = 0;
  s[128] = 0;
  s[129] = 1;
  s[130] = 0;
  s[131] = 0;
  s[132] = 0;
  s[133] = 0;
  s[134] = 0;
  s[135] = 1;
  s[136] = 0;
  s[137] = 0;
  s[138] = 1;
  s[139] = 0;
  s[140] = 1;
  s[141] = 0;
  s[142] = 0;
  s[143] = 0;
  s[144] = 1;
  s[145] = 0;
  s[146] = 0;
  s[147] = 0;
  s[148] = 1;
  s[149] = 0;
  s[150] = 0;
  s[151] = 0;
  s[152] = 1;
  s[153] = 0;
  s[154] = 0;
  s[155] = 0;
  s[156] = 0;
  s[157] = 0;
  s[158] = 0;
  s[159] = 1;
  s[160] = 0;
  s[161] = 1;
  s[162] = 0;
  s[163] = 0;
  s[164] = 0;
  s[165] = 1;
  s[166] = 0;
  s[167] = 0;
  s[168] = 0;
  s[169] = 0;
  s[170] = 1;
  s[171] = 0;
  s[172] = 0;
  s[173] = 0;
  s[174] = 0;
  s[175] = 1;
  s[176] = 0;
  s[177] = 0;
  s[178] = 1;
  s[179] = 0;
  s[180] = 0;
  s[181] = 0;
  s[182] = 0;
  s[183] = 1;
  s[184] = 0;
  s[185] = 0;
  s[186] = 1;
  s[187] = 0;
  s[188] = 0;
  s[189] = 0;
  s[190] = 0;
  s[191] = 1;
  s[192] = 0;
  s[193] = 1;
  s[194] = 0;
  s[195] = 0;
  s[196] = 0;
  s[197] = 0;
  s[198] = 1;
  s[199] = 0;
  s[200] = 0;
  s[201] = 0;
  s[202] = 1;
  s[203] = 0;
  s[204] = 0;
  s[205] = 1;
  s[206] = 0;
  s[207] = 0;
  s[208] = 0;
  s[209] = 0;
  s[210] = 1;
  s[211] = 0;
  s[212] = 0;
  s[213] = 0;
  s[214] = 1;
  s[215] = 0;
  s[216] = 0;
  s[217] = 0;
  s[218] = 0;
  s[219] = 1;
  s[220] = 0;
  s[221] = 0;
  s[222] = 0;
  s[223] = 1;
  s[224] = 0;
  s[225] = 1;
  s[226] = 0;
  s[227] = 0;
  s[228] = 0;
  s[229] = 0;
  s[230] = 1;
  s[231] = 0;
  s[232] = 1;
  s[233] = 0;
  s[234] = 0;
  s[235] = 0;
  s[236] = 0;
  s[237] = 0;
  s[238] = 0;
  s[239] = 1;
  s[240] = 0;
  s[241] = 0;
  s[242] = 0;
  s[243] = 1;
  s[244] = 0;
  s[245] = 1;
  s[246] = 0;
  s[247] = 0;
  s[248] = 0;
  s[249] = 1;
  s[250] = 0;
  s[251] = 0;
  s[252] = 0;
  s[253] = 0;
  s[254] = 1;
  s[255] = 0;
  s[256] = 0;
  s[257] = 0;
  s[258] = 0;
  s[259] = 1;
  s[260] = 0;
  s[261] = 1;
  s[262] = 0;
  s[263] = 0;
  s[264] = 0;
  s[265] = 1;
  s[266] = 0;
  s[267] = 0;
  s[268] = 0;
  s[269] = 1;
  s[270] = 0;
  s[271] = 0;
  s[272] = 0;
  s[273] = 0;
  s[274] = 1;
  s[275] = 0;
  s[276] = 1;
  s[277] = 0;
  s[278] = 0;
  s[279] = 0;
  s[280] = 0;
  s[281] = 0;
  s[282] = 1;
  s[283] = 0;
  s[284] = 0;
  s[285] = 1;
  s[286] = 0;
  s[287] = 0;
  s[288] = 0;
  s[289] = 1;
  s[290] = 0;
  s[291] = 0;
  s[292] = 1;
  s[293] = 0;
  s[294] = 0;
  s[295] = 0;
  s[296] = 1;
  s[297] = 1;
  s[298] = 1;
  s[299] = 1;
}
