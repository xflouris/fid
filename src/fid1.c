#include "fid.h"

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

double strale_fid1d_column(int m, int n,
                           double ph, double pb, double pe)
{

  double * ee = NULL;    /* pointer to the E section of the column */
  double * hh = NULL;    /* pointer to the H section of the column */
  double * bb = NULL;    /* pointer to the B section of the column */

  double diagh, diagb, diage;   /* hold values (h,b,e) of diagonal cell */
  double toph, topb, tope;      /* hold values (h,b,e) of upper cell */
  double temph, tempb, tempe;   /* temp variables when overwriting a cell */
  double h0, b0, e0;

  long i, j;

  /* reallocate matrix if necessary */
  if (3*m > ddlen)
  {
    free(dd);
    dd = xmalloc(3*m*sizeof(double), FID_ALIGNMENT_SSE);
    ddlen = 3*m;
  }
 
  /* init pointers to elements */
  hh = dd;
  bb = hh + m;
  ee = bb + m;
  
  b0 = ph*phb + pb*pbb + pe*peb;
  e0 = ph*phe + pb*pbe + pe*pee;
  h0 = ph*phh + pb*pbh + pe*peh;


  /* precompute cell (1,1) */
  hh[0] = h0;
  bb[0] = peb*e0;
  ee[0] = pbe*b0;

  /* precompute the rest of column 1 */
  for (i = 1; i < m; ++i)
  {
    ee[i] = hh[i-1]*phe + bb[i-1]*pbe + ee[i-1]*pee;
    hh[i] = peh*e0;
    e0 *= pee;
    bb[i] = peb*e0;
  }
  topb = b0;

  /* iterate through columns starting from the second */
  for (j = 1; j < n; ++j)
  {
    /* point each entry to the beginning (first element) */
    hh = dd;
    bb = hh + m;
    ee = bb + m;

    /* compute diagonal and top elements from init row */
    diagh = toph = 0;
    diagb = topb; topb = diagb*pbb;
    diage = tope = 0;

    
    /* iterate cells of a column */
    for (i = 0; i < m; ++i)
    {
      /* save curent cell before overwriitng, since it will be used as the
       * diagonal in the next round */
      temph = hh[i]; tempb = bb[i]; tempe = ee[i];

      bb[i] = hh[i]*phb + bb[i]*pbb + ee[i]*peb;
      hh[i] = diagh*phh + diagb*pbh + diage*peh;
      ee[i] = toph*phe  + topb*pbe  + tope*pee;

      /* retrieve diagonal and top for next round */
      diagh = temph; diagb = tempb; diage = tempe;
      toph  = hh[i]; topb  = bb[i]; tope  = ee[i];
    }
  }

  return (dd[m-1]*phh + dd[2*m-1]*pbh + dd[3*m-1]*peh);
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

void strale_fid1d_printcol(void)
{
  long m = ddlen / 3;
  long i;

  double * hh = dd;
  double * bb = hh+m;
  double * ee = bb+m;

  for (i=0; i<m; ++i)
    printf("(%f, %f, %f)  ", *hh++, *bb++, *ee++); 
  printf ("\n");
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

  h = strale_fid1d_column(20,
                          13,
                          1, 
                          0, 
                          0);

  printf ("-> H: %f\n", h);
  
  printf ("Format is (H,B,E) starting from row 1 until m of last column\n\n");
  strale_fid1d_printcol();

  return (EXIT_SUCCESS);
}
