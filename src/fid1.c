#include "fid.h"


static float * dd = NULL;
static long ddlen = 0;


/* 
     
     FID 1 model implementation with floats and no scaling

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

void strale_fid1f_column(int m, int n,
                         float ph,float pb, float pe,
                         float phh, float phb, float phe,
                         float pbh, float pbb, float pbe,
                         float peh, float peb, float pee)
{

  float * ee = NULL;    /* pointer to the E section of the column */
  float * hh = NULL;    /* pointer to the H section of the column */
  float * bb = NULL;    /* pointer to the B section of the column */

  float diagh, diagb, diage;   /* hold values (h,b,e) of diagonal cell */
  float toph, topb, tope;      /* hold values (h,b,e) of upper cell */
  float temph, tempb, tempe;   /* temp variables when overwriting a cell */
  float h0, b0, e0;

  long i, j;

  /* reallocate matrix if necessary */
  if (3*m > ddlen)
  {
    free(dd);
    dd = xmalloc(3*m*sizeof(float), FID_ALIGNMENT_SSE);
    ddlen = 3*m;
  }
 
  /* init pointers to elements */
  hh = dd;
  bb = hh + m;
  ee = bb + m;
  
  b0 = ph*pbh + pb*pbb + pe*pbe;
  e0 = ph*peh + pb*peb + pe*pee;
  h0 = ph*phh + pb*phb + pe*phe;


  /* precompute cell (1,1) */
  *hh   = h0;
  *bb   = pbe*e0;
  *ee++ = peb*b0;
  
  /* precompute the rest of column 1 */
  for (i = 1; i < m; ++i)
  {
    *ee = *hh * peh + *bb * peb + *ee * pee;
    ++hh; ++bb; ++ee;
    *hh = phe*e0;

    e0 *= pee;
    *bb = pbe*e0;
  }
  ++hh; ++bb;
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
    diagb = topb; topb = diagb*pbb;;
    diage = tope = 0;

    
    /* iterate cells of a column */
    for (i = 0; i < m; ++i)
    {
      /* save curent cell before overwriitng, since it will be used as the
       * diagonal in the next round */
      temph = *hh; tempb = *bb; tempe = *ee;

      *bb = *hh*pbh   + *bb*pbb   + *ee*pbe;
      *hh = diagh*phh + diagb*phb + diage*phe;
      *ee = toph*peh  + topb*peb  + tope*pee;

      /* retrieve diagonal and top for next round */
      diagh = temph; diagb = tempb; diage = tempe;
      toph  = *hh++; topb  = *bb++; tope  = *ee++;
    }
  }
}

void strale_fid1f_printcol(void)
{
  long m = ddlen / 3;
  long i;

  float * hh = dd;
  float * bb = hh+m;
  float * ee = bb+m;

  for (i=0; i<m; ++i)
    printf("(%f, %f, %f)  ", *hh++, *bb++, *ee++); 
}

int main(int argc, char * argv[])
{
  strale_fid1f_column(20, 30,
                      0.5, 0.25, 0.25,
                      0.5, 0.2, 0.3,
                      0.28, 0.32, 0.4,
                      0.38, 0.32, 0.3);

  printf ("Format is (H,B,E) starting from row 1 until m of last column\n\n");
  strale_fid1f_printcol();

  return (EXIT_SUCCESS);
}
