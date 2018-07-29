#include <math.h>

  //  Given a three-dimensional real array data[1..nn1][1..nn2][1..nn3] (where
  //  nn1 = 1 for the case of a logically two-dimensional array), this routine
  //  returns (for isign=1) the complex fast Fourier transform as two complex
  //  arrays: On output, data contains the zero and positive frequency values
  //  of the third frequency component, while speq[1..nn1][1..2*nn2]contains
  //  the Nyquist critical frequency values of the third frequency
  //  component. First (and second) frequency components are stored for zero,
  //  positive, and negative frequencies, in standard wrap- around order. See
  //  text for description of how complex values are arranged. For isign=-1,
  //  the inverse transform (times nn1*nn2*nn3/2 as a constant multiplicative
  //  factor) is performed, with output data (viewed as a real array) deriving
  //  from input data (viewed as complex) and speq. For inverse transforms on
  //  data not generated  rst by a forward transform, make sure the complex
  //  input data array satis es property (12.5.2). The dimensions nn1, nn2,
  //  nn3 must always be integer powers of 2.  

void rlft3 ( float ***data, float **speq,
	     unsigned long nn1, unsigned long nn2, unsigned long nn3,
	     int isign )
{ 
  void fourn(float data[], unsigned long nn[], int ndim, int isign); 
  void nrerror(char error_text[]); 

  unsigned long i1, i2, i3, j1, j2, j3, nn[4], ii3; 
  double wi, wr;
  float h1r, h1i, h2r, h2i;

  if ( 1+&data[nn1][nn2][nn3]-&data[1][1][1] != nn1*nn2*nn3 )
    nrerror("rlft3: problem with dimensions or contiguity of data array\n");

  float c1     = 0.5;
  float c2     = -0.5*isign;
  double theta = isign*(6.28318530717959/nn3);
  double wtemp = sin(0.5*theta);
  double wpr   = -2.0*wtemp*wtemp;
  double wpi   = sin(theta);

  nn[1] = nn1;
  nn[2] = nn2;
  nn[3] = nn3 >> 1;

  if (isign == 1)  // Case of forward transform.
    {
      // Here is where most all of the compute time is spent.

      fourn ( &data[1][1][1]-1, nn, 3, isign );

      //Extend data periodically into speq.

      for ( i1 = 1; i1 <= nn1; i1++ )
	for ( i2 = 1, j2 = 0; i2 <= nn2; i2++ )
	  {
	    speq[i1][++j2]=data[i1][i2][1];
	    speq[i1][++j2]=data[i1][i2][2];
	  }
    }
  
  for ( i1 = 1; i1 <= nn1; i1++ )
    {
      // Zero frequency is its own reflection, otherwise locate corresponding
      // negative frequency in wrap-around order.

      j1 = ( i1 != 1 ? nn1-i1+2 : 1 );

      // Initialize trigonometric recurrence.

      wr = 1.0;
      wi = 0.0;

      for ( ii3 = 1, i3 = 1; i3 <= (nn3>>2)+1; i3++, ii3+=2 )
	{
	  for ( i2 = 1; i2 <= nn2; i2++ )
	    {
	      if ( i3 == 1 )
		{ //Equation (12.3.5).
		  j2 = (i2 != 1 ? ((nn2-i2)<<1)+3 : 1);
		  h1r = c1*(data[i1][i2][1]+speq[j1][j2]);
		  h1i = c1*(data[i1][i2][2]-speq[j1][j2+1]);
		  h2i = c2*(data[i1][i2][1]-speq[j1][j2]);
		  h2r =  -c2*(data[i1][i2][2]+speq[j1][j2+1]);
		  data[i1][i2][1] = h1r+h2r;
		  data[i1][i2][2] = h1i+h2i;
		  speq[j1][j2] = h1r-h2r;
		  speq[j1][j2+1] = h2i-h1i;
		}
	      else
		{
		  j2 = (i2 != 1 ? nn2-i2+2 : 1);
		  j3 = nn3+3-(i3<<1);
		  h1r = c1*(data[i1][i2][ii3]+data[j1][j2][j3]);
		  h1i = c1*(data[i1][i2][ii3+1]-data[j1][j2][j3+1]);
		  h2i = c2*(data[i1][i2][ii3]-data[j1][j2][j3]);
		  h2r =  -c2*(data[i1][i2][ii3+1]+data[j1][j2][j3+1]);
		  data[i1][i2][ii3] = h1r+wr*h2r-wi*h2i;
		  data[i1][i2][ii3+1] = h1i+wr*h2i+wi*h2r;
		  data[j1][j2][j3] = h1r-wr*h2r+wi*h2i;
		  data[j1][j2][j3+1] =  -h1i+wr*h2i+wi*h2r;
		}
	    }
	  
	  wr = ( wtemp = wr ) * wpr - wi * wpi + wr; // Do the recurrence.
	  wi = wi * wpr + wtemp * wpi + wi;
      }
    }

  if (isign == -1) // Case of reverse transform.
    fourn( &data[1][1][1]-1, nn, 3, isign );
} 
