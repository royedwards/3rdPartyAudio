/* ------------------------------------------------------------------------- *\
   WAVE.C :

   This package provides 3 routines for generating rectangular-like,
   triangular-like and sine-like waves including specific features.

   by Christophe Schlick (10 September 1993)

   "Wave Generators for Procedural Techniques in Computer Graphics"
   in Graphics Gems V (edited by A. Paeth), Academic Press
\* ------------------------------------------------------------------------- */

#include <math.h>
#include "wave.h"

/*
** Macro functions
*/

#define ABS(a)          ((a) < 0 ? -(a) : (a))
#define FLOOR(a)        ((a) < 0 ? (int) ((a)-1.0) : (int) (a))
#define MAX(a,b)        ((a) > (b) ? (a) : (b))
#define MIN(a,b)        ((a) < (b) ? (a) : (b))

/*
** rnd : Random function (adapted from Greg Ward in Graphics Gems II)
*/

#define FRND(a)         rnd (17*(a))
#define ARND(a)         rnd (97*(a))

static double rnd (register long s)
{
  s = s << 13 ^ s;
  return (((s*(s*s*15731+789221)+1376312589) & 0X7FFFFFFF) / 2147483648.0);
}

/*
** Rwave : Rectangular-like monodimensional wave
**
** Input : t = Wave function parameter
**         s = Shape parameter (in [-1,1])
**         f = Frequency variance (in [0,1])
**         a = Amplitude variance (in [0,1])
*/

double Rwave (register double t, double s, double f, double a)
{
  register int    i, j;
  register double b;

  i = j = FLOOR (t); t -= i; j++;

  if (f) {
    a = (FRND (i) - 0.5) * f;
    b = (FRND (j) - 0.5) * f + 1.0;
    t = (t-a) / (b-a);
  }

  if (i & 1) {i++; j--; t = 1.0-t;}
  t = (s < 0.0) ? (t+s*t) / (1.0+s*t) : (s > 0.0) ? t / (1.0-s+s*t) : t;
  t = t < 0.5 ? 0.0 : 1.0;

  if (a) {
    a = ARND (i) * a * 0.5;
    b = ARND (j) * a * 0.5;
    t = a + t * (1.0-a-b);
  }
  return (t);
}

/*
** Twave : Triangular-like monodimensional wave
**
** Input : t = Wave function parameter
**         s = Shape parameter (in [-1,1])
**         f = Frequency variance (in [0,1])
**         a = Amplitude variance (in [0,1])
*/

double Twave (register double t, double s, double f, double a) {
  register int    i, j;
  register double b;

  i = j = FLOOR (t); t -= i; j++;

  if (f) {
    a = (FRND (i) - 0.5) * f;
    b = (FRND (j) - 0.5) * f + 1.0;
    if (t < a) {
      i--; j--; t++; a++;
      b = a; a = (FRND (i) - 0.5) * f;
    } else if (t > b) {
      i++; j++; t--; b--;
      a = b; b = (FRND (j) - 0.5) * f + 1.0;
    }
    t = (t-a) / (b-a);
  }

  if (i & 1) {i++; j--; t = 1.0-t;}
  t = (s < 0.0) ? (t+s*t) / (1.0+s*t) : (s > 0.0) ? t / (1.0-s+s*t) : t;

  if (a) {
    a = ARND(i) * a * 0.5;
    b = ARND(j) * a * 0.5;
    t = a + t * (1.0-a-b);
  }
  return (t);
}

/*
** Swave : sinusoidal-like monodimensional wave
**
** Input : t = Wave function parameter
**         s = Shape parameter (in [-1,1])
**         f = Frequency variance (in [0,1])
**         a = Amplitude variance (in [0,1])
*/

double Swave (register double t, double s, double f, double a)
{
  register int    i, j;
  register double b;

  i = j = FLOOR (t); t -= i; j++;

  if (f) {
    a = (FRND (i) - 0.5) * f;
    b = (FRND (j) - 0.5) * f + 1.0;
    if (t < a) {
      i--; j--; t++; a++;
      b = a; a = (FRND (i) - 0.5) * f;
    } else if (t > b) {
      i++; j++; t--; b--;
      a = b; b = (FRND (j) - 0.5) * f + 1.0;
    }
    t = (t-a) / (b-a);
  }

  if (i & 1) {i++; j--; t = 1.0-t;}
  t = (s < 0.0) ? (t+s*t) / (1.0+s*t) : (s > 0.0) ? t / (1.0-s+s*t) : t;
  t *= t * (3.0-t-t);

  if (a) {
    a = ARND(i) * a * 0.5;
    b = ARND(j) * a * 0.5;
    t = a + t * (1.0-a-b);
  }
  return (t);
}

/* ------------------------------------------------------------------------- */
