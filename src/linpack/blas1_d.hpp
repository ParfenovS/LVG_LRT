#include <cstdlib>
#include <cmath>

#define sign_x(x) (std::signbit(x) ?  -1 : 1)

using namespace std;

template <typename real_type>
inline void daxpy(const size_t & n, real_type da, real_type dx[], real_type dy[]) {
//****************************************************************************80
//
//  Purpose:
//
//    DAXPY computes constant times a vector plus a vector.
//    On output, DY[*] has been replaced by DY[*] + DA * DX[*].
//
  #if defined(__clang__)
  #pragma clang loop vectorize(enable) interleave(enable)
  #elif defined(__INTEL_COMPILER) || defined(__ICL) || defined(__ICC) || defined(__ECC)
  #pragma ivdep
  #elif defined(__GNUC__)
  #pragma GCC ivdep
  #elif defined (_MSC_VER)
  #pragma loop( ivdep )
  #endif
  for (size_t i = 0; i < n; i++) {
    dy[i] = dy[i] + da * dx[i];
  }
  return;
}

//****************************************************************************80

template <typename real_type>
inline void dscal(const size_t & n, real_type sa, real_type x[]) {
//****************************************************************************80
//
//  Purpose:
//
//    DSCAL scales a vector by a constant.
//
  #if defined(__clang__)
  #pragma clang loop vectorize(enable) interleave(enable)
  #elif defined(__INTEL_COMPILER) || defined(__ICL) || defined(__ICC) || defined(__ECC)
  #pragma ivdep
  #elif defined(__GNUC__)
  #pragma GCC ivdep
  #elif defined (_MSC_VER)
  #pragma loop( ivdep )
  #endif
  for (size_t i = 0; i < n; i++) {
    x[i] = sa * x[i];
  }
  return;
}
//****************************************************************************80

template <typename real_type>
inline size_t idamax(const size_t & n, real_type dx[]) {
  real_type dmax, temp;
  size_t value = 0;

  dmax = fabs(dx[0]);
  for (size_t i = 1; i < n; i++) {
    temp = fabs(dx[i]);
    if (dmax < temp) {
      value = i;
      dmax = temp;
    }
  }
  return value;
}
//****************************************************************************80

template <typename real_type>
inline std::pair<size_t, size_t> idamax(const size_t & n, real_type dx[], const size_t & lda, const size_t & k) {
  real_type dmax, temp;
  std::pair<size_t, size_t> value = std::make_pair(k, k);

  dmax = fabs(dx[k + k * lda]);
  for (size_t i = k; i < n; i++) {
    for (size_t j = k; j < n; j++) {
      temp = fabs(dx[i + j * lda]);
      if (dmax < temp) {
        value = std::make_pair(i, j);
        dmax = temp;
      }
    }
  }
  return value;
}
//****************************************************************************80

template <typename real_type>
real_type dnrm2 ( const size_t & n, real_type x[])

//****************************************************************************80
//
//  Purpose:
//
//    DNRM2 returns the euclidean norm of a vector.
//
//  Discussion:
//
//     DNRM2 ( X ) = sqrt ( X' * X )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Charles Lawson, Richard Hanson, 
//    David Kincaid, Fred Krogh.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double X[*], the vector whose norm is to be computed.
//
//    Output, double DNRM2, the Euclidean norm of X.
//
{
  real_type absxi;
  real_type scale = 0.0;
  real_type ssq = 1.0;
  size_t ix = 0;

  for ( size_t i = 0; i < n; i++ )
  {
    if ( x[ix] != 0.0 )
    {
      absxi = fabs ( x[ix] );
      if ( scale < absxi )
      {
        ssq = 1.0 + ssq * ( scale / absxi ) * ( scale / absxi );
        scale = absxi;
      }
      else
      {
        ssq = ssq + ( absxi / scale ) * ( absxi / scale );
      }
    }
    ix = ix + 1;
  }

  return scale * sqrt ( ssq );
}
//****************************************************************************80

template <typename real_type>
real_type ddot ( const size_t & n, real_type dx[], real_type dy[])

//****************************************************************************80
//
//  Purpose:
//
//    DDOT forms the dot product of two vectors.
//
//  Discussion:
//
//    This routine uses unrolled loops for increments equal to one.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Charles Lawson, Richard Hanson, 
//    David Kincaid, Fred Krogh.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double DX[*], the first vector.
//
//    Input, double DY[*], the second vector.
//
//    Output, double DDOT, the sum of the product of the corresponding
//    entries of DX and DY.
//
{
  real_type dtemp = 0.0;
  const size_t m = n % 5;

  for ( size_t i = 0; i < m; i++ )
  {
    dtemp = dtemp + dx[i] * dy[i];
  }

  for ( size_t i = m; i < n; i = i + 5 )
  {
    dtemp = dtemp + dx[i  ] * dy[i  ]
                  + dx[i+1] * dy[i+1]
                  + dx[i+2] * dy[i+2]
                  + dx[i+3] * dy[i+3]
                  + dx[i+4] * dy[i+4];
  }

  return dtemp;
}
//****************************************************************************80

template <typename real_type>
void drotg ( real_type *sa, real_type *sb, real_type *c, real_type *s )

//****************************************************************************80
//
//  Purpose:
//
//    DROTG constructs a Givens plane rotation.
//
//  Discussion:
//
//    Given values A and B, this routine computes
//
//    SIGMA = sign ( A ) if abs ( A ) >  abs ( B )
//          = sign ( B ) if abs ( A ) <= abs ( B );
//
//    R     = SIGMA * ( A * A + B * B );
//
//    C = A / R if R is not 0
//      = 1     if R is 0;
//
//    S = B / R if R is not 0,
//        0     if R is 0.
//
//    The computed numbers then satisfy the equation
//
//    (  C  S ) ( A ) = ( R )
//    ( -S  C ) ( B ) = ( 0 )
//
//    The routine also computes
//
//    Z = S     if abs ( A ) > abs ( B ),
//      = 1 / C if abs ( A ) <= abs ( B ) and C is not 0,
//      = 1     if C is 0.
//
//    The single value Z encodes C and S, and hence the rotation:
//
//    If Z = 1, set C = 0 and S = 1;
//    If abs ( Z ) < 1, set C = sqrt ( 1 - Z * Z ) and S = Z;
//    if abs ( Z ) > 1, set C = 1/ Z and S = sqrt ( 1 - C * C );
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 May 2006
//
//  Author:
//
//    Original FORTRAN77 version by Charles Lawson, Richard Hanson, 
//    David Kincaid, Fred Krogh.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input/output, double *SA, *SB,  On input, SA and SB are the values
//    A and B.  On output, SA is overwritten with R, and SB is
//    overwritten with Z.
//
//    Output, double *C, *S, the cosine and sine of the Givens rotation.
//
{
  real_type r;
  real_type roe;
  real_type scale;
  real_type z;

  if ( fabs ( *sb ) < fabs ( *sa ) )
  {
    roe = *sa;
  }
  else
  {
    roe = *sb;
  }

  scale = fabs ( *sa ) + fabs ( *sb );

  if ( scale == 0.0 )
  {
    *c = 1.0;
    *s = 0.0;
    r = 0.0;
  }
  else
  {
    r = scale * sqrt ( ( *sa / scale ) * ( *sa / scale )
                     + ( *sb / scale ) * ( *sb / scale ) );
    r = sign_x(roe) * r;
    *c = *sa / r;
    *s = *sb / r;
  }

  if ( 0.0 < fabs ( *c ) && fabs ( *c ) <= *s )
  {
    z = 1.0 / *c;
  }
  else
  {
    z = *s;
  }

  *sa = r;
  *sb = z;

  return;
}
