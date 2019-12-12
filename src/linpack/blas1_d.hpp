#include <cstdlib>
#include <cmath>

using namespace std;

template <typename real_type>
inline void daxpy ( const size_t & n, real_type da, real_type dx[], real_type dy[] )

//****************************************************************************80
//
//  Purpose:
//
//    DAXPY computes constant times a vector plus a vector.
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
//    Input, int N, the number of elements in DX and DY.
//
//    Input, double DA, the multiplier of DX.
//
//    Input, double DX[*], the first vector.
//
//    Input/output, double DY[*], the second vector.
//    On output, DY[*] has been replaced by DY[*] + DA * DX[*].
//
{
  size_t i;
  size_t m;

  if ( da == 0.0 )
  {
    return;
  }

  m = n % 4;

  for ( i = 0; i < m; i++ )
  {
    dy[i] = dy[i] + da * dx[i];
  }

  for ( i = m; i < n; i = i + 4 )
  {
    dy[i  ] = dy[i  ] + da * dx[i  ];
    dy[i+1] = dy[i+1] + da * dx[i+1];
    dy[i+2] = dy[i+2] + da * dx[i+2];
    dy[i+3] = dy[i+3] + da * dx[i+3];
  }

  return;
}

//****************************************************************************80

template <typename real_type>
inline void dscal ( const size_t & n, real_type sa, real_type x[] )

//****************************************************************************80
//
//  Purpose:
//
//    DSCAL scales a vector by a constant.
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
//    Input, double SA, the multiplier.
//
//    Input/output, double X[*], the vector to be scaled.
//
{
  size_t i;
  const size_t m = n % 5;

  for ( i = 0; i < m; i++ )
  {
    x[i] = sa * x[i];
  }

  for ( i = m; i < n; i = i + 5 )
  {
    x[i]   = sa * x[i];
    x[i+1] = sa * x[i+1];
    x[i+2] = sa * x[i+2];
    x[i+3] = sa * x[i+3];
    x[i+4] = sa * x[i+4];
  }

  return;
}
//****************************************************************************80

template <typename real_type>
inline size_t idamax ( const size_t & n, real_type dx[] )

//****************************************************************************80
//
//  Purpose:
//
//    IDAMAX finds the index of the vector element of maximum absolute value.
//
//  Discussion:
//
//    WARNING: This index is a 1-based index, not a 0-based index!
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
//    Input, double X[*], the vector to be examined.
//
//    Output, int IDAMAX, the index of the element of maximum
//    absolute value.
//
{
  real_type dmax;
  size_t i;
  size_t value;

  value = 0;

  if ( n == 1 )
  {
    return value;
  }

  dmax = fabs ( dx[0] );

  for ( i = 1; i < n; i++ )
  {
    if ( dmax < fabs ( dx[i] ) )
    {
      value = i;
      dmax = fabs ( dx[i] );
    }
  }

  return value;
}
