#include <cmath>
#include <cstdlib>
#include "blas1_d.hpp"

using namespace std;

template <typename real_type>
size_t dgefa ( real_type a[], const size_t & lda, const size_t & n, size_t ipvt[] )

//****************************************************************************80
//
//  Purpose:
//
//    DGEFA factors a real general matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch, 
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input/output, real_type A[LDA*N].
//    On intput, the matrix to be factored.
//    On output, an upper triangular matrix and the multipliers used to obtain
//    it.  The factorization can be written A=L*U, where L is a product of
//    permutation and unit lower triangular matrices, and U is upper triangular.
//
//    Input, int LDA, the leading dimension of A.
//
//    Input, int N, the order of the matrix A.
//
//    Output, int IPVT[N], the pivot indices.
//
//    Output, int DGEFA, singularity indicator.
//    0, normal value.
//    K, if U(K,K) == 0.  This is not an error condition for this subroutine,
//    but it does indicate that DGESL or DGEDI will divide by zero if called.
//    Use RCOND in DGECO for a reliable indication of singularity.
//
{
  size_t info;
  size_t j;
  size_t k;
  size_t l;
  real_type t;
//
//  Gaussian elimination with partial pivoting.
//
  info = 0;

  for ( k = 0; k <= n-2; k++ )
  {
//
//  Find L = pivot index.
//
    l = idamax ( n-k, a+k+k*lda ) + k;
    ipvt[k] = l;
//
//  Zero pivot implies this column already triangularized.
//
    if ( a[l+k*lda] == 0.0 )
    {
      info = k + 1;
      continue;
    }
//
//  Interchange if necessary.
//
    if ( l != k )
    {
      t = a[l+k*lda];
      a[l+k*lda] = a[k+k*lda];
      a[k+k*lda] = t;
    }
//
//  Compute multipliers.
//
    t = -1.0 / a[k+k*lda];

    dscal ( n-k-1, t, a+k+1+k*lda );
//
//  Row elimination with column indexing.
//
    for ( j = k+1; j < n; j++ )
    {
      t = a[l+j*lda];
      if ( l != k )
      {
        a[l+j*lda] = a[k+j*lda];
        a[k+j*lda] = t;
      }
      daxpy ( n-k-1, t, a+k+1+k*lda, a+k+1+j*lda );
    }

  }

  ipvt[n-1] = n-1;

  if ( a[n-1+(n-1)*lda] == 0.0 )
  {
    info = n;
  }

  return info;
}
//****************************************************************************80

template <typename real_type>
void dgesl ( real_type a[], const size_t & lda, const size_t & n, size_t ipvt[], real_type b[] )

//****************************************************************************80
//
//  Purpose:
//
//    DGESL solves a real general linear system A * X = B.
//
//  Discussion:
//
//    DGESL can solve either of the systems A * X = B or A' * X = B.
//
//    The system matrix must have been factored by DGECO or DGEFA.
//
//    A division by zero will occur if the input factor contains a
//    zero on the diagonal.  Technically this indicates singularity
//    but it is often caused by improper arguments or improper
//    setting of LDA.  It will not occur if the subroutines are
//    called correctly and if DGECO has set 0.0 < RCOND
//    or DGEFA has set INFO == 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch, 
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input, real_type A[LDA*N], the output from DGECO or DGEFA.
//
//    Input, int LDA, the leading dimension of A.
//
//    Input, int N, the order of the matrix A.
//
//    Input, int IPVT[N], the pivot vector from DGECO or DGEFA.
//
//    Input/output, real_type B[N].
//    On input, the right hand side vector.
//    On output, the solution vector.
//
{
  size_t k;
  size_t l;
  real_type t;
//
//  Solve A * X = B.
//
  for ( k = 0; k <= n-2; k++ )
  {
    l = ipvt[k];
    t = b[l];

    if ( l != k )
    {
      b[l] = b[k];
      b[k] = t;
    }

    daxpy ( n-k-1, t, a+k+1+k*lda, b+k+1 );

  }

  for ( k = n; k-- > 0; )
  {
    b[k] = b[k] / a[k+k*lda];
    t = -b[k];
    daxpy ( k, t, a+0+k*lda, b );
  }

  return;
}
