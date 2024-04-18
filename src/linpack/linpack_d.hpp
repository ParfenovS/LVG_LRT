#include <cmath>
#include <cstdlib>
#include "blas1_d.hpp"

using namespace std;

template <typename real_type>
size_t dgefa(real_type a[], const size_t &lda, const size_t &n, size_t ipvt[]) {
  size_t info;
  size_t j;
  size_t k;
  size_t l;
  real_type t;
  //
  //  Gaussian elimination with partial pivoting.
  //
  info = 0;

  for (k = 0; k < n - 1; k++) {
    //
    //  Find L = pivot index.
    //
    l = idamax(n - k, a + k + k * lda) + k;
    ipvt[k] = l;
    //
    //  Zero pivot implies this column already triangularized.
    //
    if (a[l + k * lda] == 0.0) {
      info = k + 1;
      continue;
    }
    //
    //  Interchange
    //
    t = a[l + k * lda];
    a[l + k * lda] = a[k + k * lda];
    a[k + k * lda] = t;
    //
    //  Compute multipliers.
    //
    t = - 1.0 / a[k + k * lda];

    dscal(n - k - 1, t, a + k + 1 + k * lda);
    //
    //  Row elimination with column indexing.
    //
    if (l != k) {
      for (j = k + 1; j < n; j++) {
        t = a[l + j * lda];
        a[l + j * lda] = a[k + j * lda];
        a[k + j * lda] = t;
        daxpy(n - k - 1, t, a + k + 1 + k * lda, a + k + 1 + j * lda);
      }
    } else {
      for (j = k + 1; j < n; j++) {
        t = a[l + j * lda];
        daxpy(n - k - 1, t, a + k + 1 + k * lda, a + k + 1 + j * lda);
      }
    }
  }

  ipvt[n - 1] = n - 1;
  if (a[n - 1 + (n - 1) * lda] == 0.0) info = n;

  return info;
}

template <typename real_type>
void dgesl(real_type a[], const size_t &lda, const size_t &n, size_t ipvt[], real_type b[]) {
  size_t k;
  size_t l;
  real_type t;
  //
  //  Solve A * X = B.
  //
  for (k = 0; k < n - 1; k++) {
    l = ipvt[k];
    t = b[l];
    if (l != k) {
      b[l] = b[k];
      b[k] = t;
    }
    daxpy(n - k - 1, t, a + k + 1 + k * lda, b + k + 1);
  }
  for (k = n; k-- > 0;) {
    b[k] = b[k] / a[k + k * lda];
    t = -b[k];
    daxpy(k, t, a + 0 + k * lda, b);
  }
  return;
}

template <typename real_type>
size_t dgefa(real_type a[], const size_t &lda, const size_t &n) {
  size_t info;
  size_t j;
  size_t k;
  real_type t;
  //
  //  Gaussian elimination without pivoting.
  //
  info = 0;

  for (k = 0; k < n - 1; k++) {
    //
    //  Compute multipliers.
    //
    t = - 1.0 / a[k + k * lda];

    dscal(n - k - 1, t, a + k + 1 + k * lda);
    //
    //  Row elimination
    //
    for (j = k + 1; j < n; j++) {
      t = a[k + j * lda];
      daxpy(n - k - 1, t, a + k + 1 + k * lda, a + k + 1 + j * lda);
    }
  }

  if (a[n - 1 + (n - 1) * lda] == 0.0) info = n;

  return info;
}

template <typename real_type>
void dgesl(real_type a[], const size_t &lda, const size_t &n, real_type b[]) {
  size_t k;
  real_type t;
  //
  //  Solve A * X = B.
  //
  for (k = 0; k < n - 1; k++) {
    t = b[k];
    daxpy(n - k - 1, t, a + k + 1 + k * lda, b + k + 1);
  }
  for (k = n; k-- > 0;) {
    b[k] = b[k] / a[k + k * lda];
    t = -b[k];
    daxpy(k, t, a + 0 + k * lda, b);
  }
  return;
}

template <typename real_type>
size_t dsvdc(real_type a[], const size_t &lda, const size_t &n, real_type s[], real_type e[], real_type work[])

//****************************************************************************80
//
//  Purpose:
//
//    DSVDC computes the singular value decomposition of a real rectangular matrix.
//
//  Discussion:
//
//    This routine reduces an N by N matrix A to diagonal form by orthogonal
//    transformations U and V.  The diagonal elements S(I) are the singular
//    values of A.  The columns of U are the corresponding left singular
//    vectors, and the columns of V the right singular vectors.
//
//    The form of the singular value decomposition is then
//
//      A(NxN) = U(NxN) * S(NxN) * V(NxN)'
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 May 2007
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch,
//    Pete Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input/output, double A[LDA*N].  On input, the M by N matrix whose
//    singular value decomposition is to be computed.  On output, the matrix
//    has been destroyed.  Depending on the user's requests, the matrix may
//    contain other useful information.
//
//    Input, int LDA, the leading dimension of the array A.
//    LDA must be at least N.
//
//    Input, int N, the number of columns of the matrix A.
//
//    Output, double S[MM], where MM = N.  The first
//    min(M,N) entries of S contain the singular values of A arranged in
//    descending order of magnitude.
//
//    Output, double E[MM], where MM = N, ordinarily contains zeros.
//    However see the discussion of INFO for exceptions.
//
//    Workspace, double WORK[N].
//
//    Output, int *DSVDC, status indicator INFO.
//    The singular values (and their corresponding singular vectors)
//    S(*INFO+1), S(*INFO+2),...,S(MN) are correct.  Here MN = N.
//    Thus if *INFO is 0, all the singular values and their vectors are
//    correct.  In any event, the matrix B = U' * A * V is the bidiagonal
//    matrix with the elements of S on its diagonal and the elements of E on
//    its superdiagonal.  Thus the singular values of A and B are the same.
//
{
  real_type b;
  real_type c;
  real_type cs;
  real_type el;
  real_type emm1;
  real_type f;
  real_type g;
  size_t l = 0;
  size_t info;
  size_t iter;
  size_t ls = 0;
  size_t k;
  int kase;
  const size_t maxit = 30;
  size_t mm;
  size_t mm1;
  size_t mn;
  size_t nctp1;
  size_t nrtp1;
  real_type scale;
  real_type shift;
  real_type sl;
  real_type sm;
  real_type smm1;
  real_type sn;
  real_type t;
  real_type t1;
  real_type test;
  real_type ztest;
  //
  //  Determine what is to be computed.
  //
  info = 0;
  //
  //  Reduce A to bidiagonal form, storing the diagonal elements
  //  in S and the super-diagonal elements in E.
  //
  const size_t nct = n - 1;
  const size_t nrt = n - 2;
  const size_t lu = nct;

  for (l = 1; l <= lu; l++)
  {
    //
    //  Compute the transformation for the L-th column and
    //  place the L-th diagonal in S(L).
    //
    if (l <= nct)
    {
      s[l - 1] = dnrm2(n - l + 1, a + l - 1 + (l - 1) * lda);

      if (s[l - 1] != 0.0)
      {
        if (a[l - 1 + (l - 1) * lda] != 0.0)
        {
          s[l - 1] = sign_x(a[l - 1 + (l - 1) * lda]) * fabs(s[l - 1]);
        }
        dscal(n - l + 1, 1.0 / s[l - 1], a + l - 1 + (l - 1) * lda);
        a[l - 1 + (l - 1) * lda] = 1.0 + a[l - 1 + (l - 1) * lda];
      }
      s[l - 1] = -s[l - 1];
    }

    for (size_t j = l + 1; j <= n; j++)
    {
      //
      //  Apply the transformation.
      //
      if (l <= nct && s[l - 1] != 0.0)
      {
        t = -ddot(n - l + 1, a + l - 1 + (l - 1) * lda, a + l - 1 + (j - 1) * lda) / a[l - 1 + (l - 1) * lda];
        daxpy(n - l + 1, t, a + l - 1 + (l - 1) * lda, a + l - 1 + (j - 1) * lda);
      }
      //
      //  Place the L-th row of A into E for the
      //  subsequent calculation of the row transformation.
      //
      e[j - 1] = a[l - 1 + (j - 1) * lda];
    }

    if (l <= nrt)
    {
      //
      //  Compute the L-th row transformation and place the
      //  L-th superdiagonal in E(L).
      //
      e[l - 1] = dnrm2(n - l, e + l);

      if (e[l - 1] != 0.0)
      {
        if (e[l] != 0.0)
        {
          e[l - 1] = sign_x(e[l]) * fabs(e[l - 1]);
        }
        dscal(n - l, 1.0 / e[l - 1], e + l);
        e[l] = 1.0 + e[l];
      }

      e[l - 1] = -e[l - 1];
      //
      //  Apply the transformation.
      //
      if (l + 1 <= n && e[l - 1] != 0.0)
      {
        for (size_t j = l + 1; j <= n; j++)
        {
          work[j - 1] = 0.0;
        }

        for (size_t j = l + 1; j <= n; j++)
        {
          daxpy(n - l, e[j - 1], a + l + (j - 1) * lda, work + l);
        }

        for (size_t j = l + 1; j <= n; j++)
        {
          daxpy(n - l, -e[j - 1] / e[l], work + l, a + l + (j - 1) * lda);
        }
      }
    }
  }
  //
  //  Set up the final bidiagonal matrix of order MN.
  //
  mn = n;
  nctp1 = nct + 1;
  nrtp1 = nrt + 1;

  if (nct < n)
  {
    s[nctp1 - 1] = a[nctp1 - 1 + (nctp1 - 1) * lda];
  }

  if (nrtp1 < mn)
  {
    e[nrtp1 - 1] = a[nrtp1 - 1 + (mn - 1) * lda];
  }

  e[mn - 1] = 0.0;
  //
  //  Main iteration loop for the singular values.
  //
  mm = mn;
  iter = 0;

  while (0 < mn)
  {
    //
    //  If too many iterations have been performed, set flag and return.
    //
    if (maxit <= iter)
    {
      info = mn;
      return info;
    }
    //
    //  This section of the program inspects for
    //  negligible elements in the S and E arrays.
    //
    //  On completion the variables KASE and L are set as follows:
    //
    //  KASE = 1     if S(MN) and E(L-1) are negligible and L < MN
    //  KASE = 2     if S(L) is negligible and L < MN
    //  KASE = 3     if E(L-1) is negligible, L < MN, and
    //               S(L), ..., S(MN) are not negligible (QR step).
    //  KASE = 4     if E(MN-1) is negligible (convergence).
    //
    for (size_t ll = 1; ll <= mn; ll++)
    {
      l = mn - ll;

      if (l == 0)
      {
        break;
      }

      test = fabs(s[l - 1]) + fabs(s[l]);
      ztest = test + fabs(e[l - 1]);

      if (ztest == test)
      {
        e[l - 1] = 0.0;
        break;
      }
    }

    if (l == mn - 1)
    {
      kase = 4;
    }
    else
    {
      for (size_t lls = l + 1; lls <= mn + 1; lls++)
      {
        ls = mn - lls + l + 1;

        if (ls == l)
        {
          break;
        }

        test = 0.0;
        if (ls != mn)
        {
          test = test + fabs(e[ls - 1]);
        }

        if (ls != l + 1)
        {
          test = test + fabs(e[ls - 2]);
        }

        ztest = test + fabs(s[ls - 1]);

        if (ztest == test)
        {
          s[ls - 1] = 0.0;
          break;
        }
      }

      if (ls == l)
      {
        kase = 3;
      }
      else if (ls == mn)
      {
        kase = 1;
      }
      else
      {
        kase = 2;
        l = ls;
      }
    }

    l = l + 1;
    //
    //  Deflate negligible S(MN).
    //
    if (kase == 1)
    {
      mm1 = mn - 1;
      f = e[mn - 2];
      e[mn - 2] = 0.0;

      for (size_t kk = 1; kk <= mm1; kk++)
      {
        k = mm1 - kk + l;
        t1 = s[k - 1];
        drotg(&t1, &f, &cs, &sn);
        s[k - 1] = t1;

        if (k != l)
        {
          f = -sn * e[k - 2];
          e[k - 2] = cs * e[k - 2];
        }
      }
    }
    //
    //  Split at negligible S(L).
    //
    else if (kase == 2)
    {
      f = e[l - 2];
      e[l - 2] = 0.0;

      for (k = l; k <= mn; k++)
      {
        t1 = s[k - 1];
        drotg(&t1, &f, &cs, &sn);
        s[k - 1] = t1;
        f = -sn * e[k - 1];
        e[k - 1] = cs * e[k - 1];
      }
    }
    //
    //  Perform one QR step.
    //
    else if (kase == 3)
    {
      //
      //  Calculate the shift.
      //
      scale = max(fabs(s[mn - 1]),
                  max(fabs(s[mn - 2]),
                      max(fabs(e[mn - 2]),
                          max(fabs(s[l - 1]), fabs(e[l - 1])))));

      sm = s[mn - 1] / scale;
      smm1 = s[mn - 2] / scale;
      emm1 = e[mn - 2] / scale;
      sl = s[l - 1] / scale;
      el = e[l - 1] / scale;
      b = ((smm1 + sm) * (smm1 - sm) + emm1 * emm1) / 2.0;
      c = (sm * emm1) * (sm * emm1);
      shift = 0.0;

      if (b != 0.0 || c != 0.0)
      {
        shift = sqrt(b * b + c);
        if (b < 0.0)
        {
          shift = -shift;
        }
        shift = c / (b + shift);
      }

      f = (sl + sm) * (sl - sm) - shift;
      g = sl * el;
      //
      //  Chase zeros.
      //
      mm1 = mn - 1;

      for (k = l; k <= mm1; k++)
      {
        drotg(&f, &g, &cs, &sn);

        if (k != l)
        {
          e[k - 2] = f;
        }

        f = cs * s[k - 1] + sn * e[k - 1];
        e[k - 1] = cs * e[k - 1] - sn * s[k - 1];
        g = sn * s[k];
        s[k] = cs * s[k];

        drotg(&f, &g, &cs, &sn);
        s[k - 1] = f;
        f = cs * e[k - 1] + sn * s[k];
        s[k] = -sn * e[k - 1] + cs * s[k];
        g = sn * e[k];
        e[k] = cs * e[k];
      }
      e[mn - 2] = f;
      iter = iter + 1;
    }
    //
    //  Convergence.
    //
    else if (kase == 4)
    {
      //
      //  Make the singular value nonnegative.
      //
      if (s[l - 1] < 0.0)
      {
        s[l - 1] = -s[l - 1];
      }
      //
      //  Order the singular value.
      //
      for (;;)
      {
        if (l == mm)
        {
          break;
        }

        if (s[l] <= s[l - 1])
        {
          break;
        }

        t = s[l - 1];
        s[l - 1] = s[l];
        s[l] = t;

        l = l + 1;
      }
      iter = 0;
      mn = mn - 1;
    }
  }

  return info;
}
