#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <bs.h>

#if defined(_MSC_VER)
  #define INLINE _inline  // __inline in newer versions
  #define RESTRICT __restrict
#else
  #define INLINE inline
  #define RESTRICT restrict
#endif


//-----------------------------------------------------------------------------
// debug stuff (remove later)
#include <stdio.h>

void print_a_and_b(double first[5], double last[5],
                   double *A, double  *b, int M)
{
    int i;

    printf("\nfirst: [ %f  %f  %f  %f  %f ]\n",
           first[0], first[1], first[2], first[3], first[4]);

    for (i=0; i<M; i++)
        printf("row %d : | %f  %f  %f |    | %f |\n",
               i, A[3*i+0], A[3*i+1], A[3*i+2], b[i]);

    printf("last: [ %f  %f  %f  %f  %f ]\n",
           last[0], last[1], last[2], last[3], last[4]);

}


//-----------------------------------------------------------------------------
// search functions & helpers
//
// These all return i such that x>= values[i] and x<values[i+1].
// Return -1 if x < values[0].
// Return n-1 if x >= values[n-1].
//-----------------------------------------------------------------------------

// Linear search starting from guess that
// values[start] <= x < values[start+1].
static int find_index_from(double *values, int n, double x, int start)
{
    int i;

    if (start <= -1) {
        // search down
        i = 0;
        while (i < n && x >= values[i]) i++;
        return i-1;
    }
    else if (start >= n-1) {
        // search down
        i = n - 1;
        while (i > -1 && x < values[i]) i--;
        return i;
    }
    else if (x >= values[start]) {
        // search up
        i = start + 1;
        while (i < n && x >= values[i]) i++;
        return i-1;
    }
    else {
        // search down
        i = start - 1;
        while (i > -1 && x < values[i]) i--;
        return i;
    }
}


// find index using binary search
static int find_index_binary(double *values, int n, double x)
{
    int lo, hi, mid;

    lo = 0;
    hi = n;
    mid = n/2;

    if (x < values[0]) return -1;
    if (x >= values[n-1]) return n-1;

    while (hi - lo > 1) {
        if (x >= values[mid]) lo = mid;
        else                  hi = mid;
        mid = lo + (hi - lo) / 2;
    }

    return mid;
}


static int is_monotonic(bs_array x)
{
    int i;
    int ok;
    ok = 1;
    for (i=1; i<x.size; i++) {
        ok &= (x.data[i*x.stride] >= x.data[(i-1)*x.stride]);
    }
    return ok;
}

static int min_points(bs_bctype left, bs_bctype right)
{
    // one additional point needed for each not-a-knot condition.
    return 2 + (left == BS_NOTAKNOT) + (right == BS_NOTAKNOT);
}

//-----------------------------------------------------------------------------
// knots & constants
//-----------------------------------------------------------------------------

// fill spline knots based on x array (includes padding on either
// end of array).
static double* alloc_knots(bs_array x)
{
    int N;
    double *knots;
    int i;
    
    N = x.size;
    knots = malloc((N + 5) * sizeof(double));

    // move pointer past initial two-element padding.
    knots += 2;

    // copy x into main part of knots
    for (i=0; i < N; i++) knots[i] = x.data[i * x.stride];

    // fill padded area before beginning
    knots[-2] = knots[0] - 2.0 * (knots[1] - knots[0]);
    knots[-1] = knots[0] - 1.0 * (knots[1] - knots[0]);

    // fill padded area after end.
    knots[N]   = knots[N-1] + 1.0 * (knots[N-1] - knots[N-2]);
    knots[N+1] = knots[N-1] + 2.0 * (knots[N-1] - knots[N-2]);
    knots[N+2] = knots[N-1] + 3.0 * (knots[N-1] - knots[N-2]);

    return knots;
}


static void free_knots(double *knots) {
    free(knots - 2);
}


// constants used when evaluating a spline.
// constants + 4*i is a pointer to the four constants used
// when evaluating the spline in the range knots[i] <= x < knots[i+1].
static double* alloc_constants(double *knots, int n)
{
    int i;
    double *constants;

    constants = malloc(4 * n * sizeof(double));

    for (i=0; i<n; i++) {
        constants[4*i+0] = 1.0 / ((knots[i+1] - knots[i-2]) *
                                  (knots[i+1] - knots[i-1]) *
                                  (knots[i+1] - knots[i  ]));
        
        constants[4*i+1] = 1.0 / ((knots[i+3] - knots[i  ]) *
                                  (knots[i+2] - knots[i  ]) *
                                  (knots[i+1] - knots[i  ]));

        constants[4*i+2] = 1.0 / ((knots[i+2] - knots[i-1]) *
                                  (knots[i+1] - knots[i-1]) *
                                  (knots[i+1] - knots[i  ]));

        constants[4*i+3] = 1.0 / ((knots[i+2] - knots[i-1]) *
                                  (knots[i+2] - knots[i  ]) *
                                  (knots[i+1] - knots[i  ]));
    }

    return constants;
}


//-----------------------------------------------------------------------------
// Compute the 4 basis functions that are nonzero, assuming t[i] <= x < t[i+1].
// These are: b_{3, i-3}(x), b_{3, i-2}(x), b_{3, i-1}(x), b_{3, i}(x)
//
// This is faster than computing them separately, as some parts of the
// calculation are shared. These can be derived by "manually inlining"
// recursive function calls in the formula for the basis function.
// (See tests for recursive version).
//
// consts[4*i] through consts[4*i+3] stores four constants used in calculation
// (constants are a function of knot spacings).
//
// t indicies from (i-2) to (i+3) are used.
//-----------------------------------------------------------------------------

static void b3nonzeros(double x, int i, double* RESTRICT t,
                       double* RESTRICT consts, double* RESTRICT out)
{
    double* RESTRICT c = consts + 4*i;

    double dx1 = x - t[i-2];
    double dx2 = x - t[i-1];
    double dx3 = x - t[i];
    double dx4 = t[i+1] - x;
    double dx5 = t[i+2] - x;
    double dx6 = t[i+3] - x;

    double tmp1 = dx4 * dx4 * c[0];
    double tmp2 = dx3 * dx3 * c[1];
    double tmp3 = dx2 * dx4 * c[2] + dx5 * dx3 * c[3];
    
    out[0] = dx4 * tmp1;
    out[1] = dx1 * tmp1 + dx5 * tmp3;
    out[2] = dx6 * tmp2 + dx2 * tmp3;
    out[3] = dx3 * tmp2;
}


// derivatives of previous function
static void db3nonzeros(double x, int i, double* RESTRICT t,
                        double* RESTRICT consts, double out[4])
{
    double* RESTRICT c = consts + 4*i;

    double dx1 = x - t[i-2];
    double dx2 = x - t[i-1];
    double dx3 = x - t[i];
    double dx4 = t[i+1] - x;
    double dx5 = t[i+2] - x;
    double dx6 = t[i+3] - x;
  
    double tmp1 = dx4 * c[0];
    double tmp2 = dx3 * c[1];
    double tmp3 = dx2 * c[2];
    double tmp4 = dx5 * c[3];
  
    out[0] = -3.0 * dx4 * tmp1;
  
    out[1] = ((        dx4 - 2.0 * dx1) * tmp1 +
              (-       dx4 -       dx5) * tmp3 +
              (- 2.0 * dx3 +       dx5) * tmp4 +
              dx5 * dx4 * c[2]);

    out[2] = ((-     dx3 + 2.0 * dx6) * tmp2 +
              (2.0 * dx4 -       dx2) * tmp3 +
              (      dx3 +       dx2) * tmp4
              - dx2 * dx3 * c[3]);

    out[3] = 3.0 * dx3 * tmp2;
}


// second derivatives
static void d2b3nonzeros(double x, int i, double* RESTRICT t,
                         double* RESTRICT consts, double out[4])
{
    double* RESTRICT c = consts + 4*i;

    double dx1 = x - t[i-2];
    double dx2 = x - t[i-1];
    double dx3 = x - t[i];
    double dx4 = t[i+1] - x;
    double dx5 = t[i+2] - x;
    double dx6 = t[i+3] - x;

    out[0] = 6.0 * dx4 * c[0];

    out[1] = (- 2.0 * dx4 * c[0]
              - 2.0 * (dx4 - dx1) * c[0]
              -       (dx4 - dx2) * c[2]
              +       (-dx5 - dx4) * c[2]
              -       (dx5 - dx2) * c[2]
              - 2.0 * (dx5 - dx3) * c[3]
              - 2.0 * dx5 * c[3]);

    out[2] = (- 2.0 * dx3 * c[1]
              + 2.0 * (dx6 - dx3) * c[1] 
              + 2.0 * (dx4 - dx2) * c[2]
              - 2.0 * dx2 * c[2]
              +       (dx5 - dx3) * c[3]
              -       (dx2 + dx3) * c[3]
              +       (dx5 - dx2) * c[3]);

    out[3] = 6.0 * dx3 * c[1];
}

// third derivatives
static void d3b3nonzeros(int i, double* RESTRICT consts, double out[4])
{
    double* RESTRICT c = consts + 4*i;

    out[0] = -6.0 * c[0];
    out[1] =  6.0 * (c[0] + c[2] + c[3]);
    out[2] = -6.0 * (c[1] + c[2] + c[3]);
    out[3] =  6.0 * c[1];
}


//-----------------------------------------------------------------------------
// unit basis versions of b3nonzeros and friends
// knot locations in this basis are [0, 1, ..., N-1]
// For 0 <= x < 1 return b_{-3}(x), b_{-2}(x), b_{-1}(x), b_0(x)
// For i <= x < i+1 subtract i from x first. to get b_{i-3}(x), b_{i-2}(x), ...
// (works because all the basis functions are the same with a shift,
//  so b_{-3}(x) = b_{i-3}(i+x).
//-----------------------------------------------------------------------------

static const double ONESIXTH  = 0.1666666666666666666;

static void b3unonzeros(double x, double out[4])
{
    double dx1 = x + 2.0;
    double dx2 = x + 1.0;
    double dx4 = 1.0 - x;
    double dx5 = 2.0 - x;
    double dx6 = 3.0 - x;

    double tmp1 = ONESIXTH * dx4 * dx4;
    double tmp2 = ONESIXTH * x * x;
    double tmp3 = ONESIXTH * (dx2 * dx4 + dx5 * x);

    out[0] = dx4 * tmp1;
    out[1] = dx1 * tmp1 + dx5 * tmp3;
    out[2] = dx6 * tmp2 + dx2 * tmp3;
    out[3] = x   * tmp2;
}

/*

// derivatives of previous function
static void db3unonzeros(double x, double out[4])
{
    double dx1 = x + 2.0;
    double dx2 = x + 1.0;
    double dx4 = 1.0 - x;
    double dx5 = 2.0 - x;
    double dx6 = 3.0 - x;

    double tmp1 = ONESIXTH * dx4;
    double tmp2 = ONESIXTH * x;
    double tmp3 = ONESIXTH * dx2;
    double tmp4 = ONESIXTH * dx5;

    out[0] = -3.0 * dx4 * tmp1;

    out[1] = ((        dx4 - 2.0 * dx1) * tmp1 +
              (-       dx4 -       dx5) * tmp3 +
              (- 2.0 * x +       dx5) * tmp4 +
              ONESIXTH * dx5 * dx4);

    out[2] = ((-     x + 2.0 * dx6) * tmp2 +
              (2.0 * dx4 -       dx2) * tmp3 +
              (      x +       dx2) * tmp4
              - ONESIXTH * dx2 * x);

    out[3] = 3.0 * x * tmp2;
}


// second derivatives
static void d2b3unonzeros(double x, double out[4])
{
    double dx1 = x + 2.0;
    double dx2 = x + 1.0;
    double dx4 = 1.0 - x;
    double dx5 = 2.0 - x;
    double dx6 = 3.0 - x;

    out[0] = ONESIXTH * 6.0 * dx4;

    out[1] = ONESIXTH * (- 2.0 * dx4
                         - 2.0 * (dx4 - dx1)
                         -       (dx4 - dx2)
                         +       (-dx5 - dx4)
                         -       (dx5 - dx2)
                         - 2.0 * (dx5 - x)
                         - 2.0 * dx5);

    out[2] = ONESIXTH * (- 2.0 * x
                         + 2.0 * (dx6 - x)
                         + 2.0 * (dx4 - dx2)
                         - 2.0 * dx2
                         +       (dx5 - x)
                         -       (dx2 + x)
                         +       (dx5 - dx2));

    out[3] = ONESIXTH * 6.0 * x;
}

// third derivatives
static void d3b3unonzeros(double out[4])
{
    out[0] = -1.0;
    out[1] =  3.0;
    out[2] = -3.0;
    out[3] =  1.0;
}

*/

//-----------------------------------------------------------------------------
// solve_simple()
//
// Solve A * x = b for x. The solution is stored in b.
//
// A is an almost tridiagonal n x n matrix with this form:
//
// | x x x              |
// | x x x              |
// |   x x x            |
// |        ...         |
// |            x x x   |
// |              x x x |
// |              x x x |
//
// Rows are contiguous in memory: e.g., A[0] through A[2]
// stores the second row (first row with three elements).
//
//-----------------------------------------------------------------------------

/*
static void solve_simple(double* RESTRICT A, double* RESTRICT b, int n)
{
    int i;

    // divide first row by upper left element
    double t = A[0];
    b[0] /= t;
    A[2] /= t;  
    A[1] /= t;
    A[0] = 1.0; // but not used again.

    // subtract (first element of row 1) x (row 0) from row 1
    // to eliminate first element of row 1.
    t = A[3*1+0];
    b[1]     -= t * b[0];
    A[3*1+2] -= t * A[2];
    A[3*1+1] -= t * A[1];
    A[3*1+0] = 0.0; // but not used again.

    // divide row 1 by first nonzero element, to set it to 1.
    t = A[3*1+1];
    b[1]     /= t;
    A[3*1+2] /= t;
    A[3*1+1] = 1.0; // but not used again.

    for (i=2; i<n-1; i++) {

        // subtract (first element of new row) * (previous row) from new row
        // to eliminate first element.
        t = A[3*i+0];
        b[i]        -= t * b[i-1];
        // A[3*i+2] -= t * 0.0  // no-op b/c previous row is zero.
        A[3*i+1]    -= t * A[3*(i-1)+2];
        A[3*i+0] = 0.0;  // (previous row is 1.0) but not used again.

        // divide new row by first non-zero element
        t = A[3*i+1];
        b[i]     /= t;
        A[3*i+2] /= t;
        A[3*i+1] = 1.0;
    }

    // last row is different:
    // subtract first element of last row * 3rd to last row from last row
    b[n-1]          -= A[3*(n-1)+0] * b[n-3];
    // A[3*(n-1)+2] -= A[3*(n-1)+0] * 0.0; // no-op
    A[3*(n-1)+1]    -= A[3*(n-1)+0] * A[3*(n-3)+2];
    A[3*(n-1)+0]    = 0.0;

    // subtract first non-zero element * previous row from last row
    b[n-1]       -= A[3*(n-1)+1] * b[n-2];
    A[3*(n-1)+2] -= A[3*(n-1)+1] * A[3*(n-2)+2];
    A[3*(n-1)+1] = 0.0;

    // divide row by 1st non-zero element
    b[n-1]       /= A[3*(n-1)+2];
    A[3*(n-1)+2] =  1.0;

    // back substitute
    for (i=n-2; i>0; i--) {
        b[i] -= b[i+1] * A[3*i+2];
    }

    // first row is different
    b[0] -= b[1] * A[1] + b[2] * A[2];
}
*/


//-----------------------------------------------------------------------------
// solve()
//
// Solve A * x = b for x. The solution is stored in b.
//
// A is a matrix like this:
//
// | x x x x x          |
// | x x x              |
// |   x x x            |
// |        ...         |
// |          x x x     |
// |            x x x   |
// |              x x x |
// |          x x x x x |
//
// A is stored compactly in 3*n elements, with row i corresponding to A[3*i],
// with the exception of the first and last rows which are passed in
// separately because they are too large to be stored this way.
//
// Note that the first 3 and last 3 elements of A are initially empty as
// these row values are stored in `first` and `last`. In fact the last 3
// elements of A are not used at all.
//
//-----------------------------------------------------------------------------

typedef struct {
    double *first;
    double *rows; // first + 5
    double *last; // first + 5 + 3*(M-1)
} banded_matrix;


// allocate storate for a banded matrix.
// M is total number of rows, including first and last.
static banded_matrix alloc_banded_matrix(int M)
{
    double* first = malloc((5 + 3*(M-1) + 5) * sizeof(double));
    banded_matrix m = {first, first + 5, first + 5 + 3*(M-1)};
    return m;
}


static void free_banded_matrix(banded_matrix A)
{
    free(A.first);
    A.first = NULL;
    A.rows = NULL;
    A.last = NULL;
}


static void copy_banded_matrix(banded_matrix dst, banded_matrix src, int M)
{
    size_t nbytes = (10 + 3*(M-1)) * sizeof(double);
    memcpy(dst.first, src.first, nbytes);
}


static void solve(banded_matrix mat, double* RESTRICT b, int n)
{
    int i;
    double tmp;
    double* RESTRICT first = mat.first;
    double* RESTRICT A = mat.rows;
    double* RESTRICT last = mat.last;

    // rows 1, 2, 3: divide by first non-zero
    //
    // x x x x x | y       x x x x x | y
    // x x x     | y       1 x x     | y
    //   x x x   | y  -->    1 x x   | y
    //     x x x | y           1 x x | y

    for (i=1; i<4; i++) {
        b[i]     /= A[3*i];
        A[3*i+2] /= A[3*i];
        A[3*i+1] /= A[3*i];
        A[3*i]   = 1.0;
    }

    // eliminate first two elements of first row and divide by first non-zero.
    //
    // x x x x x | y       0 0 1 x x | y
    // 1 x x     | y       1 x x     | y
    //   1 x x   | y  -->    1 x x   | y
    //     1 x x | y           1 x x | y
    b[0]     -= first[0] * b[1];
    first[2] -= first[0] * A[3*1+2];
    first[1] -= first[0] * A[3*1+1];
    first[0] = 0.0;

    b[0]     -= first[1] * b[2];
    first[3] -= first[1] * A[3*2+2];
    first[2] -= first[1] * A[3*2+1];
    first[1] = 0.0;

    b[0]     /= first[2];
    first[4] /= first[2];
    first[3] /= first[2];
    first[2] = 1.0;

    // reduce row 3
    //
    // 0 0 1 x x | y       0 0 1 x x | y
    // 1 x x     | y       1 x x     | y
    //   1 x x   | y  -->    1 x x   | y
    //     1 x x | y           0 1 x | y   
    b[3]     -= A[3*3+0] * b[0];
    A[3*3+2] -= A[3*3+0] * first[4];
    A[3*3+1] -= A[3*3+0] * first[3];
    A[3*3+0] = 0.0;

    b[3]     /= A[3*3+1];
    A[3*3+2] /= A[3*3+1];
    A[3*3+1] = 1.0;

    // permute first three rows:
    // 0 0 1 x x | y       1 x x     | y
    // 1 x x     | y         1 x x   | y
    //   1 x x   | y  -->      1 x x | y
    //     0 1 x | y           0 1 x | y
    tmp = b[0];
    b[0] = b[1];
    A[3*0+0] = A[3*1+0];
    A[3*0+1] = A[3*1+1];
    A[3*0+2] = A[3*1+2];

    b[1] = b[2];
    A[3*1+0] = A[3*2+0];
    A[3*1+1] = A[3*2+1];
    A[3*1+2] = A[3*2+2];

    b[2] = tmp;
    A[3*2+0] = first[2];
    A[3*2+1] = first[3];
    A[3*2+2] = first[4];

    // reduce rest of the middle rows
    for (i=4; i<n-1; i++) {
        b[i]     -= A[3*i+0] * b[i-1];
        A[3*i+1] -= A[3*i+0] * A[3*(i-1)+2];
        A[3*i+0] = 0.0;

        b[i]     /= A[3*i+1];
        A[3*i+2] /= A[3*i+1];
        A[3*i+1] = 1.0;
    }

    // we now have, e.g.,
    // 1 x x         | y
    //   1 x x       | y
    //     1 x x     | y  (n-5)
    //     0 1 x     | y  (n-4)
    //       0 1 x   | y  (n-3)
    //         0 1 x | y  (n-2)
    //     x x x x x | y  (n-1)

    // eliminate first element of last row using the (n-5)th row.
    b[n-1] -= last[0] * b[n-5];
    if (n-5 < 3) {
        last[2] -= last[0] * A[3*(n-5)+2];
        last[1] -= last[0] * A[3*(n-5)+1];
    }
    else {
        last[1] -= last[0] * A[3*(n-5)+2];
    }
    last[0] = 0.0;

    // eliminate second element of last row using the (n-4)th row.
    b[n-1] -= last[1] * b[n-4];
    if (n-4 < 3) {
        last[3] -= last[1] * A[3*(n-4)+2];
        last[2] -= last[1] * A[3*(n-4)+1];
    }
    else {
        last[2] -= last[1] * A[3*(n-4)+2];
    }
    last[1] = 0.0;

    // eliminate third element of last row using the (n-3)rd row.
    b[n-1] -= last[2] * b[n-3];
    if (n-3 < 3) {
        last[4] -= last[2] * A[3*(n-3)+2];
        last[3] -= last[2] * A[3*(n-3)+1];
    }
    else {
        last[3] -= last[2] * A[3*(n-3)+2];
    }
    last[2] = 0.0;

    // eliminate forth element
    b[n-1] -= last[3] * b[n-2];
    last[4] -= last[3] * A[3*(n-2)+2];
    last[3] = 0.0;

    // normalize last row
    b[n-1] /= last[4];
    last[4] = 1.0;

    // back-substitute
    for (i=n-2; i>=3; i--) {
        b[i] -= b[i+1] * A[3*i+2];
    }

    // we now have:
    // 1 x x           | y
    //   1 x x         | y
    //     1 x x       | y
    //       1         | y
    //         1       | y
    //          ...
    //
    // eliminate the remaining elements.
    b[2] -= b[3] * A[3*2+1] + b[4] * A[3*2+2];
    b[1] -= b[2] * A[3*1+1] + b[3] * A[3*1+2];
    b[0] -= b[1] * A[3*0+1] + b[2] * A[3*0+2];
}


//-----------------------------------------------------------------------------
// finding coefficients
//-----------------------------------------------------------------------------

static void notaknot_row(double *consts, int i, double row[5])
{
    int j;
    double buf[4];

    d3b3nonzeros(i-1, consts, row);
    d3b3nonzeros(i, consts, buf);
    row[4] = 0.0;
    for (j=0; j<4; j++) {
        row[j+1] -= buf[j];
    }
}


static void fill_banded_matrix(banded_matrix A, double* RESTRICT knots,
                               double* RESTRICT consts, int N,
                               bs_bctype bctypes[2])
{
    int i;
    double* RESTRICT first = A.first;
    double* RESTRICT rows = A.rows;
    double* RESTRICT last = A.last;

    // Left boundary condition
    switch (bctypes[0]) {
    case BS_DERIV1:
        db3nonzeros(knots[0], 0, knots, consts, first);
        first[3] = first[4] = 0.0;
        break;
    case BS_DERIV2:
        d2b3nonzeros(knots[0], 0, knots, consts, first);
        first[3] = first[4] = 0.0;
        break;
    case BS_NOTAKNOT:
        notaknot_row(consts, 1, first);
    }
    
    // fill rows 1 through M-1 with values of b_{i-3}, b_{i-2}, b{i-1}
    // at knot i.
    for (i=0; i<N; i++) {
        b3nonzeros(knots[i], i, knots, consts, rows + 3*(i+1));
    }

    // Right boundary condition
    switch (bctypes[1]) {
    case BS_DERIV1:
        db3nonzeros(knots[N-1], N-1, knots, consts, last);
        for (i=4; i>1; i--) last[i] = last[i-2];
        last[0] = last[1] = 0.0;
        break;
    case BS_DERIV2:
        d2b3nonzeros(knots[N-1], N-1, knots, consts, last);
        for (i=4; i>1; i--) last[i] = last[i-2];
        last[0] = last[1] = 0.0;
        break;
    case BS_NOTAKNOT:
        notaknot_row(consts, N-2, last);
    }
}

    
// Find spline coefficients along one dimension.
// knots and consts are as belong to a spline.
// A should be allocated with M = values.size + 2.
// coeffs should be size M.
static void find_1d_coefficients(banded_matrix A,
                                 bs_array values, double bcvalues[2],
                                 double* RESTRICT coeffs)
{
    int i;
    int N = values.size;
    int M = N+2;
    
    coeffs[0] = bcvalues[0];
    for (i=0; i<N; i++) {
        coeffs[i+1] = values.data[i * values.stride];
    }
    coeffs[M-1] = bcvalues[1];

    solve(A, coeffs, M);
}

//-----------------------------------------------------------------------------
// spline1d
//-----------------------------------------------------------------------------

bs_errorcode bs_spline1d_create(bs_array x, bs_array y, bs_bcs bcs,
                                bs_exts exts, bs_spline1d **out)
{
    int N, M;
    bs_spline1d* spline;
    banded_matrix A;
    bs_bctype bctypes[2] = {bcs.left.type, bcs.right.type};
    double bcvalues[2];

    *out = NULL;  // In case of error, ensure that output pointer is NULL.

    // checks
    if (x.size != y.size) return BS_SIZEMISMATCH;
    if (!is_monotonic(x)) return BS_NOTMONOTONIC;
    if (x.size < min_points(bcs.left.type, bcs.right.type))
        return BS_TOOFEWPOINTS;

    spline = malloc(sizeof(bs_spline1d));
  
    N = x.size;
    M = N + 2;

    spline->knots = alloc_knots(x);
    spline->n = N;
    spline->exts = exts;
    spline->consts = alloc_constants(spline->knots, N);
    spline->coeffs = malloc(M * sizeof(double));
    
    A = alloc_banded_matrix(M);

    fill_banded_matrix(A, spline->knots, spline->consts, N, bctypes);

    bcvalues[0] = (bcs.left.type == BS_NOTAKNOT) ? 0.0 : bcs.left.value;
    bcvalues[1] = (bcs.right.type == BS_NOTAKNOT) ? 0.0 : bcs.right.value;

    find_1d_coefficients(A, y, bcvalues, spline->coeffs);

    free_banded_matrix(A);

    *out = spline;
    return BS_OK;
}


void bs_spline1d_free(bs_spline1d* spline)
{
    if (spline != NULL) {
        free_knots(spline->knots);
        free(spline->consts);
        free(spline->coeffs);
        free(spline);
    }
}


bs_errorcode bs_spline1d_eval(bs_spline1d *spline, bs_array x, bs_array out)
{
    int i;
    int j;
    double xval;
    double b3vals[4];

    // for first index, it could be anywhere, so use binary search
    i = find_index_binary(spline->knots, spline->n, x.data[0]);

    for (j=0; j<x.size; j++) {
        xval = x.data[j*x.stride];
        i = find_index_from(spline->knots, spline->n, xval, i);

        // index outside left boundary
        if (i == -1) {
            switch (spline->exts.left.type) {
            case BS_EXTRAPOLATE:
                i = 0;
                break;
            case BS_CONSTANT:
                i = 0;
                xval = spline->knots[0];
                break;
            case BS_VALUE:
                out.data[j * out.stride] = spline->exts.left.value;
                continue;
            case BS_RAISE:
                return BS_DOMAINERROR;
            }
        }

        // index outside right boundary
        else if (i == spline->n - 1) {
            switch (spline->exts.right.type) {
            case BS_EXTRAPOLATE:
                i = spline->n - 2;
                break;
            case BS_CONSTANT:
                i = spline->n-2;
                xval = spline->knots[spline->n-1];
                break;
            case BS_VALUE:
                out.data[j * out.stride] = spline->exts.right.value;
                continue;
            case BS_RAISE:
                return BS_DOMAINERROR;
            }
        }

        // if we get this far, we're either extrapolating or xval is in range.
        b3nonzeros(xval, i, spline->knots, spline->consts, b3vals);
        out.data[j*out.stride] = (spline->coeffs[i]   * b3vals[0] +
                                  spline->coeffs[i+1] * b3vals[1] +
                                  spline->coeffs[i+2] * b3vals[2] +
                                  spline->coeffs[i+3] * b3vals[3]);
    }

    return BS_OK;
}


//-----------------------------------------------------------------------------
// spline2d
//-----------------------------------------------------------------------------

// get the i-th bc value
static double get_bcarray_value(bs_bcarray bc, int i)
{
    return (bc.type == BS_NOTAKNOT) ? 0.0 : bc.data[i * bc.stride];
}

static int bcarray_size_match(bs_bcarray bc, int size)
{
    // if not-a-knot, value is not used, so size is OK.
    return (bc.type == BS_NOTAKNOT) ? 1 : (bc.size == size);
}


bs_errorcode bs_spline2d_create(bs_array x, bs_array y, bs_array2d z,
                                bs_bcarray_pair xbcs, bs_bcarray_pair ybcs,
                                bs_exts xexts, bs_exts yexts,
                                bs_spline2d **out)
{
    int i, j;
    int nx, mx, ny, my;
    bs_bctype bctypes[2] = {ybcs.left.type, ybcs.right.type};
    double bcvalues[2];
    bs_spline2d* spline;
    double *coeffs;
    banded_matrix A, Awork;
    bs_array zslice, coeffs_slice;
    double *buf;
    
    *out = NULL;  // In case of error, ensure that output pointer is NULL.

    if ((x.size != z.sizes[0]) || (y.size != z.sizes[1]))
        return BS_SIZEMISMATCH;
    if (!is_monotonic(x) || !is_monotonic(y))
        return BS_NOTMONOTONIC;
    if ((x.size < min_points(xbcs.left.type, xbcs.right.type)) ||
        (y.size < min_points(ybcs.left.type, ybcs.right.type)))
        return BS_TOOFEWPOINTS;
    
    // check if boundary condition sizes match x and y sizes.
    if (!(bcarray_size_match(xbcs.left, y.size) &&
          bcarray_size_match(xbcs.right, y.size) &&
          bcarray_size_match(ybcs.left, x.size) &&
          bcarray_size_match(ybcs.left, x.size)))
        return BS_BCSIZEMISMATCH;

    spline = malloc(sizeof(bs_spline2d));
  
    nx = x.size;
    mx = nx + 2;

    ny = y.size;
    my = ny + 2;

    spline->xknots = alloc_knots(x);
    spline->xconsts = alloc_constants(spline->xknots, nx);
    spline->nx = nx;
    spline->xexts = xexts;

    spline->yknots = alloc_knots(y);
    spline->yconsts = alloc_constants(spline->yknots, ny);
    spline->ny = ny;
    spline->yexts = yexts;

    coeffs = malloc(mx * my * sizeof(double));
    
    // find coefficients along y (fast axis)
    A     = alloc_banded_matrix(my);
    Awork = alloc_banded_matrix(my);

    fill_banded_matrix(A, spline->yknots, spline->yconsts, spline->ny,
                       bctypes);

    for (i=0; i<nx; i++) {
        zslice.data = z.data + z.strides[0]*i;
        zslice.size = z.sizes[1];
        zslice.stride = z.strides[1];

        bcvalues[0] = get_bcarray_value(ybcs.left, i);
        bcvalues[1] = get_bcarray_value(ybcs.right, i);

        copy_banded_matrix(Awork, A, my);
        find_1d_coefficients(Awork, zslice, bcvalues, coeffs+(i*my));
    }

    free_banded_matrix(A);
    free_banded_matrix(Awork);

    // find coefficients along x (slow axis)
    A     = alloc_banded_matrix(mx);
    Awork = alloc_banded_matrix(mx);

    bctypes[0] = xbcs.left.type;
    bctypes[1] = xbcs.right.type;
    fill_banded_matrix(A, spline->xknots, spline->xconsts, spline->nx,
                       bctypes);

    buf = malloc(mx * sizeof(double));
    for (i=0; i<my; i++) {
        // for this slice in constant y, the target values are the
        // `nx` coefficients we just found. They are strided in
        // `coeffs` by `my`.
        copy_banded_matrix(Awork, A, mx);
        coeffs_slice.data = coeffs + i;
        coeffs_slice.size = nx;
        coeffs_slice.stride = my;

        bcvalues[0] = get_bcarray_value(xbcs.left, i);
        bcvalues[1] = get_bcarray_value(xbcs.right, i);

        find_1d_coefficients(Awork, coeffs_slice, bcvalues, buf);
        
        // the results in `buf` are contiguous in x, but we need to
        // copy them back into the coefficients array strided.
        for (j=0; j<mx; j++) coeffs[i+my*j] = buf[j];
    }

    free_banded_matrix(A);
    free_banded_matrix(Awork);
    free(buf);

    spline->coeffs = coeffs;
    *out = spline;

    return BS_OK;
}


void bs_spline2d_free(bs_spline2d* spline)
{
    if (spline != NULL) {
        free_knots(spline->xknots);
        free(spline->xconsts);
        free_knots(spline->yknots);
        free(spline->yconsts);
        free(spline->coeffs);
        free(spline);
    }
}


bs_errorcode bs_spline2d_eval(bs_spline2d *spline, bs_array x, bs_array y,
                              bs_array2d out)
{
    // for first index, it could be anywhere, so use binary search
    int i = find_index_binary(spline->xknots, spline->nx, x.data[0]);
    int j0 = find_index_binary(spline->yknots, spline->ny, y.data[0]);
    int j, k, l;
    int my = spline->ny + 2; // for indexing coeffs.
    double xb3vals[4];
    double yb3vals[4];
    double xval, yval;

    for (k=0; k<x.size; k++) {
        xval = x.data[k*x.stride];
        i = find_index_from(spline->xknots, spline->nx, xval, i);

        // index outside left boundary
        if (i == -1) {
            switch (spline->xexts.left.type) {
            case BS_EXTRAPOLATE:
                i = 0;
                break;
            case BS_CONSTANT:
                i = 0;
                xval = spline->xknots[0];
                break;
            case BS_VALUE:
                for (l=0; l<y.size; l++) {
                    out.data[k * out.strides[0] + l * out.strides[1]] =
                        spline->xexts.left.value;
                }
                continue;
            case BS_RAISE:
                return BS_DOMAINERROR;
            }
        }

        // index outside right boundary
        else if (i == spline->nx - 1) {
            switch (spline->xexts.right.type) {
            case BS_EXTRAPOLATE:
                i = spline->nx - 2;
                break;
            case BS_CONSTANT:
                i = spline->nx - 2;
                xval = spline->xknots[spline->nx-1];
                break;
            case BS_VALUE:
                for (l=0; l<y.size; l++) {
                    out.data[k * out.strides[0] + l * out.strides[1]] =
                        spline->xexts.right.value;
                }
                continue;
            case BS_RAISE:
                return BS_DOMAINERROR;
            }
        }

        // get basis function values for x coordinate
        b3nonzeros(xval, i, spline->xknots, spline->xconsts, xb3vals);

        // x value is in range (or extrapolating); loop over y values:
        j = j0;
        for (l=0; l<y.size; l++) {
            yval = y.data[l*y.stride];
            j = find_index_from(spline->yknots, spline->ny, yval, j);

            // index outside left boundary
            if (j == -1) {
                switch (spline->yexts.left.type) {
                case BS_EXTRAPOLATE:
                    j = 0;
                    break;
                case BS_CONSTANT:
                    j = 0;
                    yval = spline->yknots[0];
                    break;
                case BS_VALUE:
                    out.data[k * out.strides[0] + l * out.strides[1]] =
                        spline->yexts.left.value;
                    continue;
                case BS_RAISE:
                    return BS_DOMAINERROR;
                }
            }

            // index outside right boundary
            else if (j == spline->ny - 1) {
                switch (spline->yexts.right.type) {
                case BS_EXTRAPOLATE:
                    j = spline->ny - 2;
                    break;
                case BS_CONSTANT:
                    j = spline->ny - 2;
                    yval = spline->yknots[spline->ny-1];
                    break;
                case BS_VALUE:
                    out.data[k * out.strides[0] + l * out.strides[1]] =
                        spline->yexts.right.value;
                    continue;
                case BS_RAISE:
                    return BS_DOMAINERROR;
                }
            }

            // get basis function values for y coordinate
            b3nonzeros(yval, j, spline->yknots, spline->yconsts, yb3vals);

            out.data[k * out.strides[0] + l * out.strides[1]] =
                (spline->coeffs[(i  )*my+j]   * xb3vals[0] * yb3vals[0] +
                 spline->coeffs[(i  )*my+j+1] * xb3vals[0] * yb3vals[1] +
                 spline->coeffs[(i  )*my+j+2] * xb3vals[0] * yb3vals[2] +
                 spline->coeffs[(i  )*my+j+3] * xb3vals[0] * yb3vals[3] +

                 spline->coeffs[(i+1)*my+j]   * xb3vals[1] * yb3vals[0] +
                 spline->coeffs[(i+1)*my+j+1] * xb3vals[1] * yb3vals[1] +
                 spline->coeffs[(i+1)*my+j+2] * xb3vals[1] * yb3vals[2] +
                 spline->coeffs[(i+1)*my+j+3] * xb3vals[1] * yb3vals[3] +
                 
                 spline->coeffs[(i+2)*my+j]   * xb3vals[2] * yb3vals[0] +
                 spline->coeffs[(i+2)*my+j+1] * xb3vals[2] * yb3vals[1] +
                 spline->coeffs[(i+2)*my+j+2] * xb3vals[2] * yb3vals[2] +
                 spline->coeffs[(i+2)*my+j+3] * xb3vals[2] * yb3vals[3] +

                 spline->coeffs[(i+3)*my+j]   * xb3vals[3] * yb3vals[0] +
                 spline->coeffs[(i+3)*my+j+1] * xb3vals[3] * yb3vals[1] +
                 spline->coeffs[(i+3)*my+j+2] * xb3vals[3] * yb3vals[2] +
                 spline->coeffs[(i+3)*my+j+3] * xb3vals[3] * yb3vals[3]);
        }
    }

    return BS_OK;
}

//-----------------------------------------------------------------------------
// uspline1d
//-----------------------------------------------------------------------------
static void fill_banded_matrix_u(banded_matrix A, int N, double didx,
                                 bs_bctype bctypes[2])
{
    int i, j;

    // values, first and second derivatives of b3_{i-3}, b3_{i-2}, b3_{i-1}
    // at knot i.
    const double b3vals[3] = {1.0/6.0, 2.0/3.0, 1.0/6.0};
    const double db3vals[3] = {-0.5 * didx, 0.0, 0.5 * didx};
    const double d2b3vals[3] = {didx * didx, -2.0 * didx * didx, didx * didx};
    const double notaknot_row[5] = {-1.0, 4.0, -6.0, 4.0, -1.0};

    // Left boundary condition
    switch (bctypes[0]) {
    case BS_DERIV1:
        for (i=0; i<3; i++) A.first[i] = db3vals[i];
        A.first[3] = A.first[4] = 0.0;
        break;
    case BS_DERIV2:
        for (i=0; i<3; i++) A.first[i] = d2b3vals[i];
        A.first[3] = A.first[4] = 0.0;
        break;
    case BS_NOTAKNOT:
        for (i=0; i<5; i++) A.first[i] = notaknot_row[i];
    }

    // rows
    for (i=0; i<N; i++) {
        for (j=0; j<3; j++) A.rows[3*(i+1)+j] = b3vals[j];
    }

    // Right boundary condition
    switch (bctypes[1]) {
    case BS_DERIV1:
        A.last[0] = A.last[1] = 0.0;
        for (i=0; i<3; i++) A.last[i+2] = db3vals[i];
        break;
    case BS_DERIV2:
        A.last[0] = A.last[1] = 0.0;
        for (i=0; i<3; i++) A.last[i+2] = d2b3vals[i];
        break;
    case BS_NOTAKNOT:
        for (i=0; i<5; i++) A.last[i] = notaknot_row[i];
    }
}
    


bs_errorcode bs_uspline1d_create(bs_range x, bs_array y, bs_bcs bcs,
                                 bs_exts exts, bs_uspline1d **out)
{
    bs_uspline1d* spline;
    int i, N, M;
    double didx;
    double* RESTRICT b;
    banded_matrix A;
    bs_bctype bctypes[2] = {bcs.left.type, bcs.right.type};
    
    if (y.size < min_points(bcs.left.type, bcs.right.type))
        return BS_TOOFEWPOINTS;

    spline = malloc(sizeof(bs_uspline1d));

    N = y.size;
    M = N + 2;
    didx = (N - 1) / (x.max - x.min); // equal to  1 / (step size)

    spline->x = x;
    spline->didx = didx;

    spline->n = N;
    spline->exts = exts;

    b = malloc(M * sizeof(double));

    A = alloc_banded_matrix(M);
    
    fill_banded_matrix_u(A, N, didx, bctypes);

    // fill RHS
    b[0] = (bcs.left.type == BS_NOTAKNOT)? 0.0: bcs.left.value;
    b[M-1] = (bcs.right.type == BS_NOTAKNOT)? 0.0: bcs.right.value;
    for (i=0; i<N; i++) {
        b[i+1] = y.data[i * y.stride];
    }

    solve(A, b, M);
    free_banded_matrix(A);
    
    spline->coeffs = b;
    *out = spline;
    return BS_OK;
}


void bs_uspline1d_free(bs_uspline1d* spline)
{
    if (spline != NULL) {
        free(spline->coeffs);
        free(spline);
    }
}


bs_errorcode bs_uspline1d_eval(bs_uspline1d *spline, bs_array x, bs_array out)
{
    int i, j;
    double xval;
    double xfloor;
    double b3vals[4];
    for (j=0; j<x.size; j++) {
        // translate x onto unit basis
        xval = (x.data[j*x.stride] - spline->x.min) * spline->didx;
        xfloor = floor(xval);
        i = (int)xfloor;

        // index outside left boundary
        if (i < 0) {
            switch (spline->exts.left.type) {
            case BS_EXTRAPOLATE:
                i = 0;
                xfloor = 0.0;
                break;
            case BS_CONSTANT:
                i = 0;
                xfloor = 0.0;
                xval = 0.0;
                break;
            case BS_VALUE:
                out.data[j * out.stride] = spline->exts.left.value;
                continue;
            case BS_RAISE:
                return BS_DOMAINERROR;
            }
        }

        // index outside right boundary
        else if (i >= spline->n-1) {
            switch (spline->exts.right.type) {
            case BS_EXTRAPOLATE:
                i = spline->n - 2;
                xfloor = i;
                break;
            case BS_CONSTANT:
                i = spline->n - 2;
                xfloor = i;
                xval = xfloor+1.0;
                break;
            case BS_VALUE:
                out.data[j * out.stride] = spline->exts.right.value;
                continue;
            case BS_RAISE:
                return BS_DOMAINERROR;
            }
        }

        // if we get this far, we're evaluating the spline
        b3unonzeros(xval - xfloor, b3vals);
        out.data[j*out.stride] = (spline->coeffs[i]   * b3vals[0] +
                                  spline->coeffs[i+1] * b3vals[1] +
                                  spline->coeffs[i+2] * b3vals[2] +
                                  spline->coeffs[i+3] * b3vals[3]);
    }

    return BS_OK;
}
