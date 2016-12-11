#ifndef BS_H
#define BS_H

//-----------------------------------------------------------------------------
// Error codes
//-----------------------------------------------------------------------------

typedef enum {
  BS_OK           = 0,
  BS_OUTOFMEMORY  = 1,
  BS_DOMAINERROR  = 2,
  BS_NOTMONOTONIC = 3,
  BS_SIZEMISMATCH = 4,
  BS_BCSIZEMISMATCH = 5,
  BS_TOOFEWPOINTS = 6,
} bs_errorcode;


//-----------------------------------------------------------------------------
// Input data types
//-----------------------------------------------------------------------------

typedef struct {
    double *data;
    int size;
    int stride;
} bs_array;

typedef struct {
    double *data;
    int sizes[2];
    int strides[2];
} bs_array2d;

typedef struct {
    double min; // inclusive
    double max; // inclusive
} bs_range;

//-----------------------------------------------------------------------------
// Boundary conditions
//-----------------------------------------------------------------------------

typedef enum {BS_DERIV1, BS_DERIV2, BS_NOTAKNOT} bs_bctype;

typedef struct {
  bs_bctype type;
  double value;
} bs_bc;

typedef struct {
  bs_bc left;
  bs_bc right;
} bs_bcs;

// array boundary conditions for spline2d.
typedef struct {
    bs_bctype type;
    double *data;
    int size;
    int stride;
} bs_bcarray;

typedef struct {
    bs_bcarray left;
    bs_bcarray right;
} bs_bcarray_pair;

//-----------------------------------------------------------------------------
// out-of-domain behavior ("extension")
//-----------------------------------------------------------------------------

typedef enum {BS_EXTRAPOLATE, BS_CONSTANT, BS_VALUE, BS_RAISE} bs_exttype;

typedef struct {
  bs_exttype type;
  double value;
} bs_ext;

typedef struct {
  bs_ext left;
  bs_ext right;
} bs_exts;

//-----------------------------------------------------------------------------
// 1-d splines
//-----------------------------------------------------------------------------

typedef struct {
  double *knots;
  double *consts;
  double *coeffs;
  int n;
  bs_exts exts;
} bs_spline1d;

bs_errorcode bs_spline1d_create(bs_array x, bs_array y, bs_bcs bcs,
                                bs_exts exts, bs_spline1d **out);
bs_errorcode bs_spline1d_eval(bs_spline1d *spline, bs_array x, bs_array out);
void         bs_spline1d_free(bs_spline1d *spline);


typedef struct {
    bs_range x;
    double didx;
    double *coeffs;
    int n;
    bs_exts exts;
} bs_uspline1d;

bs_errorcode bs_uspline1d_create(bs_range x, bs_array y,
                                 bs_bcs bcs, bs_exts exts, bs_uspline1d **out);
bs_errorcode bs_uspline1d_eval(bs_uspline1d *spline, bs_array x, bs_array out);
void         bs_uspline1d_free(bs_uspline1d *spline);


//-----------------------------------------------------------------------------
// 2-d splines
//-----------------------------------------------------------------------------

typedef struct {
    double *xknots;
    double *xconsts;
    double *yknots;
    double *yconsts;
    double *coeffs;
    int nx;
    int ny;
    bs_exts xexts;
    bs_exts yexts;
} bs_spline2d;

bs_errorcode bs_spline2d_create(bs_array x, bs_array y, bs_array2d z,
                                bs_bcarray_pair xbcs, bs_bcarray_pair ybcs,
                                bs_exts xexts, bs_exts yexts,
                                bs_spline2d **out);
bs_errorcode bs_spline2d_eval(bs_spline2d *spline, bs_array x, bs_array y, bs_array2d out);
void         bs_spline2d_free(bs_spline2d *spline);


#endif
