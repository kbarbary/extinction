cdef extern from "bs.h":
    ctypedef enum bs_errorcode:
        BS_OK           = 0
        BS_OUTOFMEMORY  = 1
        BS_DOMAINERROR  = 2
        BS_NOTMONOTONIC = 3
        BS_SIZEMISMATCH = 4
        BS_BCSIZEMISMATCH = 5
        BS_TOOFEWPOINTS = 6
        
    ctypedef struct bs_array:
        double *data
        int size
        int stride

    ctypedef struct bs_array2d:
        double *data
        int sizes[2]
        int strides[2]

    ctypedef struct bs_range:
        double min
        double max

    ctypedef enum bs_bctype:
        BS_DERIV1
        BS_DERIV2
        BS_NOTAKNOT
        
    ctypedef struct bs_bc:
        bs_bctype type
        double value

    ctypedef struct bs_bcs:
        bs_bc left
        bs_bc right

    ctypedef enum bs_exttype:
        BS_EXTRAPOLATE
        BS_CONSTANT
        BS_VALUE
        BS_RAISE

    ctypedef struct bs_ext:
        bs_exttype type
        double value

    ctypedef struct bs_exts:
        bs_ext left
        bs_ext right
        
    ctypedef struct bs_spline1d:
        double *knots
        double *coeffs
        int n

    ctypedef struct bs_bcarray:
        bs_bctype type
        double *data
        int size
        int stride

    ctypedef struct bs_bcarray_pair:
        bs_bcarray left
        bs_bcarray right
        
    bs_errorcode bs_spline1d_create(bs_array x, bs_array y, bs_bcs bcs, bs_exts exts, bs_spline1d **out)
    bs_errorcode bs_spline1d_eval(bs_spline1d *spline, bs_array x, bs_array out)
    void bs_spline1d_free(bs_spline1d *spline)

    ctypedef struct bs_spline2d:
        double *xknots
        double *xconsts
        double *yknots
        double *yconsts
        double *coeffs
        int nx
        int ny
        bs_exts xexts
        bs_exts yexts

    bs_errorcode bs_spline2d_create(bs_array x, bs_array y, bs_array2d z,
                                    bs_bcarray_pair xbcs, bs_bcarray_pair ybcs,
                                    bs_exts xexts, bs_exts yexts,
                                    bs_spline2d **out)
    bs_errorcode bs_spline2d_eval(bs_spline2d *spline, bs_array x, bs_array y, bs_array2d out)
    void         bs_spline2d_free(bs_spline2d *spline)

    
    ctypedef struct bs_uspline1d:
        bs_range x
        double didx
        double *coeffs
        int n
        bs_exts exts

    bs_errorcode bs_uspline1d_create(bs_range x, bs_array y,
                                     bs_bcs bcs, bs_exts exts,
                                     bs_uspline1d **out)
    bs_errorcode bs_uspline1d_eval(bs_uspline1d *spline, bs_array x,
                                   bs_array out)
    void bs_uspline1d_free(bs_uspline1d *spline)

    
#------------------------------------------------------------------------------
# helpers

class DomainError(Exception):
    """Raised when spline input(s) are outside spline boundaries."""
    pass

cdef int assert_ok(bs_errorcode code) except -1:
    """raise an exception if the error code is not OK"""
    if code == BS_OK:
        return 0
    elif code == BS_OUTOFMEMORY:
        raise MemoryError()
    elif code == BS_DOMAINERROR:
        raise DomainError()
    elif code == BS_NOTMONOTONIC:
        raise ValueError("input array(s) not monotonically increasing")
    elif code == BS_SIZEMISMATCH:
        raise ValueError("input array size mismatch")
    elif code == BS_BCSIZEMISMATCH:
        raise ValueError("boundary condition size mismatch")
    elif code == BS_TOOFEWPOINTS:
        raise ValueError("Too few points in input array. required: 2 + (1 for each not-a-knot condition)")
    else:
        raise Exception("unrecognized error")


cdef inline bs_array to_bs_array(double[:] x):
    return bs_array(&x[0], x.shape[0], x.strides[0]//sizeof(double))


cdef inline bs_array2d to_bs_array2d(double[:, :] x):
    cdef bs_array2d out
    out.data = &x[0, 0]
    out.sizes[0] = x.shape[0]
    out.sizes[1] = x.shape[1]
    out.strides[0] = x.strides[0] // sizeof(double)
    out.strides[1] = x.strides[1] // sizeof(double)
    return out
