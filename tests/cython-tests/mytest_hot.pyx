from hot cimport eval_hot
from libc.stdio cimport printf


cdef size_t THREAD = 1
cdef double E = 1.0
cdef double mu = 1.0
cdef double srcCellParams = 1.0
cdef double *const VEC = &(srcCellParams)
cdef void *const data = NULL

@cytest
def test_eval_hot_runs():
    assert eval_hot(THREAD, E, mu, VEC, data) == 0.0

