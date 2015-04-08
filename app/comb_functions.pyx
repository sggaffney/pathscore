__author__ = 'sgg'

"""
list of patient objects:
patient.n_mutated
patient.is_mutated

could become int array and boolean array

write function get_patient_arrays that converts Patient list to
is_mutated boolean array and n_mutations int array


"""
import numpy as np
cimport numpy as np

# DTYPE_t -> float64
ctypedef np.int_t DTYPE_ti
ctypedef np.float64_t DTYPE_td


def get_pway_likelihood_cython(int G,
                        int pway_size,
                        int n_patients,
                        np.ndarray[DTYPE_ti, ndim=1] n_mutated_array,
                        np.ndarray[DTYPE_ti, ndim=1] is_mutated_array):


    cdef DTYPE_td Gd = np.float64(G)
    cdef np.ndarray[DTYPE_td, ndim=1] prob_array = np.zeros(n_patients,
                                                            dtype=np.float64)
    cdef int patient_no
    cdef DTYPE_ti n_patient
    # cdef DTYPE_tb is_mutated
    cdef double p_no_mut, p
    # iterate over patients
    # for patient in self.pway.patients:
    for patient_no in xrange(n_patients):

        n_patient = n_mutated_array[patient_no]

        if G - pway_size >= n_patient:
            p_no_mut = get_p_no_mutations_cython(Gd, pway_size, n_patient)
            if is_mutated_array[patient_no]>0:
                p = 1 - p_no_mut
            else:
                p = p_no_mut
        else:
            # patient has too many hits to get zero pathway mutations
            p = 1  # observation (of mutation in pathway) certain
        prob_array[patient_no] = np.log(p)
        # prob_list.append(log(p))
    # prob_array = array(prob_list)
    return prob_array.sum()

def get_p_no_mutations_cython(DTYPE_td G, int pway_size, DTYPE_ti n):
    """G is background genome size, x is pathway size, n is number of
    genes mutated in patient."""
    cdef DTYPE_td prob = 1
    cdef int i
    for i in xrange(n):
        prob = prob * (1 - pway_size/(G-i))
    # PREVIOUS prob = exp(math.log(comb(G - x, n, exact=True)) - math.log(comb(G, n, exact=True)))
    return prob
