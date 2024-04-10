import ctypes as ct
from numpy.ctypeslib import ndpointer
import numpy as np


def add_one(number):
    return number + 1


fort_lib = ct.CDLL('_libggcm.cpython-311-darwin.so')

fort_lib.calculate_sum.argtypes = [
    ndpointer(dtype=ct.c_float),
    ct.POINTER(ct.c_float),
    ct.c_int]
fort_lib.calculate_sum.restype = None


def my_sum(arr):
    arr = np.ascontiguousarray(arr, dtype=ct.c_float)
    sum = ct.c_float(0.0)
    n = arr.size
    fort_lib.calculate_sum(arr, sum, n)
    return sum.value


# print(my_mean([1, 2, 3, 4, 5]))
