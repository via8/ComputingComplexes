import os
import numpy

from kaucherpy.kaucher import Kaucher as Interval

from core import core

_PATH_DATA = "data"
_PATH_A_INF = os.path.join(_PATH_DATA, "A_inf.txt")
_PATH_A_SUP = os.path.join(_PATH_DATA, "A_sup.txt")
_PATH_B_INF = os.path.join(_PATH_DATA, "b_inf.txt")
_PATH_B_SUP = os.path.join(_PATH_DATA, "b_sup.txt")
_A_INF = numpy.loadtxt(_PATH_A_INF)
_A_SUP = numpy.loadtxt(_PATH_A_SUP)
_B_INF = numpy.loadtxt(_PATH_B_INF)
_B_SUP = numpy.loadtxt(_PATH_B_SUP)
_M = len(_A_INF)
_N = len(_B_INF)
assert(len(_A_SUP) == _M and len(_B_SUP) == _N)
for i in range(0, _N):
    assert(len(_A_INF[i]) == _N and len(_A_SUP[i]) == _N)

_A = numpy.array([
    [
        Interval(_A_INF[i][j], _A_SUP[i][j])
        for j in range(_N)
    ] for i in range(_M)
])
_B = numpy.array([
    Interval(_B_INF[j], _B_SUP[j])
    for j in range(_N)
])


def main():
    initial = core.calc_initial(_A, _B)
    print("initial:")
    print(initial)
    result = core.subdiff_newton(_A, _B, initial)
    print("result:")
    print(result)


if __name__ == "__main__":
    main()
