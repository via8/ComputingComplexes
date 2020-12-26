import os
import numpy

from kaucherpy.kaucher import Kaucher as Interval

from core import core

_PATH_DATA = "data"
# _PATH_A_INF = os.path.join(_PATH_DATA, "A_inf.txt")  # task 1
# _PATH_A_SUP = os.path.join(_PATH_DATA, "A_sup.txt")  # task 1
# _PATH_B_INF = os.path.join(_PATH_DATA, "b_inf.txt")  # task 1
# _PATH_B_SUP = os.path.join(_PATH_DATA, "b_sup.txt")  # task 1
_PATH_A_INF = os.path.join(_PATH_DATA, "phi12_A_3.txt")  # task 2
_PATH_A_SUP = os.path.join(_PATH_DATA, "phi12_A_3.txt")  # task 2
_PATH_B_INF = os.path.join(_PATH_DATA, "phi12_b_3.txt")  # task 2
_PATH_B_SUP = os.path.join(_PATH_DATA, "phi12_b_3.txt")  # task 2
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
        Interval(_A_INF[i][j], _A_SUP[i][j])  # task 1
        for j in range(_N)
    ] for i in range(_M)
])
_B = numpy.array([
    # Interval(_B_INF[j], _B_SUP[j])  # task 1
    Interval(_B_INF[j] - 0.15, _B_SUP[j] + 0.15)  # task 2
    for j in range(_N)
])


def main():
    initial = core.calc_initial(_A, _B)
    print("initial:")
    for k in range(len(initial)):
        print(initial[k])
    result = core.subdiff_newton(_A, _B, initial)
    print("result:")
    for k in range(len(result)):
        print(result[k])
    print()

    # task 2 loop
    with open("phi12_result_3.txt", "w") as file:
        for k in range(len(result)):
            lower = min(result[k].lower, result[k].upper)
            upper = max(result[k].lower, result[k].upper)
            file.write("%2.18f %2.18f\n" % (lower, upper))


if __name__ == "__main__":
    main()
