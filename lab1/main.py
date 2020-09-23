import numpy
from kaucherpy import Kaucher as interval

def _read_input():
	while(True):
		print("enter float eps from [0, 1) (type 'q' to exit):", end="")
		inp = input().strip()
		if inp == "q" or inp == "Q":
			print("Exiting...")
			exit()

		try:
			eps = float(inp)
		except(Exception):
			print("ERROR: unable to read eps")
			continue

		if not 0 <= eps < 1:
			print("ERROR: specify value for the proper interval")
			continue
		break


	while(True):
		print("enter integer n >= 2 (type 'q' to exit):", end="")
		inp = input()
		if inp == "q" or inp == "Q":
			print("Exiting...")
			exit()

		try:
			n = int(inp)
		except(Exception):
			print("ERROR: unable to read n")
			continue

		if not n >= 2:
			print("ERROR: specify value for the proper interval")
			continue
		break

	return [eps, n]

if __name__ == "__main__":
	eps, n = _read_input()

	# problem 1 (Main Theorem)
	elemDefault = interval(1.0 - eps, 1.0 + eps)
	elemSpecial = interval(1.1 - eps, 1.1 + eps)
	resultInterval = elemDefault * elemDefault - elemSpecial * elemDefault
	print("problem 1:")
	print("eps = " + str(eps))
	print("result determinant interval:" + str(resultInterval))
	if resultInterval.lower <= 0 <= resultInterval.upper:
		print(str(resultInterval.lower) + " <= 0 <= " + str(resultInterval.upper))
		print("initial matrix IS special")
	else:
		print("0 is not included")
		print("initial matrix IS NOT special\n\n")

	# problem 2 (Beck)
	elemDiagonal = interval(1.0, 1.0)
	elemNonDiagonal = interval(0.0, eps)

	midDiagonal = (elemDiagonal.lower + elemDiagonal.upper) / 2.0
	midNonDiagonal = (elemNonDiagonal.lower + elemNonDiagonal.upper) / 2.0

	radDiagonal = (elemDiagonal.upper - elemDiagonal.lower) / 2.0
	radNonDiagonal = (elemNonDiagonal.upper - elemNonDiagonal.lower) / 2.0

	midMatrix = numpy.zeros([n, n])
	radMatrix = numpy.zeros([n, n])
	for i in range(n):
		for j in range(n):
			if i == j:
				midMatrix[i][j] = midDiagonal
				radMatrix[i][j] = radDiagonal
			else:
				midMatrix[i][j] = midNonDiagonal
				radMatrix[i][j] = radNonDiagonal

	midInverseMatrix = abs(numpy.linalg.inv(midMatrix))
	eigenvals = numpy.linalg.eigvals(numpy.matmul(midInverseMatrix, radMatrix))
	spectralRadius = max(abs(eigenvals))
	print("problem 2:")
	print("eps = " + str(eps) + " n = " + str(n))
	print("result spectral radius:" + str(spectralRadius))
	if spectralRadius >= 1.0:
		print(str(spectralRadius) + " >= 1")
		print("initial matrix IS special")
	else:
		print(str(spectralRadius) + " < 1")
		print("initial matrix IS NOT special")

