import numpy
import sys
from kaucherpy import Kaucher as interval
from itertools import permutations

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
	print("\n")
	print("problem 1:")
	print("eps = " + str(eps))
	print("result determinant interval:" + str(resultInterval))
	if resultInterval.lower <= 0 <= resultInterval.upper:
		print(str(resultInterval.lower) + " <= 0 <= " + str(resultInterval.upper))
		print("initial matrix IS special")
	else:
		print("0 is not included")
		print("initial matrix IS NOT special\n\n")
	print("\n")

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
	print("\n")

	# test
	if len(sys.argv) < 2:
		exit()
	
	if "--test" not in sys.argv and "-T" not in sys.argv:
		exit()

	def count_inversions(arr):
		count = 0
		n = len(arr)
		for i in range(n):
			for j in range(i + 1, n):
				if arr[i] > arr[j]:
					count += 1
		return count

	def determinant_by_definition(matrix):
		if matrix is None or \
			not hasattr(matrix, "__len__") or \
			not hasattr(matrix[0], "__len__") or \
			len(matrix) != len(matrix[0]):
			return None
		n = len(matrix)
		# n! is quite large number
		if n > 7:
			return None
		perms = permutations(range(n), n)
		sum = interval(0, 0)
		for perm in perms:
			#addend = interval(1, 1)
			addend = matrix[0][perm[0]]
			#for i in range(n):
			for i in range(1, n):
				addend *= matrix[i][perm[i]]
				if count_inversions(perm) % 2 == 0:
					sum += addend
				else:
					sum -= addend
		return sum

	A1 = [[elemDefault, elemDefault], [elemSpecial, elemDefault]]
	print("problem 1: testing determinant value...")
	print(determinant_by_definition(A1))
	print("\n")

	A2 = [[interval(0, 0)] * n for i in range(n)]
	for i in range(n):
		for j in range(n):
			if i == j:
				A2[i][j] = elemDiagonal
			else:
				A2[i][j] = elemNonDiagonal

	print("problem 2: testing determinant value...")
	print(determinant_by_definition(A2))
	print("\n")
