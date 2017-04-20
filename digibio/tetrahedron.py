import csv
import numpy
import scipy
import sys
from numpy import linalg
import math

def main():
	filename = sys.argv[1]
	with open(filename, 'rb') as csvfile:
		reader = csv.reader(csvfile, delimiter=',')
		data = list(reader)

		target = open("method_2_hist_"+filename, 'w+')

		j=0
		while j < len(data):
        	#finds the 3x3 matrix of vectors given 4 points
			matrix = [[float(data[j][0])-float(data[j+3][0]), float(data[j+1][0])-float(data[j+3][0]), float(data[j+2][0])-float(data[j+3][0])], [float(data[j][1])-float(data[j+3][1]), float(data[j+1][1])-float(data[j+3][1]), float(data[j+2][1])-float(data[j+3][1])],[float(data[j][2])-float(data[j+3][2]), float(data[j+1][2])-float(data[j+3][2]), float(data[j+2][2])-float(data[j+3][2])]]
			
			volume = (1.0/6.0)*math.fabs(linalg.det(matrix))

			#target.write(data[j][3]+": "+str(volume)+"\n")
			target.write(str(volume)+"\n")

			j+=4

		target.close()

if __name__ == '__main__':
	main()