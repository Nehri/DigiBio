import csv
import numpy as np
import scipy
import sys
from numpy import linalg
import math

def dist(vec):
	return math.sqrt((vec[0]**2)+(vec[1]**2)+(vec[2]**2))
'''
def dihedralAngle(vecA, vecB):
	numerator = vecA[0]*vecB[0] + vecA[1]*vecB[1] + vecA[2]*vecB[2]
	denominator = dist(vecA)*dist(vecB)

	return math.degrees(math.acos(numerator/denominator))

def dihedralAngle(vecA, vecB, vecC):
	part1 = np.dot(np.cross(np.cross(vecA, vecB), np.cross(vecB, vecC)), np.linalg.norm(vecB))
	part2 = np.dot(np.cross(vecA, vecB), np.cross(vecB, vecC))

	return np.arctan2(part1, part2)
'''
def dihedralAngle2(vecA, vecB, vecC):
	normalA = np.cross(vecA, vecB)/np.dot(vecA, vecB)
	normalB = np.cross(vecC, vecB)/np.dot(vecC, vecB)

	numerator = np.dot(normalA, normalB)
	denominator = np.linalg.norm(normalA) * np.linalg.norm(normalB)

	return math.degrees(math.acos(-numerator/denominator))

def main():
	filename = sys.argv[1]
	with open(filename, 'rb') as csvfile:
		reader = csv.reader(csvfile, delimiter=',')
		data = list(reader)
		target = open("ang_diff_"+filename, 'w+')

		j=0

		while j < len(data):
			#Takes the next 5 points (N, CA, CB, CG1, CG2)
			tp1 = data[j]
			n = (float(tp1[0]), float(tp1[1]), float(tp1[2]))
			tp2 = data[j+1]
			ca = (float(tp2[0]), float(tp2[1]), float(tp2[2]))
			tp3 = data[j+2]
			cb = (float(tp3[0]), float(tp3[1]), float(tp3[2]))
			tp4 = data[j+3]
			cg1 = (float(tp4[0]), float(tp4[1]), float(tp4[2]))
			tp5 = data[j+4]
			cg2 = (float(tp5[0]), float(tp5[1]), float(tp5[2]))
				
			#Forms the three relevant vectors (CA_N, CB_CG1, and CB_CG2) for
			#calculating the two planes
			vecA = np.array([ca[0]-n[0], ca[1]-n[1], ca[2]-n[2]])
			vecB = np.array([cb[0]-ca[0], cb[1]-ca[1], cb[2]-ca[2]])
			vecC = np.array([cg1[0]-cb[0], cg1[1]-cb[1], cg1[2]-cb[2]])
			vecD = np.array([cg2[0]-cb[0], cg2[1]-cb[1], cg2[2]-cb[2]])

			#Calculates the angle between the planes N_CA_CB and CA_CB_CG1
			angleA = dihedralAngle2(vecA, vecB, vecC)
			#Calculates the angle between the planes N_CA_CB and CA_CB_CG2
			angleB = dihedralAngle2(vecA, vecB, vecD)

			#target.write(data[j][3]+": \n\tAngle (CG1) = "+str(angleA)+"\n\tAngle (CG2) = "+str(angleB)+"\n")
			target.write(str(math.fabs(angleA-angleB))+"\n")

			j+=5
			
		target.close()

if __name__ == '__main__':
	main()