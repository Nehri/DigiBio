import csv
import numpy
import scipy
import sys
from numpy import linalg
import math

def angleSAS(a,b,c):
	return math.acos((pow(b,2)+pow(c,2)-pow(a,2))/(2*b*c))

def main():
	filename = sys.argv[1]
	with open(filename, 'rb') as csvfile:
		reader = csv.reader(csvfile, delimiter=',')
		data = list(reader)

		target = open("method_3_hist_"+filename, 'w+')

        j=0
        while j < len(data):
        	#finds the 3x3 matrix of vectors given 4 points
        	tp1 = data[j]
        	P1 = numpy.array((float(tp1[0]), float(tp1[1]), float(tp1[2])))
        	tp2 = data[j+2]
        	P2 = numpy.array((float(tp2[0]), float(tp2[1]), float(tp2[2])))
        	tp3 = data[j+3]
        	P3 = numpy.array((float(tp3[0]), float(tp3[1]), float(tp3[2])))

        	tp_cen = data[j+1]
        	center = numpy.array((float(tp_cen[0]), float(tp_cen[1]), float(tp_cen[2])))

        	d1_2 = numpy.linalg.norm(P1-P2)
        	d1_3 = numpy.linalg.norm(P1-P3)
        	d2_3 = numpy.linalg.norm(P2-P3)

        	dcen_1 = numpy.linalg.norm(center-P1)
        	dcen_2 = numpy.linalg.norm(center-P2)
        	dcen_3 = numpy.linalg.norm(center-P3)

        	angle1 = angleSAS(d1_2, dcen_1, dcen_2)
        	angle2 = angleSAS(d1_3, dcen_1, dcen_3)
        	angle3 = angleSAS(d2_3, dcen_2, dcen_3)

        	sumAngles = angle1 + angle2 + angle3

        	diffAngles = math.fabs(sumAngles - (2*math.pi))
        	#target.write(data[j][3]+": "+str(diffAngles)+"\n")
		target.write(str(diffAngles)+"\n")

        	j+=4

    	target.close()

if __name__ == '__main__':
	main()