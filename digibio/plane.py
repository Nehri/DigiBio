import csv
import numpy as np
import scipy.optimize
import functools
import sys

def plane(x, y, params):
    a = params[0]
    b = params[1]
    c = params[2]
    z = a*x + b*y + c
    return z

def distance(params, points):
    dist = 0
    for (x,y,z) in points:
        plane_z = plane(x, y, params)
        diff = abs(plane_z - z)
        dist += diff**2
    return dist

def main():
	filename = sys.argv[1]
	with open(filename, 'rb') as csvfile:
		reader = csv.reader(csvfile, delimiter=',')
		data = list(reader)
		
		target = open("method_1_hist_"+filename, 'w+')
		
		j=0
		while j < len(data):
			#calculates the center of mass for each chain
			tp1 = data[j]
        		P1 = (float(tp1[0]), float(tp1[1]), float(tp1[2]))
			tp2 = data[j+1]
			P2 = (float(tp2[0]), float(tp2[1]), float(tp2[2]))
        		tp3 = data[j+2]
        		P3 = (float(tp3[0]), float(tp3[1]), float(tp3[2]))
        		tp4 = data[j+3]
        		P4 = (float(tp4[0]), float(tp4[1]), float(tp4[2]))

			PS = [P1, P2, P3, P4]

			centerX = (P1[0]+P2[0]+P3[0]+P4[0])/4.0
			centerY = (P1[1]+P2[1]+P3[1]+P4[1])/4.0
			centerZ = (P1[2]+P2[2]+P3[2]+P4[2])/4.0
			initialGuess = [centerX, centerY, centerZ]
		
			function = functools.partial(distance, points=PS)
			result = scipy.optimize.minimize(function, initialGuess)

			#target.write(data[j][3]+": "+str(distance(result.x,PS))+"\n")
			target.write(str(distance(result.x,PS))+"\n")


			j+=4
		
		target.close()

if __name__ == '__main__':
	main()
