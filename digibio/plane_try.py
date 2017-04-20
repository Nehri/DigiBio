import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import functools

def plane(x, y, params):
    a = params[0]
    b = params[1]
    c = params[2]
    z = a*x + b*y + c
    return z

def distance(params, points):
    result = 0
    for (x,y,z) in points:
        plane_z = plane(x, y, params)
        diff = abs(plane_z - z)
        result += diff**2
    return result

points = [(1.1,2.1,8.1),
          (3.2,4.2,8.0),
          (5.3,1.3,8.2),
          (3.4,2.4,8.3)]

fun = functools.partial(distance, points=points)
params0 = [0, 0, 0]
res = scipy.optimize.minimize(fun, params0)

print(distance(res.x,points))
