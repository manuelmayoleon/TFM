from numpy.core.function_base import linspace
import pandas as pd
import numpy as np
import scipy . stats as ss
import math
import matplotlib . mlab as mlab
import matplotlib . pyplot as plt
from scipy import optimize
import matplotlib
matplotlib.rcParams['text.usetex'] = True
from matplotlib.transforms import (
    Bbox, TransformedBbox, blended_transform_factory)
from sympy import *

al= Symbol("a",positive=True)
eps=Symbol("e",positive=True)
g= (4*al*eps**2.0 +12*(1-al))/((3*al+1)*eps**2.0)
# M = Matrix([[1, 0, 1], [2, -1, 3], [4, 3, 2]])
M=Matrix( [[ -1 + al - ((4*al)/(12))*eps**2.0 , (3*al +1)*eps**2.0/12 ],
     [ ((1+al)/4-g/6)*eps**2 , -(1+al)*eps**2.0/(6*g)] ] )
# print(M)

# lamda=M.eigenvals()
# print(lamda)
# data = list(lamda.items())
# # print(data)
# lamda_array=np.array(data)
# print(lamda_array)

# l = symbols('l')
# p = M.charpoly(l)



# s=solveset(p.as_expr(),l)
# # cojer cada autovalor con .args[0] y .args[1]
# print(s.args[0])

# lam1=s.args[0]
# data = list(s.items())
# # print(data)
# lamda_array=np.array(data)
# print(lamda_array)


def lamda1(a,e):
    return -np.sqrt(-16*(a*e**2.0 - 3*a + 3)*(9*a**3*e**4.0 - 3*a**3*e**6.0 + 72*a**2*e**2.0 - 45*a**2*e**4.0 + 4*a**2*e**6.0 - 144*a*e**2.0 + 39*a*e**4.0 - a*e**6.0 + 72*e**2.0 - 3*e**4.0) + (-48*a**2*e**2.0 + 11*a**2*e**4.0 + 72*a**2 + 48*a*e**2.0 + 4*a*e**4.0 - 144*a + e**4.0 + 72)**2)/(48*(a*e**2.0 - 3*a + 3)) + (48*a**2*e**2.0 - 11*a**2*e**4.0 - 72*a**2 - 48*a*e**2.0 - 4*a*e**4.0 + 144*a - e**4.0 - 72)/(48*(a*e**2.0 - 3*a + 3))
def lamda2(a,e):            
        return np.sqrt(-16*(a*e**2.0 - 3*a + 3)*(9*a**3*e**4.0 - 3*a**3*e**6.0 + 72*a**2*e**2.0 - 45*a**2*e**4.0 + 4*a**2*e**6.0 - 144*a*e**2.0 + 39*a*e**4.0 - a*e**6.0 + 72*e**2.0 - 3*e**4.0) + (-48*a**2*e**2.0 + 11*a**2*e**4.0 + 72*a**2 + 48*a*e**2.0 + 4*a*e**4.0 - 144*a + e**4.0 + 72)**2)/(48*(a*e**2.0 - 3*a + 3)) + (48*a**2*e**2.0 - 11*a**2*e**4.0 - 72*a**2 - 48*a*e**2.0 - 4*a*e**4.0 + 144*a - e**4.0 - 72)/(48*(a*e**2.0 - 3*a + 3))
                
                
alfa=np.linspace(0.0,1.0)

epsilon=0.5

fig, ax = plt.subplots(1,1)

plt.plot(alfa,lamda1(alfa,epsilon),label="$\lambda_1 $")
plt.plot(alfa,lamda2(alfa,epsilon),label="$\lambda_2 $")
plt.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)
plt.xlabel ( r' $\alpha$ ', fontsize=20)
plt.ylabel ( r'  ',fontsize=20)
plt.legend(loc=0,fontsize=20)

plt.show()