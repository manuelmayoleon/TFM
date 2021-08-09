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
g= (12*(1-al) +(5*al-1)*eps**2.0)/((3*al+1)*eps**2.0)
# M = Matrix([[1, 0, 1], [2, -1, 3], [4, 3, 2]])
M=Matrix( [[ -1 + al - ((5*al-1)/(12))*eps**2.0 , (3*al +1)*eps**2.0/12 ],
     [ ((1+al)/2-g/3)*eps**2 , -(1+al)*eps**2.0/(3*g)] ] )
# print(M)

lamda=M.eigenvals()
data = list(lamda.items())
# print(data)
lamda_array=np.array(data)
# print(lamda_array)
l = symbols('l')
p = M.charpoly(l)



s=solveset(p.as_expr(),l)
# cojer cada autovalor con .args[0] y .args[1]
print(s.args[0])

lam1=s.args[0]
# data = list(s.items())
# # print(data)
# lamda_array=np.array(data)
# print(lamda_array)
print(lam1.replace("sqrt", "np.sqrt"))


def lamda1(a,e):
    return -np.sqrt(-24*(5*a*e**2.0 - 12*a - e**2.0 + 12)**3*(12*a**3*e**4.0 -  5*a**3*e**6.0 + 96*a**2*e**2.0 - 76*a**2*e**4.0 + 11*a**2*e**6.0 -  192*a*e**2.0 + 84*a*e**4.0 - 7*a*e**6.0 + 96*e**2.0 - 20*e**4.0 + e**6.0) +  (2160*a**3*e**2.0 - 1044*a**3*e**4.0 + 185*a**3*e**6.0 - 1728*a**3 - 4752*a**2*e**2.0 + 1212*a**2*e**4.0 -   7*a**2*e**6.0 + 5184*a**2 + 3024*a*e**2.0 - 252*a*e**4.0 + 19*a*e**6.0 - 5184*a - 432*e**2.0 + 84*e**4.0 -   5*e**6.0 + 1728)**2)/(24*(5*a*e**2.0 - 12*a - e**2.0 + 12)**2) + (-2160*a**3*e**2.0 + 1044*a**3*e**4.0 - 185*a**3*e**6.0 + 1728*a**3 + 4752*a**2*e**2.0 - 1212*a**2*e**4.0 + 7*a**2*e**6.0 - 5184*a**2 - 3024*a*e**2.0 + 252*a*e**4.0 - 19*a*e**6.0 + 5184*a + 432*e**2.0 - 84*e**4.0 + 5*e**6.0 - 1728)/(24*(5*a*e**2.0 - 12*a - e**2.0 + 12)**2)
def lamda2(a,e):            
        return np.sqrt(-24*(5*a*e**2.0 - 12*a - e**2.0 + 12)**3*(12*a**3*e**4.0 - 5*a**3*e**6.0 + 96*a**2*e**2.0 - 76*a**2*e**4.0 + 11*a**2*e**6.0 -
                192*a*e**2.0 + 84*a*e**4.0 - 7*a*e**6.0 + 96*e**2.0 - 20*e**4.0 + e**6.0) + (2160*a**3*e**2.0 - 1044*a**3*e**4.0 + 185*a**3*e**6.0 - 1728*a**3 - 
                4752*a**2*e**2.0 + 1212*a**2*e**4.0 - 7*a**2*e**6.0 + 5184*a**2 + 3024*a*e**2.0 - 252*a*e**4.0 + 19*a*e**6.0 - 5184*a - 432*e**2.0 + 84*e**4.0 - 
                5*e**6.0 + 1728)**2)/(24*(5*a*e**2.0 - 12*a - e**2.0 + 12)**2) + (-2160*a**3*e**2.0 + 1044*a**3*e**4.0 - 185*a**3*e**6.0 + 1728*a**3 + 4752*a**2*e**2.0 -
                1212*a**2*e**4.0 + 7*a**2*e**6.0 - 5184*a**2 - 3024*a*e**2.0 + 252*a*e**4.0 - 19*a*e**6.0 + 5184*a + 432*e**2.0 - 84*e**4.0 + 5*e**6.0 - 1728)/(24*(5*a*e**2.0 - 12*a - e**2.0 + 12)**2)
                
                
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