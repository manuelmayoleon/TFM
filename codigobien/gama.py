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

alfa_exp=np.array([0.95,0.90,0.85,0.80,0.70,0.60])


temp_1=np.genfromtxt("temperaturas_0.95_0.50.txt" ,names= ["y","z"])
temp_2=np.genfromtxt("temperaturas_0.90_0.50.txt" ,names= ["y","z"])
temp_3=np.genfromtxt("temperaturas_0.85_0.50.txt" ,names= ["y","z"])
temp_4=np.genfromtxt("temperaturas_0.80_0.50.txt" ,names= ["y","z"])
temp_5=np.genfromtxt("temperaturas_0.70_0.50.txt" ,names= ["y","z"])
def gama(a,e):
 return  (4*a*e**2.0 +12*(1-a))/((3*a+1)*e**2.0)

gamma=np.zeros(6);
gamma[0]=temp_1["z"][len(temp_1["z"])-1]/temp_1["y"][len(temp_1["y"])-1]
gamma[1]=temp_2["z"][len(temp_2["z"])-10]/temp_2["y"][len(temp_2["y"])-10]
gamma[2]=temp_3["z"][len(temp_3["z"])-100]/temp_3["y"][len(temp_3["y"])-100]
gamma[3]=temp_4["z"][len(temp_4["z"])-100]/temp_4["y"][len(temp_4["y"])-100]
gamma[4]=temp_5["z"][len(temp_5["z"])-100]/temp_5["y"][len(temp_5["y"])-100]


alfa=np.linspace(0.6,1.0)

alfa_exp=np.array([0.95,0.90,0.85,0.80,0.70,0.60])

print(alfa_exp)


plt.plot(alfa,gama(alfa,0.5),label="$\epsilon =0.5 $")
plt.plot(alfa,gama(alfa,0.2), linestyle='--',label="$\epsilon =0.2 $")
plt.plot(alfa_exp[0], gamma[0],marker="o")
plt.plot(alfa_exp[1], gamma[1],marker="o")
plt.plot(alfa_exp[2], gamma[2],marker="o")
plt.plot(alfa_exp[3], gamma[3],marker="o")
plt.plot(alfa_exp[4], gamma[4],marker="o")
plt.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)
plt.xlabel ( r' $\alpha$ ', fontsize=20)
plt.ylabel ( r' $\gamma$  ',fontsize=20)
plt.legend(loc=0,fontsize=20)

plt.show()