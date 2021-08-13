import numpy as np
import pandas as pd
from numba import jit
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['text.usetex'] = True
from matplotlib.transforms import (
    Bbox, TransformedBbox, blended_transform_factory)


# temp= pd.read_csv("temperaturas_0.95_0.50.txt" ,header=None,sep='\s+' ,names= ["y","z"])
# tiempo= pd.read_csv("tiemposdecol_0.95.txt",names=["t"])

x0=1.0
y0=5.0
a=0
# b=int(tiempo["t"][len(tiempo)-1])
b=40000
h=0.2

# print(temp)
# print(tiempo)


def f(t1,t2,t):
    alfa=0.95
    epsilon=0.5
    rho=0.010
    # return np.sqrt(np.pi)*(1+alfa)*epsilon*rho*np.sqrt(t1)*( -(1-alfa)*t1+epsilon**2.0*(-(5*alfa-1)*t1 +(3*alfa+1)*t2 )/12)
    return 2.0*(1+alfa)*epsilon*rho*np.sqrt(t1)*( -(1-alfa)*t1+epsilon**2.0*(-4.0*alfa*t1 +(3.0*alfa+1.0)*t2 )/12)/np.sqrt(np.pi)
    # return 4.0*epsilon**3.0*rho*np.sqrt(tx)*( ty -tx )/(3.0*np.sqrt(np.pi))

def g(t1,t2,t):
    alfa=0.95
    epsilon=0.5
    rho=0.010
    vp=0.001
    # return   2.0*np.sqrt(np.pi)*(1+alfa)*epsilon**3.0*rho*np.sqrt(t1)*( 0.5*(1+alfa)*t1-t2)/3.0 +2.0*vp*t2/epsilon
    return   2.0*(1+alfa)*epsilon**3.0*rho*np.sqrt(t1)*( 0.5*(1+alfa)*t1-t2)/(3.0*np.sqrt(np.pi)) +2.0*vp*t2/epsilon
    
    # return -4.0*epsilon**3.0*rho*np.sqrt(tx)*( ty -tx )/(3.0*np.sqrt(np.pi))

# plt.plot(tiempo["t"],temp["z"],color='C0',label="$T_z$ MD")
# plt.plot(tiempo["t"],temp["y"],color='C1',label="$T_y$ MD")      
plt.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)
plt.xlabel ( r' $t(T_0/m\sigma^2)^{1/2}$ ', fontsize=20)
plt.ylabel ( r' $T$ ',fontsize=20)

plt.title ( r' \textbf {Solución de las ecuaciones para las temperaturas}  ',fontsize=30)


@jit
def runge_kutta_system(f, g, x0, y0, a, b, h):
    t = np.arange(a, b + h, h)
    n = len(t)
    x = np.zeros(n)
    y = np.zeros(n)
    x[0] = x0
    y[0] = y0
    for i in range(n - 1):
        k1 = h * f(x[i], y[i], t[i])
        l1 = h * g(x[i], y[i], t[i])
        k2 = h * f(x[i] + k1 / 2, y[i] + l1 / 2, t[i] + h / 2)
        l2 = h * g(x[i] + k1 / 2, y[i] + l1 / 2, t[i] + h / 2)
        k3 = h * f(x[i] + k2 / 2, y[i] + l2 / 2, t[i] + h / 2)
        l3 = h * g(x[i] + k2 / 2, y[i] + l2 / 2, t[i] + h / 2)
        k4 = h * f(x[i] + k3, y[i] + l3, t[i] + h)
        l4 = h * g(x[i] + k3, y[i] + l3, t[i] + h)
        x[i + 1] = x[i] + (1 / 6) * (k1 + 2 * k2 + 2 * k3 + 2 * k4)
        y[i + 1] = y[i] + (1 / 6) * (l1 + 2 * l2 + 2 * l3 + 2 * l4)
    plt.plot(t, x,color='C6',label='$T_y$ ')
    plt.plot(t, y,color='C8',label='$T_z$')
    plt.legend(loc=0,fontsize=20)
    plt.show()
# np.seterr('raise')
runge_kutta_system(f,g,x0,y0,a,b,h)