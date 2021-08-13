# !/ usr / bin / env python3
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





# r=np.genfromtxt("pos.txt",names=["ry","rz"])
# v = np.genfromtxt("velocidad.txt", names=["vy","vz"])
# rinicial=np.genfromtxt("posiciones_init.txt",names=["ry","rz"])
# vinit=np.genfromtxt("velocidad_init.txt", names=["vy","vz"])
# temp=np.genfromtxt("temperaturas.txt" ,names= ["y","z"])
# tiempo=np.genfromtxt("tiemposdecol.txt",names=["t"])

r= pd.read_csv("pos_0.95.txt",header=None,sep='\s+',names=["ry","rz"])
v =  pd.read_csv("velocidad_0.95.txt" ,header=None,sep='\s+' , names=["vy","vz"])
rinicial= pd.read_csv("posiciones_init.txt",header=None,sep='\s+',names=["ry","rz"])
vinit= pd.read_csv("velocidad_init.txt",header=None,sep='\s+', names=["vy","vz"])
temp= pd.read_csv("temperaturas_0.95_0.50.txt",header=None,sep='\s+' ,names= ["y","z"])
tiempo= pd.read_csv("tiemposdecol_0.95.txt",header=None,sep='\s+',names=["t"])






# Esta parte simplemente sirve para poder plotear a la vez
# el fit y el histograma
num_bins =100


fig , ax = plt.subplots (1 ,1)





plt.xlabel ( r' $v_i$ ', fontsize=20)
plt.ylabel ( r' Frecuencia ',fontsize=20)

plt.title ( r' \textbf {Histograma de la velocidad en el eje y}  ',fontsize=30)
# plt.xlim (1 ,9)

# Hacer el histograma 
n,bins,patches = ax.hist(v['vy'],num_bins,density ='true',facecolor ='C0',edgecolor='white',label='$v_y$ ')
# ,edgecolor='yellow'
n2,bins2,patches2 = ax.hist(v['vz'],num_bins,density ='true',facecolor ='C1',edgecolor='white',alpha=0.8,label='$v_z$ ')
plt.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)



N=1000

def gaussian(x,mu,sig):
    return np.exp(-np.power(x-mu,2.)/(2*sig))/np.sqrt(2*np.pi*sig)
x1 =np.linspace(min(v['vy']),max(v['vy']),N)
x2=np.linspace(min(v['vz']),max(v['vz']),N)
#Hacer el fiting de la ley de potencias

# # def fit_func (x ,a , b ) :
# #     return a * x **( b )

params , params_covariance = optimize.curve_fit(gaussian,bins[1:],n,method='dogbox')

params2 , params_covariance2 = optimize.curve_fit(gaussian,bins2[1:],n2,method='dogbox')
# # print('param')
# # print (params)
plt.plot(x1,gaussian(x1,params[0],params[1]),color='C3',label='Ajuste Gaussiano de $v_y$')
plt.plot(x2,gaussian(x2,params2[0],params2[1]),color='C4',label='Ajuste Gaussiano de $v_z$')
plt.legend(loc=0,fontsize=20)
plt.savefig ('histogram_dist.pdf',format ='pdf')

# def fit_func (x ,a , b ) :
#     return a * x+ b 
ax1 = plt.subplots(1,1)
plt.plot(tiempo["t"],temp["z"],label="$T_z$")
plt.plot(tiempo["t"],temp["y"],label="$T_y$")

plt.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)
plt.xlabel ( r' $t$ ', fontsize=20)
plt.ylabel ( r' $T$ ',fontsize=20)
plt.legend(loc=0,fontsize=20)
# plt.xlim(0,1500)
# plt.ylim(0.0,1.0)


ax2 = plt.subplots(1,1)

plt.plot(r["ry"][0:len(r["ry"])-1],r["rz"][0:len(r["ry"])-1], "o")

plt.plot(rinicial["ry"],rinicial["rz"], "o")
# params , params_covariance = optimize.curve_fit(fit_func,corr['tau'],np.log(corr['corr']),method='dogbox')
# print('param')
# print (params)
# print(params_covariance)
# print(np.exp(params[0]))
# ax1 = plt.subplots(1,1)
# # plt.yscale("log")
# tau=np.linspace(0,2.5,100)
# plt.plot(corr['tau'],np.log(corr['corr']),'o')
# # plt.plot(tau,fit_func(tau,params[0],params[1]))
# plt.xlim(-7500,7500)
# plt.ylim(0.5,1.0)
plt.xlabel ( r' $y$ ', fontsize=20)
plt.ylabel ( r' $z$ ',fontsize=20)
# plt.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)
plt.title ( r' \textbf {Posiciones de las partículas confinadas entre placas}  ',fontsize=30)
plt.savefig ('posiciones.pdf',format ='pdf',dpi=1200)

# # plt.legend(loc=0,fontsize=20)
plt.show()
