#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt

path_to_data='../data/radial_density/'


blue='#2571ff'
red='#ff0232'
orange='#ffb927'
green='#008744'
green2='#00ee00'

average_320k_triangle=np.loadtxt(path_to_data+'12,0/RD_320K_nvt.txt')
average_320k_square=np.loadtxt(path_to_data+'13,0/RD_320K_nvt.txt')
average_320k_penta_1=np.loadtxt(path_to_data+'14,0/RD_320K_nvt.txt')
average_320k_penta_2=np.loadtxt(path_to_data+'15,0/RD_320K_nvt.txt')
average_320k_hexa=np.loadtxt(path_to_data+'16,0/RD_320K_nvt.txt')


from scipy.interpolate import make_interp_spline

x=average_320k_triangle[:,0]
y=average_320k_triangle[:,1]
triangle_spline = make_interp_spline(x, y)

X_triangle = np.linspace(0, x.max(), 500)
Y_triangle = triangle_spline(X_triangle)
########################################################
x=average_320k_square[:,0]
y=average_320k_square[:,1]
square_spline = make_interp_spline(x, y)
 

X_square = np.linspace(0, x.max(), 500)
Y_square = square_spline(X_square)
########################################################

x=average_320k_penta_1[:,0]
y=average_320k_penta_1[:,1]
penta1_spline = make_interp_spline(x, y)
 

X_penta_1 = np.linspace(0, x.max(), 500)
Y_penta_1 = penta1_spline(X_penta_1)
########################################################

x=average_320k_penta_2[:,0]
y=average_320k_penta_2[:,1]
penta2_spline = make_interp_spline(x, y)
 

X_penta_2 = np.linspace(0, x.max(), 500)
Y_penta_2 = penta2_spline(X_penta_2)
########################################################

x=average_320k_hexa[2::5,0]
y=average_320k_hexa[2::5,1]
hexa_spline = make_interp_spline(x, y)
 
X_hexa = np.linspace(0, x.max(), 200)
Y_hexa = hexa_spline(X_hexa)


fig, ax = plt.subplots(figsize=(10,5))

plt.plot(X_triangle,Y_triangle,ls='solid',lw=3.5,c=blue)
plt.plot(X_square,Y_square,ls='solid',lw=3.5,c=red)
plt.plot(X_penta_1,Y_penta_1,ls='solid',lw=3.5,c=orange)
plt.plot(X_penta_2,Y_penta_2,ls='solid',lw=3.5,c=green)
plt.plot(X_hexa,Y_hexa,ls='solid',lw=2.5,c=green2)

plt.xlim([0,5])
plt.yticks([])
plt.ylabel(r'$\rho(r)$ [au]',fontsize=24,family='arial',labelpad=6)
plt.xlabel('Radial distance $r$ [$\mathrm{\AA}$]',fontsize=24,family='arial')

x=[0,1.0,2.0,3.0,4.0,5.0]
plt.xticks(x,fontsize=22)
ax.tick_params(direction='in',length=8)

from matplotlib.lines import Line2D
legend_elements = [Line2D([0], [0], marker='',mfc=blue,mec='black',mew=1.5,markersize=8,color=blue, ls='solid',lw=2.5, label=r'd ~ 9.5 $\mathrm{\AA}$'),
                   Line2D([0], [0], marker='',mfc=red,mec='black',mew=1.5,markersize=8,color=red, ls='solid',lw=2.5, label=r'd ~ 10.2 $\mathrm{\AA}$'),
                   Line2D([0], [0], marker='',mfc=orange,mec='black',mew=1.5,markersize=8,color=orange, ls='solid',lw=2.5, label=r'd ~ 11.0 $\mathrm{\AA}$'),
                   Line2D([0], [0], marker='',mfc=green,mec='black',mew=1.5,markersize=8,color=green, ls='solid',lw=2.5, label=r'd ~ 11.8 $\mathrm{\AA}$'),
                   Line2D([0], [0], marker='',mfc=green2,mec='black',mew=1.5,markersize=8,color=green2, ls='solid',lw=2.5, label=r'd ~ 12.5 $\mathrm{\AA}$'),
                  ]

plt.legend(handles=legend_elements,fontsize=20,shadow=False,fancybox=False,edgecolor='black')
plt.tight_layout()
plt.show()
