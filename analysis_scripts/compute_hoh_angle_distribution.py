#!/usr/bin/env python
# coding: utf-8

import matplotlib.pyplot as plt
import numpy as np
from ase.io import read


dist_OH_12 = []
for i in ['triangle_tip4p','triangle']:
    geo = read(i+'.xyz')
    n = geo.get_global_number_of_atoms()
    
    dist_OH = []
    for j in np.arange(0,n,3):
        a = geo.get_angles(np.array([[j+1,j,j+2]]),True)
        dist_OH.append(a[0])
    
    dist_OH_12.append([i,np.array(dist_OH)])

dist_OH_13 = []
for i in ['square_tip4p','square']:
    geo = read(i+'.xyz')
    n = geo.get_global_number_of_atoms()
    
    dist_OH = []
    for j in np.arange(0,n,3):
        a = geo.get_angles(np.array([[j+1,j,j+2]]),True)
        dist_OH.append(a[0])

    dist_OH_13.append([i,np.array(dist_OH)])


dist_OH_14 = []
for i in ['penta_1_tip4p','penta_1']:
    geo = read(i+'.xyz')
    n = geo.get_global_number_of_atoms()
    
    dist_OH = []
    for j in np.arange(0,n,3):
        a = geo.get_angles(np.array([[j+1,j,j+2]]),True)
        dist_OH.append(a[0])

    dist_OH_14.append([i,np.array(dist_OH)])


dist_OH_15 = []
for i in ['penta_2_tip4p','penta_2']:
    geo = read(i+'.xyz')
    n = geo.get_global_number_of_atoms()
    
    dist_OH = []
    for j in np.arange(0,n,3):
        a = geo.get_angles(np.array([[j+1,j,j+2]]),True)
        dist_OH.append(a[0])

    dist_OH_15.append([i,np.array(dist_OH)])


dist_OH_16 = []
for i in ['hexa_tip4p','hexa']:
    geo = read(i+'.xyz')
    n = geo.get_global_number_of_atoms()
    
    dist_OH = []
    for j in np.arange(0,n,3):
        a = geo.get_angles(np.array([[j+1,j,j+2]]),True)
        dist_OH.append(a[0])

    dist_OH_16.append([i,np.array(dist_OH)])


fig ,ax = plt.subplots(1,5,figsize=(20,5))

ax[0].hist(dist_OH_12[1][1], lw=4., bins=15, density=True, label='MLP')
ax[0].axvline(104.52, c='black', ls='--', lw=2.5,label='TIP4P')
ax[0].set_title('ice (3,1) @ d $\sim $ 9.5 $\mathrm{ \AA}$',fontsize=20)

ax[1].hist(dist_OH_13[1][1], lw=4., bins=15, density=True,label='MLP')
ax[1].axvline(104.52, c='black', ls='--', lw=2.5)
ax[1].set_title('ice (4,0) @ d $\sim $ 10.2 $\mathrm{ \AA}$',fontsize=20)

ax[2].hist(dist_OH_14[1][1], lw=4., bins=15, density=True, label='MLP')
ax[2].axvline(104.52, c='black', ls='--', lw=2.5)
ax[2].set_title('ice (5,0) @ d $\sim $ 11.0 $\mathrm{ \AA}$',fontsize=20)

ax[3].hist(dist_OH_15[1][1], lw=4., bins=15, density=True, label='MLP')
ax[3].axvline(104.52, c='black', ls='--', lw=2.5)
ax[3].set_title('ice (5,0) @ d $\sim$ 11.8 $\mathrm{ \AA}$',fontsize=20)

ax[4].hist(dist_OH_16[1][1], lw=4., bins=15, density=True, label='MLP')
ax[4].axvline(104.52, c='black', ls='--', lw=2.5)
ax[4].set_title('ice (6,0) @ d $\sim $ 12.5 $\mathrm{ \AA}$',fontsize=20)

ax[0].set_yticks([])
ax[1].set_yticks([])
ax[2].set_yticks([])
ax[3].set_yticks([])
ax[4].set_yticks([])
ax[0].tick_params(direction='out', length=6, width=1, grid_alpha=0.5, labelsize=20)
ax[1].tick_params(direction='out', length=6, width=1, grid_alpha=0.5, labelsize=20)
ax[2].tick_params(direction='out', length=6, width=1, grid_alpha=0.5, labelsize=20)
ax[3].tick_params(direction='out', length=6, width=1, grid_alpha=0.5, labelsize=20)
ax[4].tick_params(direction='out', length=6, width=1, grid_alpha=0.5, labelsize=20)



ax[0].set_ylabel(r'P($\theta_{\mathrm{HOH}})$ [au]',fontsize=24)
ax[2].set_xlabel(r'$\theta_{\mathrm{HOH}}$ $[\degree]$',fontsize=24)
ax[0].set_xlim([100,114])
ax[1].set_xlim([100,114])
ax[2].set_xlim([100,114])
ax[3].set_xlim([100,114])
ax[4].set_xlim([100,114])
ax[0].set_xticks([101,104,107,110,113])
ax[1].set_xticks([101,104,107,110,113])
ax[2].set_xticks([101,104,107,110,113])
ax[3].set_xticks([101,104,107,110,113])
ax[4].set_xticks([101,104,107,110,113])

ax[0].legend(fontsize=14, edgecolor='black')
plt.tight_layout()
plt.show()
