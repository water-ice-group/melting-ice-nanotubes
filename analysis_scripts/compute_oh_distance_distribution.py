#!/usr/bin/env python
# coding: utf-8


import MDAnalysis as mda
import matplotlib.pyplot as plt
import numpy as np
from ase.io import read
from MDAnalysis.lib.nsgrid import FastNS


dist_OH_12 = []
for i in ['triangle_tip4p','triangle']:
    u = mda.Universe(i+'.xyz')
    geo = read(i+'.xyz')
    z = geo.cell[2][2]

    dist_OH = []
    for ts in u.trajectory:
        ts.dimensions = [30,30,z,90,90,90]
        gridsearch = FastNS(1.1, ts.positions, box=ts.dimensions, pbc=True)
        a = gridsearch.self_search()
        dist_OH.append(a.get_pair_distances())

    dist_OH_12.append([i,np.concatenate(dist_OH)])


dist_OH_13 = []
for i in ['square_tip4p','square']:
    u = mda.Universe(i+'.xyz')
    geo = read(i+'.xyz')
    z = geo.cell[2][2]

    dist_OH = []
    for ts in u.trajectory:
        ts.dimensions = [30,30,z,90,90,90]
        gridsearch = FastNS(1.1, ts.positions, box=ts.dimensions, pbc=True)
        a = gridsearch.self_search()
        dist_OH.append(a.get_pair_distances())

    dist_OH_13.append([i,np.concatenate(dist_OH)])


dist_OH_14 = []
for i in ['penta_1_tip4p','penta_1']:
    u = mda.Universe(i+'.xyz')
    geo = read(i+'.xyz')
    z = geo.cell[2][2]

    dist_OH = []
    for ts in u.trajectory:
        ts.dimensions = [30,30,z,90,90,90]
        gridsearch = FastNS(1.1, ts.positions, box=ts.dimensions, pbc=True)
        a = gridsearch.self_search()
        dist_OH.append(a.get_pair_distances())

    dist_OH_14.append([i,np.concatenate(dist_OH)])


dist_OH_15 = []
for i in ['penta_2_tip4p','penta_2']:
    u = mda.Universe(i+'.xyz')
    geo = read(i+'.xyz')
    z = geo.cell[2][2]

    dist_OH = []
    for ts in u.trajectory:
        ts.dimensions = [30,30,z,90,90,90]
        gridsearch = FastNS(1.1, ts.positions, box=ts.dimensions, pbc=True)
        a = gridsearch.self_search()
        dist_OH.append(a.get_pair_distances())

    dist_OH_15.append([i,np.concatenate(dist_OH)])


dist_OH_16 = []
for i in ['hexa_tip4p','hexa']:
    u = mda.Universe(i+'.xyz')
    geo = read(i+'.xyz')
    z = geo.cell[2][2]

    dist_OH = []
    for ts in u.trajectory:
        ts.dimensions = [30,30,z,90,90,90]
        gridsearch = FastNS(1.1, ts.positions, box=ts.dimensions, pbc=True)
        a = gridsearch.self_search()
        dist_OH.append(a.get_pair_distances())

    dist_OH_16.append([i,np.concatenate(dist_OH)])


fig ,ax = plt.subplots(1,5,figsize=(20,5))

ax[0].hist(dist_OH_12[1][1], lw=2., bins=20, density=True,label='MLP')
ax[0].axvline(0.9572, c='black', ls='--', lw=2.5,label='TIP4P')
ax[0].set_title('ice (3,1) @ d $\sim $ 9.5 $\mathrm{ \AA}$',fontsize=20)

ax[1].hist(dist_OH_13[1][1], lw=2., bins=20, density=True,label='MLP')
ax[1].axvline(0.9572, c='black', ls='--', lw=2.5)
ax[1].set_title('ice (4,0) @ d $\sim $ 10.2 $\mathrm{ \AA}$',fontsize=20)

ax[2].hist(dist_OH_14[1][1], lw=2., bins=20, density=True,label='MLP')
ax[2].axvline(0.9572, c='black', ls='--', lw=2.5)
ax[2].set_title('ice (5,0) @ d $\sim $ 11.0 $\mathrm{ \AA}$',fontsize=20)

ax[3].hist(dist_OH_15[1][1], lw=2., bins=20, density=True,label='MLP')
ax[3].axvline(0.9572, c='black', ls='--', lw=2.5)
ax[3].set_title('ice (5,0) @ d $\sim$ 11.8 $\mathrm{ \AA}$',fontsize=20)


ax[4].hist(dist_OH_16[1][1], lw=2., bins=20, density=True,label='MLP')
ax[4].axvline(0.9572, c='black', ls='--', lw=2.5)
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



ax[0].set_ylabel('P(d$_{\mathrm{OH}}$) [au]',fontsize=24)
ax[2].set_xlabel('d$_{\mathrm{OH}} [\mathrm{\AA}]$',fontsize=24)
ax[0].set_xlim([0.94,0.995])
ax[1].set_xlim([0.94,0.995])
ax[2].set_xlim([0.94,0.995])
ax[3].set_xlim([0.94,0.995])
ax[4].set_xlim([0.94,0.995])


ax[0].legend(fontsize=12, frameon=True,edgecolor='black')
plt.tight_layout()
plt.show()
