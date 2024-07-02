#!/usr/bin/env python
# coding: utf-8


import matplotlib.pyplot as plt
import numpy as np

d_cnt120=np.loadtxt('../data/diffusion_coefficient/Dz_cnt120_nvt_vacf')
d_cnt130=np.loadtxt('../data/diffusion_coefficient/Dz_cnt130_nvt_vacf')
d_cnt140=np.loadtxt('../data/diffusion_coefficient/Dz_cnt140_nvt_vacf')
d_cnt150=np.loadtxt('../data/diffusion_coefficient/Dz_cnt150_nvt_vacf')
d_cnt160=np.loadtxt('../data/diffusion_coefficient/Dz_cnt160_nvt_vacf')

d_bulk_nvt=np.loadtxt('../data/diffusion_coefficient/Dz-bulk-nvt') # bulk simulation

exp=np.loadtxt('../data/diffusion_coefficient/Dexp') # bulk experiment

blue='#2571ff'
red='#ff0232'
orange='#ffb927'
green='#008744'
green2='#00ee00'

def  get_error_bootstrap(data):
    from scipy.stats import bootstrap 
    data = (data,)  # samples must be in a sequence
    rng = np.random.default_rng()
    res = bootstrap(data, np.std, confidence_level=0.9,random_state=rng)
    return res.confidence_interval[1]


radius_120=2.35
radius_130=2.75
radius_140=3.15 
radius_150=3.45
radius_160=3.85


nmol_120=710
nmol_130=720
nmol_140=700
nmol_150=762
nmol_160=768


def get_density(data,radius,nmol):
    volume=np.zeros(len(data))
    rho=np.zeros(len(data))
    mass=18 * 1.6 *nmol 
    volume= (np.pi * radius**2)*data
    rho= mass/volume
    err=get_error_bootstrap(rho)
    return [rho.mean(),err]

def flatten_nested_list(input_list):
    flattened_list = []
    for item in input_list:
        if isinstance(item, list):
            flattened_list.extend(flatten_nested_list(item))
        else:
            flattened_list.append(item)
    return flattened_list


#### analysis vs temperature
## CNT(12,0)
density_120=[]
for temp in [250,260,270,275,280,285,290,300,310,320,340]:
    data=np.loadtxt('../data/density/12,0/length.'+str(temp))
    density_120.append([temp, get_density(data,radius_120,nmol_120)])

flattened_list = [flatten_nested_list(sublist) for sublist in density_120]
density_120=np.array(flattened_list)

## CNT(13,0)
density_130=[]
for temp in [250,260,270,275,280,285,290,300,320,340]:
    data=np.loadtxt('../data/density/13,0/length.'+str(temp))
    density_130.append([temp, get_density(data,radius_130,nmol_130)])

flattened_list = [flatten_nested_list(sublist) for sublist in density_130]
density_130=np.array(flattened_list)

## CNT(14,0)
density_140=[]
for temp in [250,260,270,275,280,285,290,300,320,340]:
    data=np.loadtxt('../data/density/14,0/length.'+str(temp))
    density_140.append([temp, get_density(data,radius_140,nmol_140)])

flattened_list = [flatten_nested_list(sublist) for sublist in density_140]
density_140=np.array(flattened_list)

## CNT(15,0)
density_150=[]
for temp in [250,260,275,280,285,300,310,320,330,340]:
    data=np.loadtxt('../data/density/15,0/length.'+str(temp))
    density_150.append([temp, get_density(data,radius_150,nmol_150)])

flattened_list = [flatten_nested_list(sublist) for sublist in density_150]
density_150=np.array(flattened_list)

## CNT(16,0)
density_160=[]
for temp in [250,260,270,280,290,300,310,320,340]:
    data=np.loadtxt('../data/density/16,0/length.'+str(temp))
    density_160.append([temp, get_density(data,radius_160,nmol_160)])

flattened_list = [flatten_nested_list(sublist) for sublist in density_160]
density_160=np.array(flattened_list)


#### UP density 

#### analysis vs temperature
## CNT(12,0)
updensity_120=[]
for temp in [250,260,270,275,280,285,290,300,310,320,340]:
    data=np.loadtxt('../data/density/12,0/length.'+str(temp))
    updensity_120.append([temp, get_density(data,0.05+radius_120,nmol_120)])

flattened_list = [flatten_nested_list(sublist) for sublist in updensity_120]
updensity_120=np.array(flattened_list)

## CNT(13,0)
updensity_130=[]
for temp in [250,260,270,275,280,285,290,300,320,340]:
    data=np.loadtxt('../data/density/13,0/length.'+str(temp))
    updensity_130.append([temp, get_density(data,0.05+radius_130,nmol_130)])

flattened_list = [flatten_nested_list(sublist) for sublist in updensity_130]
updensity_130=np.array(flattened_list)

## CNT(14,0)
updensity_140=[]
for temp in [250,260,270,275,280,285,290,300,320,340]:
    data=np.loadtxt('../data/density/14,0/length.'+str(temp))
    updensity_140.append([temp, get_density(data,0.05+radius_140,nmol_140)])

flattened_list = [flatten_nested_list(sublist) for sublist in updensity_140]
updensity_140=np.array(flattened_list)

## CNT(15,0)
updensity_150=[]
for temp in [250,260,275,280,285,300,310,320,330,340]:
    data=np.loadtxt('../data/density/15,0/length.'+str(temp))
    updensity_150.append([temp, get_density(data,0.05+radius_150,nmol_150)])

flattened_list = [flatten_nested_list(sublist) for sublist in updensity_150]
updensity_150=np.array(flattened_list)

## CNT(16,0)
updensity_160=[]
for temp in [250,260,270,280,290,300,310,320,340]:
    data=np.loadtxt('../data/density/16,0/length.'+str(temp))
    updensity_160.append([temp, get_density(data,0.05+radius_160,nmol_160)])

flattened_list = [flatten_nested_list(sublist) for sublist in updensity_160]
updensity_160=np.array(flattened_list)#### UP density 

#### analysis vs temperature
## CNT(12,0)
updensity_120=[]
for temp in [250,260,270,275,280,285,290,300,310,320,340]:
    data=np.loadtxt('../data/density/12,0/length.'+str(temp))
    updensity_120.append([temp, get_density(data,0.05+radius_120,nmol_120)])

flattened_list = [flatten_nested_list(sublist) for sublist in updensity_120]
updensity_120=np.array(flattened_list)

## CNT(13,0)
updensity_130=[]
for temp in [250,260,270,275,280,285,290,300,320,340]:
    data=np.loadtxt('../data/density/13,0/length.'+str(temp))
    updensity_130.append([temp, get_density(data,0.05+radius_130,nmol_130)])

flattened_list = [flatten_nested_list(sublist) for sublist in updensity_130]
updensity_130=np.array(flattened_list)

## CNT(14,0)
updensity_140=[]
for temp in [250,260,270,275,280,285,290,300,320,340]:
    data=np.loadtxt('../data/density/14,0/length.'+str(temp))
    updensity_140.append([temp, get_density(data,0.05+radius_140,nmol_140)])

flattened_list = [flatten_nested_list(sublist) for sublist in updensity_140]
updensity_140=np.array(flattened_list)

## CNT(15,0)
updensity_150=[]
for temp in [250,260,275,280,285,300,310,320,330,340]:
    data=np.loadtxt('../data/density/15,0/length.'+str(temp))
    updensity_150.append([temp, get_density(data,0.05+radius_150,nmol_150)])

flattened_list = [flatten_nested_list(sublist) for sublist in updensity_150]
updensity_150=np.array(flattened_list)

## CNT(16,0)
updensity_160=[]
for temp in [250,260,270,280,290,300,310,320,340]:
    data=np.loadtxt('../data/density/16,0/length.'+str(temp))
    updensity_160.append([temp, get_density(data,0.05+radius_160,nmol_160)])

flattened_list = [flatten_nested_list(sublist) for sublist in updensity_160]
updensity_160=np.array(flattened_list)


#### LOWs density 

#### analysis vs temperature
## CNT(12,0)
lowdensity_120=[]
for temp in [250,260,270,275,280,285,290,300,310,320,340]:
    data=np.loadtxt('../data/density/12,0/length.'+str(temp))
    lowdensity_120.append([temp, get_density(data,-0.05+radius_120,nmol_120)])

flattened_list = [flatten_nested_list(sublist) for sublist in lowdensity_120]
lowdensity_120=np.array(flattened_list)

## CNT(13,0)
lowdensity_130=[]
for temp in [250,260,270,275,280,285,290,300,320,340]:
    data=np.loadtxt('../data/density/13,0/length.'+str(temp))
    lowdensity_130.append([temp, get_density(data,-0.05+radius_130,nmol_130)])

flattened_list = [flatten_nested_list(sublist) for sublist in lowdensity_130]
lowdensity_130=np.array(flattened_list)

## CNT(14,0)
lowdensity_140=[]
for temp in [250,260,270,275,280,285,290,300,320,340]:
    data=np.loadtxt('../data/density/14,0/length.'+str(temp))
    lowdensity_140.append([temp, get_density(data,-0.05+radius_140,nmol_140)])

flattened_list = [flatten_nested_list(sublist) for sublist in lowdensity_140]
lowdensity_140=np.array(flattened_list)

## CNT(15,0)
lowdensity_150=[]
for temp in [250,260,275,280,285,300,310,320,330,340]:
    data=np.loadtxt('../data/density/15,0/length.'+str(temp))
    lowdensity_150.append([temp, get_density(data,-0.05+radius_150,nmol_150)])

flattened_list = [flatten_nested_list(sublist) for sublist in lowdensity_150]
lowdensity_150=np.array(flattened_list)

## CNT(16,0)
lowdensity_160=[]
for temp in [250,260,270,280,290,300,310,320,340]:
    data=np.loadtxt('../data/density/16,0/length.'+str(temp))
    lowdensity_160.append([temp, get_density(data,-0.05+radius_160,nmol_160)])

flattened_list = [flatten_nested_list(sublist) for sublist in lowdensity_160]
lowdensity_160=np.array(flattened_list)



fig,ax=plt.subplots(2,1,figsize=(10,10))

ax[1].errorbar(d_bulk_nvt[:,0],10**(9)*d_bulk_nvt[:,1],10**(9)*3*d_bulk_nvt[:,2],fmt='o',ls='dotted',lw=2.5,ms=12,label='bulk$^{MLP}$',color='black',mfc='grey',mec='black',mew=1.5)
ax[1].fill_between(d_bulk_nvt[:,0],10**(9)*d_bulk_nvt[:,1]-10**(9)*3*d_bulk_nvt[:,2],10**(9)*d_bulk_nvt[:,1]+10**(9)*3*d_bulk_nvt[:,2],alpha=.3,color='black')


ax[1].plot(exp[:15,0],exp[:15,1],'o',ls='dotted',lw=2.5,ms=12,label='bulk$^{exp}$',color='gray',mfc='none')

ax[1].errorbar(d_cnt120[:,0],10**9*d_cnt120[:,1],2*d_cnt120[:,2]*10**(9),ls='dashed',lw=3.,fmt='^',ms=10,mfc=blue,label='d $\sim 9.5$ $\mathrm{\AA}$',color=blue,mec='black')
ax[1].fill_between(d_cnt120[:,0],10**9*d_cnt120[:,1]-2*d_cnt120[:,2]*10**(9),10**9*d_cnt120[:,1]+2*d_cnt120[:,2]*10**(9),alpha=.3,color=blue)


ax[1].errorbar(d_cnt130[:,0],10**(9)*d_cnt130[:,1],2*d_cnt130[:,2]*10**(9),ls='dashed',lw=3.,fmt='s',ms=10,mfc=red,label='d $\sim 10.2$ $\mathrm{\AA}$',color=red,mec='black')
ax[1].fill_between(d_cnt130[:,0],10**9*d_cnt130[:,1]-2*d_cnt130[:,2]*10**(9),10**9*d_cnt130[:,1]+2*d_cnt130[:,2]*10**(9),alpha=.3,color=red)


ax[1].errorbar(d_cnt140[:,0],10**(9)*d_cnt140[:,1],2*d_cnt140[:,2]*10**(9),ls='-',lw=3.,fmt='p',ms=10,mfc=orange,label='d $\sim 11.0$ $\mathrm{\AA}$',color=orange,mec='black')
ax[1].fill_between(d_cnt140[:,0],10**9*d_cnt140[:,1]-2*d_cnt140[:,2]*10**(9),10**9*d_cnt140[:,1]+2*d_cnt140[:,2]*10**(9),alpha=.3,color=orange)


ax[1].errorbar(d_cnt150[:,0],10**(9)*d_cnt150[:,1],2*d_cnt150[:,2]*10**(9),ls='-',lw=3.,fmt='p',ms=10,mfc=green,label='d $\sim 11.8$ $\mathrm{\AA}$',color=green,mec='black')
ax[1].fill_between(d_cnt150[:,0],10**9*d_cnt150[:,1]-2*d_cnt150[:,2]*10**(9),10**9*d_cnt150[:,1]+2*d_cnt150[:,2]*10**(9),alpha=.3,color=green)


ax[1].errorbar(d_cnt160[:,0],10**(9)*d_cnt160[:,1],2*d_cnt160[:,2]*10**(9),ls='--',lw=3.,fmt='h',ms=10,mfc=green2,label='d $\sim 12.5$ $\mathrm{\AA}$',color=green2,mec='black')
ax[1].fill_between(d_cnt160[:,0],10**9*d_cnt160[:,1]-2*d_cnt160[:,2]*10**(9),10**9*d_cnt160[:,1]+2*d_cnt160[:,2]*10**(9),alpha=.3,color=green2)


ax[0].errorbar(density_120[:,0],density_120[:,1],density_120[:,2],fmt='^',mec='black',mfc=blue,c=blue,ms=10,ls='dashed',lw=2,capsize=6)
ax[0].fill_between(density_120[:,0], lowdensity_120[:,1], updensity_120[:,1] ,alpha=0.3, facecolor=blue)
ax[0].errorbar(density_130[:,0],density_130[:,1],density_130[:,2],fmt='s',mec='black',mfc=red,c=red,ms=10,ls='dashed',lw=2,capsize=6)
ax[0].fill_between(density_130[:,0], lowdensity_130[:,1], updensity_130[:,1] ,alpha=0.3, facecolor=red)
ax[0].errorbar(density_140[:,0],density_140[:,1],density_140[:,2],fmt='p',mec='black',mfc=orange,c=orange,ms=10,ls='dashed',lw=2,capsize=6)
ax[0].fill_between(density_140[:,0], lowdensity_140[:,1], updensity_140[:,1] ,alpha=0.3, facecolor=orange)
ax[0].errorbar(density_150[:,0],density_150[:,1],density_150[:,2],fmt='p',mec='black',mfc=green,c=green,ms=10,ls='dashed',lw=2,capsize=6)
ax[0].fill_between(density_150[:,0], lowdensity_150[:,1], updensity_150[:,1] ,alpha=0.3, facecolor=green)
ax[0].errorbar(density_160[:,0],density_160[:,1],density_160[:,2],fmt='h',mec='black',mfc=green2,c=green2,ms=10,ls='dashed',lw=2,capsize=6)
ax[0].fill_between(density_160[:,0], lowdensity_160[:,1], updensity_160[:,1] ,alpha=0.3, facecolor=green2)

ax[0].set_xticks([250,260,270,280,290,300,310,320,330,340])
ax[0].set_ylabel('Density [g/cm$^3$]',fontsize=24,labelpad=8)

ax[1].set_xticks(np.arange(240,350,10))
ax[1].set_xlim([245,345])
ax[1].yaxis.get_offset_text().set_fontsize(20)

ax[1].tick_params(length=8,labelsize=20)
ax[1].set_ylabel('$D_z$ [$10^{-9} \mathrm{m^2 s^{-1}}$]',fontsize=24,fontname='arial',labelpad=8)
ax[1].set_xlabel('Temperature [K]',fontsize=24,fontname='arial',labelpad=8)
ax[1].legend(fontsize=17,loc='upper left',fancybox=False,shadow=False,edgecolor='black')

ax[0].grid(ls='dotted',drawstyle='steps',alpha=.4)
ax[1].grid(ls='dotted',drawstyle='steps',alpha=.4)

plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9,top=0.9, wspace=0.4,hspace=0.3)
plt.tight_layout()
fig.subplots_adjust(hspace=0.3)

plt.show()
