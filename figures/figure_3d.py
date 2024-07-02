#!/usr/bin/env python
# coding: utf-8


import matplotlib.pyplot as plt
import numpy as np
blue='#2571ff'
red='#ff0232'
orange='#ffb927'
green='#008744'
green2='#00ee00'

path_to_data='../data/hbonds'

def  get_error_bootstrap(data):
    from scipy.stats import bootstrap 
    data = (data,)  # samples must be in a sequence
    rng = np.random.default_rng()
    res = bootstrap(data, np.std, confidence_level=0.9,random_state=rng)
    return res.confidence_interval[1]

### load data 

hb_12_320k=np.loadtxt(path_to_data+'/12,0/HBONDS_NVT_320K.txt')
hb_13_320k=np.loadtxt(path_to_data+'/13,0/HBONDS_NVT_320K.txt')
hb_14_320k=np.loadtxt(path_to_data+'/14,0/HBONDS_NVT_320K.txt')
hb_15_320k=np.loadtxt(path_to_data+'/15,0/HBONDS_NVT_320K.txt')
hb_16_320k=np.loadtxt(path_to_data+'/16,0/HBONDS_NVT_320K.txt')

###### average data
n=1500

av_12_320k=[2*np.mean(hb_12_320k[-n:,1]),2*get_error_bootstrap(hb_12_320k[-n:,1])]
av_13_320k=[2*np.mean(hb_13_320k[-n:,1]),2*get_error_bootstrap(hb_13_320k[-n:,1])]
av_14_320k=[2*np.mean(hb_14_320k[-n:,1]),2*get_error_bootstrap(hb_14_320k[-n:,1])]
av_15_320k=[2*np.mean(hb_15_320k[-n:,1]),2*get_error_bootstrap(hb_15_320k[-n:,1])]
av_16_320k=[2*np.mean(hb_16_320k[-n:,1]),2*get_error_bootstrap(hb_16_320k[-n:,1])]

data_320=[]
data_320.append([4.74,av_12_320k[0],av_12_320k[1]])
data_320.append([5.1,av_13_320k[0],av_13_320k[1]])
data_320.append([5.54,av_14_320k[0],av_14_320k[1]])
data_320.append([5.88,av_15_320k[0],av_15_320k[1]])
data_320.append([6.25,av_16_320k[0],av_16_320k[1]])
data_320=np.array(data_320)


fig,ax= plt.subplots(figsize=(7,7))

ax.errorbar(data_320[:,0]*2,data_320[:,1],2*data_320[:,2],ls='--',lw=1.5,label='300K',c='black')

ax.plot(data_320[0,0]*2,data_320[0,1],marker='^',lw=1.5,ms=16,c='black',mfc=blue,mec='black',mew=2.)
ax.plot(data_320[1,0]*2,data_320[1,1],marker='s',lw=1.5,ms=16,c='black',mfc=red,mec='black',mew=2.)
ax.plot(data_320[2,0]*2,data_320[2,1],marker='p',lw=1.5,ms=16,c='black',mfc=orange,mec='black',mew=2.)
ax.plot(data_320[3,0]*2,data_320[3,1],marker='p',lw=1.5,ms=16,c='black',mfc=green,mec='black',mew=2.)
ax.plot(data_320[4,0]*2,data_320[4,1],marker='h',lw=1.5,ms=16,c='black',mfc=green2,mec='black',mew=2.)



plt.xlim([4.7*2,6.30*2])
plt.yticks(fontsize=28)
ax.tick_params(length=8)
plt.xticks(np.arange(9.4,13.0,0.4),fontsize=28,rotation=45)
plt.ylabel('H-Bonds per water',fontsize=30,fontname='arial',labelpad=8)
plt.xlabel('Diameter [$\mathrm{\AA}$]',fontsize=30,fontname='arial',labelpad=6)
plt.grid(ls='dotted',drawstyle='steps',alpha=.4)
plt.tight_layout()
plt.show()
