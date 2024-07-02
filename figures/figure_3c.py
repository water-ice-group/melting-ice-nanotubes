#!/usr/bin/env python
# coding: utf-8

import matplotlib.pyplot as plt
import numpy as np

path_to_data='../data/diffusion_coefficient'

d_cnt120=np.loadtxt(path_to_data+'/Dz_cnt120_nvt_vacf')
d_cnt130=np.loadtxt(path_to_data+'/Dz_cnt130_nvt_vacf')
d_cnt140=np.loadtxt(path_to_data+'/Dz_cnt140_nvt_vacf')
d_cnt150=np.loadtxt(path_to_data+'/Dz_cnt150_nvt_vacf')
d_cnt160=np.loadtxt(path_to_data+'/Dz_cnt160_nvt_vacf')

d_bulk_nvt=np.loadtxt(path_to_data+'/Dz-bulk-nvt')


blue='#2571ff'
red='#ff0232'
orange='#ffb927'
green='#008744'
green2='#00ee00'

#### diffusion vs diameter

for i in range(len(d_cnt120)):
    if(d_cnt120[i,0]==320):
        d_cnt120_320=[d_cnt120[i,1],d_cnt120[i,2]]
        
for i in range(len(d_cnt130)):
    if(d_cnt130[i,0]==320):
        d_cnt130_320=[d_cnt130[i,1],d_cnt130[i,2]]
        
for i in range(len(d_cnt140)):
    if(d_cnt140[i,0]==320):
        d_cnt140_320=[d_cnt140[i,1],d_cnt140[i,2]]

        
for i in range(len(d_cnt150)):
    if(d_cnt150[i,0]==320):
        d_cnt150_320=[d_cnt150[i,1],d_cnt150[i,2]] 

for i in range(len(d_cnt160)):
    if(d_cnt160[i,0]==320):
        d_cnt160_320=[d_cnt160[i,1],d_cnt160[i,2]]       



diff_320=[]
diff_320.append([4.75*2,d_cnt120_320[0]])
diff_320.append([5.1*2,d_cnt130_320[0]])
diff_320.append([5.5*2,d_cnt140_320[0]])
diff_320.append([5.9*2,d_cnt150_320[0]])
diff_320.append([6.25*2,d_cnt160_320[0]])
diff_320=np.array(diff_320)


fig, ax = plt.subplots(figsize=(7,7))

plt.errorbar(x=4.75*2,y=d_cnt120_320[0]*10.**(9),yerr=2*d_cnt120_320[1]*10.**(9),fmt='^',ms=18,mfc=blue,mec='black',c='black',mew=2)
plt.errorbar(x=5.1*2,y=d_cnt130_320[0]*10.**(9),yerr=2*d_cnt130_320[1]*10.**(9),fmt='s',ms=18,mfc=red,mec='black',c='black',mew=2)
plt.errorbar(x=5.5*2,y=d_cnt140_320[0]*10.**(9),yerr=2*d_cnt140_320[1]*10.**(9),fmt='p',ms=18,mfc=orange,mec='black',c='black',mew=2)
plt.errorbar(x=5.9*2,y=d_cnt150_320[0]*10.**(9),yerr=2*d_cnt150_320[1]*10.**(9),fmt='p',ms=18,mfc=green,mec='black',c='black',mew=2)
plt.errorbar(x=6.25*2,y=d_cnt160_320[0]*10.**(9),yerr=2*d_cnt160_320[1]*10.**(9),fmt='h',ms=18,mfc=green2,mec='black',c='black',label='Ice (6,0)',mew=2)

plt.plot(diff_320[:,0],diff_320[:,1]*10.**(9),c='black',ls='--',lw=1.5,label='320K')

plt.axhline(d_bulk_nvt[-1:,1]*10.**(9),lw=2,c='grey',ls='dashed')

plt.xlim([4.7*2,6.3*2])
plt.xticks(np.arange(9.4,12.8,0.4),fontsize=28,rotation=45)
plt.yticks(fontsize=28)
ax.yaxis.get_offset_text().set_fontsize(20)
ax.tick_params(length=8)
plt.ylabel(r'$D_z$ [$10^{-9} \mathrm{m^2 s^{-1}}$]',fontsize=30,fontname='arial',labelpad=8)
plt.xlabel(r'Diameter [$\mathrm{\AA}$]',fontsize=30,fontname='arial',labelpad=8)
plt.grid(ls='dotted',drawstyle='steps',alpha=.4)
plt.tight_layout()
plt.show()
