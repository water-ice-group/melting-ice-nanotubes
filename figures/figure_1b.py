#!/usr/bin/env python
# coding: utf-8


import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

path_to_data='../data/melting'


#data

raman=np.loadtxt(path_to_data+'/exp_raman')
photo=np.loadtxt(path_to_data+'/exp_photo')
xrd=np.loadtxt(path_to_data+'/exp_xrd')
tip4p=np.loadtxt(path_to_data+'/md_tip4p')
reaxff=np.loadtxt(path_to_data+'/md_reaxff')

# colors
triangle='#2571ff'
square='#ff0232'
penta_1='#ffb927'
penta_2='#008744'
hexa='#00ee00'

exp ='black'


fig,ax=plt.subplots(figsize=(10,10))

ff1=plt.plot(tip4p[:4,0], tip4p[:4,1],marker='o',ms=10,c='grey',label='MD TIP4P',alpha=0.9,mec='black')
plt.plot(tip4p[4:11,0], tip4p[4:11,1],marker='o',ms=10,c='grey',alpha=0.9,mec='black')
plt.plot(tip4p[11:17,0], tip4p[11:17,1],marker='o',ms=10,c='grey',alpha=0.9,mec='black')
plt.plot(tip4p[17:25,0], tip4p[17:25,1],marker='o',ms=10,c='grey',alpha=0.9,mec='black')
plt.plot(tip4p[25:,0], tip4p[25:,1],marker='o',ms=10,c='grey',alpha=0.9,mec='black')

ff2=plt.plot(reaxff[:,0], reaxff[:,1],marker='*',ms=12,c='grey',ls='none',label='MD ReaxFF',alpha=0.9,mec='black')



p1=plt.errorbar(4.75*2,287.5,yerr=5,fmt='^',ms=14,capsize=10,lw=3.,mec='black',c=triangle,mew=2.)
p2=plt.errorbar(5.1*2,287.5,yerr=5,fmt='s',ms=14,capsize=10,lw=3.,mec='black',c=square,mew=2.)
p3=plt.errorbar(5.5*2,295,yerr=5,fmt='p',ms=14,capsize=10,lw=3.,mec='black',c=penta_1,mew=2.)
p4=plt.errorbar(5.9*2,305,yerr=5,fmt='p',ms=14,capsize=10,lw=3.,mec='black',c=penta_2,mew=2.)
p5=plt.errorbar(6.25*2,295,yerr=5,fmt='h',ms=14,capsize=10,lw=3.,mec='black',c=hexa,mew=2.)


eb_raman=plt.errorbar(raman[:,0]*10, raman[:,1]/2.+raman[:,2]/2., yerr=raman[:,2]-raman[:,1],marker='s',ms=10,c=exp,ls='none',mfc='none',capsize=8,label='Raman',alpha=1)
eb_raman[-1][0].set_linestyle('dotted')



eb_pl=plt.errorbar(photo[:,0], photo[:,1],yerr= [ ( photo[:,1]-photo[:,2]), (photo[:,3]-photo[:,1])  ], marker='D',ms=10,c=exp,ls='none',mfc='none',label='PL',capsize=8,alpha=1)
eb_pl[-1][0].set_linestyle('dashed')

eb_xrd=plt.errorbar(xrd[:-1,0], xrd[:-1,1], xerr=(xrd[:-1,3]-xrd[:-1,2])/2., yerr=xrd[:-1,4], marker='>',ms=10,c=exp,ls='none',mfc='none',capsize=8,label='XRD',alpha=1)
eb_xrd[-1][0].set_linestyle('solid')


plt.fill_between(np.arange(8,16,1),270-5,270+5,ls='dashed',color='cornflowerblue',alpha=.2)



ax.yaxis.set_major_locator(MultipleLocator(20))
ax.yaxis.set_minor_locator(MultipleLocator(10))

ax.tick_params(axis='both', which='major', length=10)
ax.tick_params(axis='both', which='minor', length=8)


plt.ylim([195,460])
plt.xticks(np.arange(9.4,14.2,0.8),fontsize=20)
plt.xlim([9.2,12.65])

ax.tick_params(direction='in',length=8)
plt.yticks(fontsize=20)
plt.legend(fontsize=18,shadow=True)
plt.ylabel('Temperature [K]',fontname='Arial',fontsize=24,labelpad=8)
plt.xlabel('Diameter [$\mathrm{\AA}$]',fontname='Arial',fontsize=24,labelpad=8)
#plt.grid(ls='--',alpha=0.4)
plt.tight_layout()

from matplotlib.legend_handler import HandlerTuple
import matplotlib.lines as mlines
tip4p_line = mlines.Line2D([], [], color='grey', marker='o',
                          markersize=10, label='TIP4P',mec='black',alpha=.9)
reaxff_line = mlines.Line2D([], [], color='none',markerfacecolor='grey',  marker='*',
                          markersize=12, label='ReaxFF',mec='black',alpha=.9)
raman_line = mlines.Line2D([], [], color=exp,mec=exp,markerfacecolor='none',  marker='s',
                          markersize=10, label='Raman',alpha=1,ls='dotted')
pl_line = mlines.Line2D([], [], color=exp,mec=exp,markerfacecolor='none',  marker='D',
                          markersize=10, label='PL',alpha=1,ls='dashed')
xrd_line = mlines.Line2D([], [], color=exp,mec=exp,markerfacecolor='none',  marker='>',
                          markersize=10, label='XRD',alpha=1,ls='solid')


plt.legend([(p1, p2, p3, p4,p5),tip4p_line,reaxff_line,raman_line,xrd_line,pl_line],['MD MLP','MD TIP4P','MD ReaxFF','Raman','XRD','PL'],scatterpoints=1, numpoints=1, handler_map={tuple: HandlerTuple(ndivide=None)},fontsize=24,shadow=False,edgecolor='black',fancybox=False)

plt.grid(ls='dotted',drawstyle='steps',alpha=.4)
plt.show()
