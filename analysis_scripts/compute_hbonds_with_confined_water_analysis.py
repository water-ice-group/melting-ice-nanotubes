#!/usr/bin/env python
# coding: utf-8

get_ipython().run_line_magic('matplotlib', 'widget')
get_ipython().run_line_magic('load_ext', 'nb_black')

import os
import sys
from tqdm.notebook import tqdm

from ase import io
import MDAnalysis as mdanalysis

sys.path.append("/home/fd385/codes/confined-water-analysis-master/")
from confined_water import analysis
from confined_water import hydrogen_bonding


import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import figure
from statistics import mean 



# define path to trajectory directory
path = "./"
# specify topology name if necessary (multiple pdb files)

topology_name = "init" # topology read from in init.pdb

nwat=700 # specify the number of water molecules


# create instance of analysis.Simulation class
simulation = analysis.Simulation(path)

# read in positions based on topology file
simulation.read_in_simulation_data(
    read_positions=True, topology_file_name=topology_name
)

# set preferred sampling times and define time between frames
simulation.set_sampling_times(
    start_time=1.*10**6, end_time=-1, frame_frequency=1, time_between_frames=100
)

# set which directions are periodic
simulation.set_pbc_dimensions("z")

# perform analysis
simulation.set_up_hydrogen_bonding_analysis()

# results are accessible like this:
simulation.hydrogen_bonding[0].dataframe

df = simulation.hydrogen_bonding[0].dataframe

start_time=int(1*10**6)
end_time=2000000 
time_between_frames=100
frame_times = list(range(start_time,end_time, time_between_frames)) #in fs

frame_times_conv = []
for item in frame_times:
    frame_times_conv.append(item /1000)


curr_h_bonds = []
for i in frame_times:
    num_hb = len(simulation.hydrogen_bonding[0].dataframe["Time"][df.Time==i])/nwat
    curr_h_bonds.append(num_hb)

#compute running average 
step=100
running_average=np.zeros(len(curr_h_bonds))

for i in range(len(running_average)):
    if(i<step):
        running_average[i]=curr_h_bonds[i]
    else:
        for k in range(step):
            running_average[i] += curr_h_bonds[i-k]/step


fig = plt.figure()
ax = plt.axes()

ax.plot(frame_times_conv, curr_h_bonds, c='black', linewidth=2, label="instant value")
ax.plot(frame_times_conv, running_average, c='red', linewidth=2, label="running average")
ax.set_xlabel(r'Time [ps]')
ax.set_ylabel(r"$N_{HB}$")
ax.label_outer()

legend_3 = fig.legend(bbox_to_anchor = (1.0, 1.05),loc="upper right", ncol=2)
legend_3.get_frame().set_alpha(0)
legend_3.get_frame().set_facecolor((0, 0, 0, 0))

fig.tight_layout()


import matplotlib.pyplot as plt
import seaborn as sns

# matplotlib histogram
plt.hist(curr_h_bonds, color = 'blue', edgecolor = 'black',
         bins = int(200/5))
plt.xlabel('Number of hydrogen bonds')
plt.ylabel('Count')



################################
####### SAVING HBONDS INFO ######
#################################
file = open("HBONDS.txt", "w+")
for i in range(len(curr_h_bonds)):
    file.write(str(frame_times_conv[i])+'   '+ str(curr_h_bonds[i])+'   '+str(running_average[i])+' \n')
file.close()
