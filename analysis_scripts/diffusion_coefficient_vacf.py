#!/usr/bin/env python
# coding: utf-8

# # Import necessary packages and definitions

# In[1]:


#get_ipython().run_line_magic('load_ext', 'nb_black')
#get_ipython().run_line_magic('matplotlib', 'widget')

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator
import numpy as np
import os
import pickle
import sys
import MDAnalysis
import scipy
from scipy import stats

import pdb
from tqdm.notebook import tqdm
import math
from ase import io
import pandas as pd


# Let us start by defining some variables which are important for our trajectory.

# In[5]:


#########################################################################################
# THIS CELL NEEDS TO BE ADJUSTED
#########################################################################################
# path to extxyz trajectory file
path = "/scratch/fd385/work/water+carbon/MELTING/self-diffusivity/DIFFUSION_NVT-velocities/16,0/250/velocities.xyz"
# the time between to frames in the file in femtoseconds (here 4fs)
time_between_frames = 4


# Let us also create a general class where we can save the objects and data to. We use MDAnalysis to read in the data and create a so-called universe.

# In[6]:


class CNTSimulation:
    """
    Class to store the trajectory and computed data of a given CNT.
    """

    def __init__(
        self, path_to_trajectory: str, time_between_frames: float, name: str = "cnt"
    ):
        """
        Reads in the trajectory and initialises everything we need.
        Args:
            path_to_trajectory (str): Path to trajectory file.
            time_between_frames (float): Time between saved frames in the trajectory in femtoseconds.
            name (str): Name of the system, for instance 12,0 when looping over many tubes or temperatures.
        """

        self.time_between_frames = time_between_frames
        self.name = name
        # we read in the data via MDAnalysis
        self.universe = MDAnalysis.Universe(path_to_trajectory)
        # checking which elements are in the system
        self.species_in_system = np.unique(self.universe.atoms.names)
        # we also define some attributes
        self.vacf_dict = {}  # to save anything related to vacf computation


# # Read in

# Using our defined class we can read in the trajectory and initialise an instance of the class.

# In[7]:


#########################################################################################
# THIS CELL NEEDS TO BE ADJUSTED (name of the trajectory)
#########################################################################################
simulation = CNTSimulation(
    path_to_trajectory=path, time_between_frames=time_between_frames, name="(14,0)"
)


# That's it already, we can now move on to the computation of the diffusion coefficient based on the MSD.

# # Compute Diffusion coefficient based on MSD

# This involves the following steps:
# 1. Compute the VACF in z direction. We split the trajectory into blocks to compute the statistical error.
# 2. Get the plateau value for each block.
# 3. Do bootstrapping with the computed values to extract mean and std of the distribution which corresponds to our prediction and uncertainty estimate.

# ## 1. Compute VACF in z direction for some blocks

# Here is the function which computes the VACF in z direction for a given number of blocks. There is no need to change anything in there. We discuss the arguments we can pass below.

# In[8]:


def compute_velocity_autocorrelation_z(
    simulation,
    species: list = ["O"],
    correlation_time: float = 5000.0,
    number_of_blocks: int = 2,
    start_frame: int = None,
    end_frame: int = None,
    last_n_frames: int = None,
    frame_frequency: int = 1,
):
    """
    Compute mean squared displacement (MSD) for given species.
    Arguments:
        simulation (CNTSimulation): Instance of our class.
        species (list) : Elements considered for MSD. We do it with only oxygen (for COM diffusion let me know).
        correlation_time (float): Time (in fs) for which we will trace the movement of the atoms.
        number_of_blocks (int): Number of blocks used for block average of MSD.
        start_frame (int) : Start frame for analysis (optional).
        end_frame (int) : End frame for analysis (optional).
        last_n_frames (int): Alternatively to given start and end time we can also just say use the last n frames.
        frame_frequency (int): Take every nth frame only (optional). When used with n_last_frames it takes every ith frame of the last_n_frames.
    Returns:

    """


    # check if all species are a subset of species in system
    if not set(species).issubset(simulation.species_in_system):
        raise KeyNotFound(f"At least on of the species specified is not in the system.")

    # get information about sampling either from given arguments or previously set
    if last_n_frames:
        start_frame = -last_n_frames

    else:
        start_frame = start_frame if start_frame else 0
        end_frame = end_frame if end_frame else simulation.universe.trajectory.n_frames

    frame_frequency = frame_frequency

    # convert correlation_time to correlation_frames taken into account the time between frames and
    # the frame frequency
    number_of_correlation_frames = int(
        correlation_time / simulation.time_between_frames / frame_frequency
    )

    # define universe to be used here (basically a copy)
    tmp_universe = simulation.universe

    # check if correlation time can be obtained with current trajectory:
    number_of_samples = (
        int(last_n_frames / frame_frequency)
        if last_n_frames
        else int((end_frame - start_frame))
    )
    if number_of_correlation_frames >= number_of_samples:
        raise UnphysicalValue(
            f" You want to compute a correlation based on {number_of_correlation_frames} frames."
            f"However, the provided trajectory will only be analysed for {number_of_samples} frames.",
            f" Please adjust your correlation or sampling times or run longer trajectories.",
        )

    # convert list to string for species and select atoms of these species:
    # to get water diffusion we need to look at the movement of the center of mass of the molecule
    # if this is the case we simply need to trace allocate only one thrid of the number of atoms
    selected_species_string = " ".join(species)
    atoms_selected = tmp_universe.select_atoms(f"name {selected_species_string}")

    number_of_tracers = (
        int(len(atoms_selected) / 3)
        if selected_species_string == "O H"
        else len(atoms_selected)
    )

    # allocate array for all velocities of all selected atoms for all frames sampled
    saved_velocities_atoms_selected = np.zeros(
        (number_of_samples, number_of_tracers, 3)
    )

    # allocate array for VACF, length of number_of_correlation_frames
    VACF = np.zeros(number_of_correlation_frames)

    # allocate array for number of samples per correlation frame
    number_of_samples_correlated = np.zeros(number_of_correlation_frames)

    # allocate array for blocks of VACF for statistical error analysis
    VACF_block = np.zeros((number_of_blocks, number_of_correlation_frames))

    # define how many samples are evaluated per block
    number_of_samples_per_block = math.ceil(number_of_samples / number_of_blocks)
    index_current_block_used = 0
    previous_global_number_of_samples_correlated = 0

    # make sure that each block can reach full correlation time
    if number_of_samples_per_block < number_of_correlation_frames:
        raise UnphysicalValue(
            f" Your chosen number of blocks ({number_of_blocks}) is not allowed as:",
            f"samples per block ({number_of_samples_per_block}) < correlation frames {number_of_correlation_frames}.",
            f"Please reduce the number of blocks or run longer trajectories.",
        )
        
    # Loop over trajectory to sample all velocities of selected atoms
    for count_frames, frames in enumerate(
       tqdm(
            (tmp_universe.trajectory[start_frame:end_frame])[::frame_frequency]
        )
    ):

        # compute center of mass velocities of selected atoms, which will be substracted afterwards
        center_of_mass_selected_atoms = atoms_selected.center_of_mass()

        # fill array with positions
        saved_velocities_atoms_selected[count_frames] = (
            atoms_selected.positions - center_of_mass_selected_atoms
        )

    # Loop over saved velocities
    for frame, velocities_per_frames in enumerate(
        tqdm(saved_velocities_atoms_selected)
    ):

        # compute last frame sampled, i.e. usually frame+correlation frames
        last_correlation_frame = frame + number_of_correlation_frames
        if last_correlation_frame > number_of_samples - 1:
            last_correlation_frame = number_of_samples

        # define variable to save how many frames where used for correlation
        number_of_frames_correlated = last_correlation_frame - frame

        # increment which correlation frames were sampled
        number_of_samples_correlated[0:number_of_frames_correlated] += 1

        # compute autocorrelation function per frame, yet in all directions

        VACF_per_frame = np.mean(
            #             np.sum(
            saved_velocities_atoms_selected[frame, :, 2]
            * saved_velocities_atoms_selected[frame:last_correlation_frame, :, 2],
            #                 axis=2,
            #             ),
            axis=1,
        )

        VACF[0:number_of_frames_correlated] += VACF_per_frame
        # to get insight on the statistical error we compute block averages
        VACF_block[
            index_current_block_used, 0:number_of_frames_correlated
        ] += VACF_per_frame

        # close block when number of samples per block are reached
        if (
            frame + 1 >= (index_current_block_used + 1) * number_of_samples_per_block
            or frame + 1 == number_of_samples
        ):
            # initialise with 0
            number_of_samples_correlated_per_block = 0
            # check how many samples per frame were taken for this block
            if index_current_block_used == 0:
                # in first block this corresponds to the global number of samples correlated
                number_of_samples_correlated_per_block = number_of_samples_correlated
            else:

                # in all others we just need to get the difference between current and previous global samples
                number_of_samples_correlated_per_block = (
                    number_of_samples_correlated
                    - previous_global_number_of_samples_correlated
                )

            # average current block
            VACF_block[index_current_block_used, :] = (
                VACF_block[index_current_block_used, :]
                / number_of_samples_correlated_per_block
            )

            # define previous global number of samples
            previous_global_number_of_samples_correlated = (
                number_of_samples_correlated.copy()
            )

            # increment index to move to next block
            index_current_block_used += 1

        time = (
            np.arange(number_of_correlation_frames)
            * simulation.time_between_frames
            * frame_frequency
        )

        
        simulation.vacf_dict = {}
        simulation.vacf_dict = {
            "time": time,
            "vacf": VACF / number_of_samples_correlated,
            "vacf_blocks": VACF_block,
        }


# Let's define some parameters for the computation of the VACF. These include:
# - **Correlation time**: For how long would we like to compute the VACF for, i.e. what is the longest time interval we are looking at. 10 ps should be more than sufficient and usually (bulk water) up to 5 ps is used for fitting the diffusion coefficient.
# - **Number of blocks**: In how many blocks would we like to split the trajectory in the computation. The more the better for the error estimate, however, we will later use bootstrapping to get an even better estimate based on these plots so I use usually between 2 and 4. Generally, the fewer blocks the more accurate the respective results.

# In[9]:


#########################################################################################
# THIS CELL NEEDS TO BE ADJUSTED (name of the trajectory)
#########################################################################################
corr_time = 10000
number_of_blocks = 2
start_frame = 0  # we start with the first frame, by not defining the end frame it will just go from start to end


# With these parameters we can run the function to compute the MSD. Per system there will be to loading bars to first substract the center of mass for all frames and compute the MSD.

# In[10]:


compute_velocity_autocorrelation_z(
            simulation=simulation,
            species=["O"],
            correlation_time=corr_time,
            number_of_blocks=number_of_blocks,
            start_frame= start_frame
        )


# The results are saved to the attribute *vacf_dict* which is a dictionary with the following keys:
# - *time*: An array with the time intervals of the computation.
# - *vacf*: The vacf computed from the whole trajectory (I usually don't report this value).
# - *vacf_blocks*: The vacf computed for each block individually. I usually use these computed values to do bootstrapping and get a reliable estimate of the mean and std of the msd. This will be done in section 2.
# 
# We can now plot the VACF as a function of time.

# In[17]:


fig, ax = plt.subplots(1, 1)

ax.plot(
    simulation.vacf_dict["time"] / 1000,
    simulation.vacf_dict["vacf"],
    label="Average over whole trajectory",
)
for i in np.arange(number_of_blocks):
    ax.plot(
        simulation.vacf_dict["time"] / 1000,
        simulation.vacf_dict["vacf_blocks"][i],
        label=f"Block {i}",
    )


ax.legend()
ax.set_xlabel(r"Time [ps]")
ax.set_ylabel(r"VACF")
plt.tight_layout()
plt.savefig('vacf_1.pdf')


# I would say this looks already fairly good. However, let's also have a look at what the integrated vacfoefficient looks like by integrating over the VACF. We assume that the velocity is in atomic units and the conversion factor might need adjusting.

# In[18]:


#########################################################################################
# THIS CELL NEEDS TO BE ADJUSTED (if not a.u.)
#########################################################################################

# conversion factor from atomic units to SI
velocity_unit_conversion = 5.2917721e-11 / 2.4188843e-17


# In[19]:


fig, ax = plt.subplots(1, 1)

# blocks
int_vacf_blocks = scipy.integrate.cumtrapz(
        simulation.vacf_dict["vacf_blocks"] * velocity_unit_conversion ** 2,
        dx=1 * simulation.time_between_frames * 1e-15,
        axis=1,
        initial=0.0,
    )

# average
int_vacf = scipy.integrate.cumtrapz(
        simulation.vacf_dict["vacf"] * velocity_unit_conversion ** 2,
        dx=1 * simulation.time_between_frames * 1e-15,
        initial=0.0,
    )

ax.plot(
    simulation.vacf_dict["time"] / 1000,
    int_vacf,
    label="Average over whole trajectory",
)
for i in range(len(int_vacf_blocks)):
    ax.plot(
        simulation.vacf_dict["time"] / 1000,
        int_vacf_blocks[i],
        label=f"Block {i}",
    )

ax.legend()
ax.set_xlabel(r"Time [ps]")
ax.set_ylabel(r"Integrated VACF [$m^2/s$]")


# Okay this looks fairly good. Unsurprisingly, with increasing time the uncertainty rises. There exist several strategies on how to extract the diffusion coefficient, such as fitting an exponential to the entire curve and taking the plateau value or fitting a straight line to parts of the curve. Here, we will do the latter, but we can discuss alternatives.

# ## 2 Extract diffusion coefficient

# The only parameters we need to specify for our line fit is the range in time we want to fit to. As we said above, it looks like between 5 and 10 ps looks good. However, you can also explore other ranges even though with longer time the error and uncertainty gets larger. So, play a bit around with it.

# In[20]:


#########################################################################################
# THIS CELL NEEDS TO BE ADJUSTED (times)
#########################################################################################
# define these values in femtoseconds here
lower_time = 5000
upper_time = 10000


# Now we compute the diffusion based on the following function:

# In[21]:


from sklearn.metrics import r2_score


def compute_diffusion_coefficient_based_on_VACF(
    measured_int_vacf: np.array,
    measured_time: np.array,
    start_time_fit: float = None,
    end_time_fit: float = None,
    dimension: int = 1,
):
    """
    Computes diffusion coefficient based on MSD.
    Args:
        measured_int_vacf (np.array): The array of the computed MSD.
        measured_time (np.array): The array with the corresponding times where the MSD was measured at. In femtoseconds!
        start_time_fit (float): Start time used for the fit. In femtoseconds!
        end_time_fit (float): End time used for the fit. In femtoseconds!
        dimension (int): MSD was computed for how many dimensions. Used for normalisation.
    """

    # determine interval between measurements in fs:
    time_interval_measurements = measured_time[1] - measured_time[0]

    # determine start and end time for fit:
    start_time_fit = start_time_fit if start_time_fit else (0.2 * measured_time[-1])
    end_time_fit = end_time_fit if end_time_fit else measured_time[-1]

    # convert times into frames
    start_frame_fit = math.ceil(start_time_fit / time_interval_measurements)
    end_frame_fit = int(end_time_fit / time_interval_measurements + 1)

    # get the mean within the range
    mean_for_range = np.mean(
        measured_int_vacf[start_frame_fit:end_frame_fit],
    )

    # check the std score
    std = np.std(measured_int_vacf[start_frame_fit:end_frame_fit])

    if std > 0.01 * mean_for_range:
        print(
            f"WARNING: The std of {std} is {std/mean_for_range*100:1.2f}% of the measured mean.",
            f"We recommend to increase the run time of the trajectory or change the fitting settings.",
        )

    return mean_for_range


# Let's call it for every block and we save it all to the diffusion dictionary.

# In[22]:


time = simulation.vacf_dict["time"]
# diffusion coefficient for each block
diff_coeff_blocks = np.asarray(
    [
        compute_diffusion_coefficient_based_on_VACF(
            measured_int_vacf=int_vacf_blocks[i],
            measured_time=time,
            start_time_fit=lower_time,
            end_time_fit=upper_time,
        )
        for i in range(len(int_vacf_blocks))
    ]
)
# because we are at it let's also fit it for the average over the whole trajectory
diff_coeff_mean = compute_diffusion_coefficient_based_on_VACF(
    measured_int_vacf=int_vacf,
    measured_time=time,
    start_time_fit=lower_time,
    end_time_fit=upper_time,
)

# let us compare the values
print(f"Block 1: {diff_coeff_blocks[0]:1.2e}")
print(f"Block 2: {diff_coeff_blocks[1]:1.2e}")
print(f"Whole trajectory: {diff_coeff_mean:1.2e}")


# Let's plot the quality of these 'fits'.

# In[23]:


fig, ax = plt.subplots(1, 1)

# average over full trajectory
ax.plot(
    simulation.vacf_dict["time"] / 1000,
    int_vacf,
    label="Average over whole trajectory",
)
ax.hlines(diff_coeff_mean,0,10, ls='--')

for i in range(len(int_vacf_blocks)):
    sc = ax.plot(
        simulation.vacf_dict["time"] / 1000,
        int_vacf_blocks[i],
        label=f"Block {i}",
    )
    ax.hlines(diff_coeff_blocks[i],0,10, ls= '--', color = sc[0].get_color())
    

ax.legend()
ax.set_xlabel(r"Time [ps]")
ax.set_ylabel(r"Integrated VACF [$m^2/s$]")
plt.tight_layout()
plt.savefig('vacf_2.pdf')


# Okay block1 is a bit off, we can discuss these but let's move on for now with the bootstrapping as we did for the MSD.

# It seems like we have a good first estimate and the mean lies between the blocks. Also the fits seem sufficiently accurate as we did not get a warning that the R2 score was too low. So you could either now report the mean from here and take the blocks for the std but I recommend doing a 3rd step and do bootsrapping.

# ## 3 Uncertainty estimation via bootstrapping

# Here is the bootstrap function which resamples from the computed diffusion coefficients a lot of times, 10000 should be sufficient though. It then computes the mean and std of the resampled distribution. To be fair, with only 2 blocks there is little gain in it as there are only 3 different combinations (AA, BB, AB). However, with more blocks this gets better.

# In[24]:


def bootstrap(data, num_resamples=10000):
    # Initialize an array to store resampled means
    resampled_means = np.empty(num_resamples)
    for i in range(num_resamples):
        # Resample with replacement
        resampled_data = np.random.choice(data, size=len(data), replace=True)

        # Compute the mean of the resampled data
        resampled_means[i] = np.mean(resampled_data)

    # Calculate statistics from the resampled means
    mean_mean = np.mean(resampled_means)
    std_mean = np.std(resampled_means)

    return mean_mean, std_mean


# Now, let's call the function

# In[25]:


boots_mean, boots_std = bootstrap(diff_coeff_blocks)


# In[26]:


print(f"Mean diffusion coefficient: {boots_mean:1.2e}")
print(f"Std diffusion coefficient: {boots_std:1.2e}")


# So, we see we get the same value as for the expectation value as above. In the paper I usually report 3 times the std as the error.

# In[27]:


with open('Dz','w') as f:
    f.write("%1.3e  %1.3e\n"%(boots_mean,boots_std))
f.close()


# In[ ]:




