import os
import sys
import numpy as np

from ase.calculators.socketio import SocketClient
from ase.io import read
from ase.calculators.calculator import Calculator, all_changes
from copy import deepcopy

sys.path.append('/rds/project/rds-Uqezk0eGY00/dellapia/CALCULATORS/Confining_Potential_Calculator')
from confining_potential_calculator_cnt import ConfiningPotentialMorseCalculator


# Define atoms object
print ("Reading atoms object.")
atoms = read("init.pdb", 0)

# Set ASE calculator #################

morse=np.loadtxt('./morse_15_0.txt') # confining potential paramters
confcalc = ConfiningPotentialMorseCalculator(morse[0]/1000,morse[1],morse[2], morse[3]/1000,5.88,0) # confining potential parameters + diameter

atoms.set_calculator(confcalc)
# Create Client
# unix
port = 11111
host = "cf"
print ("Setting up socket.")
client = SocketClient(unixsocket=host)
print ("Running socket.")
client.run(atoms)
print ("setting up calculator")
