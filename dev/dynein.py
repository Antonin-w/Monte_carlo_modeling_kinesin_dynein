import matplotlib.pyplot as plt
import numpy as np
from random import uniform
from random import randint
import matplotlib.animation as animation

# Creation of the matrices

Nt = 100000*5 # Number of timepoints
x = np.zeros(Nt) # Matrix with positions of the motor on the MT 
t = np.zeros(Nt) # Matrix with timepoints
s = np.zeros(Nt) # Matrix with kinesin head description (1: kinesin binding an ATP molecule
                                                    #    0: kinesin without ATP bound)

# Variables 

delta_t = 2e-4 # length of a time step
step = 8e-9 # taking 8-nm steps
Ktrap = 7e-6 # optical trap stiffness
F0 = 0.7e-11 # stalling force
Kcat0 = 55 # load-dependent rate constant
kB = 1.38064852e-23 # boltzmann constant
T = 300 # temperature (Kelvin)
alpha = 0.3 
beta = 0.7
ATP = 1e-9 # ATP concentration
d0 = 6e-9
Psyn0 = 0.23

# rate constants for binding
Kon1 = 4e5
Kon2 = 4e5
Kon3 = Kon2 / 4
Kon4 = Kon2 / 6

Koff1 = 1e-1
Koff2 = 250e-1
Koff3 = Koff2
Koff4 = Koff3

F = F0 

# probability of unbinding ATP
Poff1 = Koff1*delta_t
Poff2 = Koff2*delta_t
Poff3 = Koff3*delta_t
Poff4 = Koff4*delta_t

# Dynein

def bind_unbind(s, ADP_released: False, F):
    
    # probability of binding ATP
    Pon1 = Kon1*ATP*delta_t 
    Pon2 = Kon2*np.exp((F*d0)/(kB*T))*ATP*delta_t 
    Pon3 = Kon3*np.exp((F*d0)/(kB*T))*ATP*delta_t 
    Pon4 = Kon4*np.exp((F*d0)/(kB*T))*ATP*delta_t 

    if not s in [0, 1, 2, 3]:
        raise Exception("s has to be between 0 and 3")

    p = uniform(0,1)

    if not ADP_released:

        if s == 0:
            return 1 if p <= Pon1 else 0

        elif s == 1:
            if p <= Poff1:
                return 0
            elif p > Poff1 and p <= (Poff1 + Pon2):
                return 2
            else:
                return 1 

        elif s == 2:
            if p <= Poff2:
                return 1
            elif p > Poff2 and p <= (Poff2 + Pon3):
                return 3
            else:
                return 2

        else:
            if p <= Poff3:
                return 2
            elif p > Poff3 and p <= (Poff3 + Pon4):
                return 4
            else:
                return 3 
    
    # If site 1 is empty
    else:

        if s == 0:
            return 1 if p <= Pon1 else 0

        elif s == 1:
            if p <= Poff2:
                return 0
            elif p > Poff2 and p <= (Poff2 + Pon1):
                return 2
            else:
                return 1 

        elif s == 2:
            if p <= Poff3:
                return 1
            elif p > Poff3 and p <= (Poff3 + Pon1):
                return 3
            else:
                return 2

        else:
            if p <= Poff4:
                return 2
            elif p > Poff4 and p <= (Poff4 + Pon1):
                return 4
            else:
                return 3 


# Values Pon & Poff
# 8e-08 0.002027246208393239 0.0005068115520983098 0.0003378743680655399
# 2e-05 0.005 0.005 0.005
# really low

def size_step(s):
    s_to_size = {0: 32e-9, 1: 32e-9, 2: 24e-9, 3: 16e-9, 4: 8e-9}
    return s_to_size[s]

def hydrolysis_step(s):
    # Hydrolysis?
    p = uniform(0,1)
    a = 1 if s > 1 else 1/100
    step = size_step(s)
    Pcat = a*(Kcat0*np.exp((-alpha*F*step)/(kB*T))) * delta_t

    if p <= Pcat:
        # Reverse hydrolysis?
        p = uniform(0,1)
        Psyn = Psyn0*np.exp((beta*F*step)/(kB*T))
        if p <= Psyn:
            return s, x
        else:
            return s-1, x+step
    else:
        return s, x



