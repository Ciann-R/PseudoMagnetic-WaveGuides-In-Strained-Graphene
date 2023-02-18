#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d import Axes3D
from scipy import constants
#-----------------------------------Graphene Geometry-------------------------------------#
a = 1.42*10E-10 #Lattice constant
acc = a/np.sqrt(3)
Unitcell_x = np.array([0, a/(2*np.sqrt(3)), a*np.sqrt(3)/2, 2*a/np.sqrt(3) ])
Unitcell_y = np.array([0.5*a, 0, 0, 0.5*a])
Unitcell_s = np.array([1,-1,1,-1])
a1 = np.array([np.sqrt(3)*a, 0])
a2 = np.array([0, a])
#-----------------------------------------------------------------------------------------#
csfont = {'fontname':'Times New Roman'}
fig, ax = plt.subplots()
ax.set_aspect('equal')
ax.set_title("4-sublattice Unit Cell", fontsize=20, **csfont)
plt.xlabel('Armchair Direction', fontsize=17.5, **csfont)
plt.ylabel('Zig Zag Direction', fontsize=17.5, **csfont)
ax.plot(Unitcell_x, Unitcell_y, 'o--',c='purple', markersize=20)
ax.plot(Unitcell_x[1], Unitcell_y[1], 'o--',c='black', markersize=20)
ax.plot(Unitcell_x[3], Unitcell_y[3], 'o--',c='black', markersize=20)
ax.arrow(0, 0, a1[0], a1[1],fc='black',ec='black', width=0.2*10E-10)
ax.arrow(0, 0, a2[0], a2[1],fc='black',ec='black', width=0.2*10E-10)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.show()

