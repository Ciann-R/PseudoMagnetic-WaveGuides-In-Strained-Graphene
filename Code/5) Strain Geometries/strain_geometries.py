#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#-------------------------------------StrainParameters----------------------------------#

A = 10E-10
x0 = 1E-8
y0 = 1.2E-8

y00 = 3.8E-8#2.2E-8

sigma = 0.35E-8 #7.5E-9

cutoff = 20

#Creates array of 0's with identical dimensions as xcoords
GaussianBubble = np.zeros_like(xcoords)
GaussianRipple = np.zeros_like(xcoords)
SinusoidalRipple = np.zeros_like(4*xcoords)

#Gaussian Bubble
R_GB = np.sqrt((xcoords - x0)**2 + (ycoords -y0)**2)

#Gaussian Ripple
R_GR = np.sqrt((ycoords - x0)**2)
R_GRR = np.sqrt((ycoords - y00)**2)

inside_GB = np.where(R_GB < cutoff)
inside_GR = np.where(R_GR < cutoff)
#inside_SR = np.where(R_SR < cutoff)


#------------------------------------StrainFunctions------------------------------------#

GaussianBubble[inside_GB] =((1.0 / (2 * np.pi * sigma**2)) * A * np.exp(- R_GB[inside_GB]**2
                                                                          / (2 * sigma**2)))


GaussianRipple[inside_GR] =((1.0 / (2 * np.pi * sigma**2)) * np.exp(- R_GR[inside_GR]**2
                                                                      / (2 * sigma**2)))


SinusoidalRipple[inside_GR] =((1.0 / (2 * np.pi * sigma**2)) * np.exp(- R_GR[inside_GR]**2
                                                                        / (2 * sigma**2)))
        + ((1.0 / (2 * np.pi * sigma**2)) * np.exp(- R_GRR[inside_GR]**2 / (2 * sigma**2)))


#--------------------------------PlottingStrainGeometry---------------------------------#

fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(projection='3d')
ax.scatter(xcoords, ycoords, {strain-function}, s=5, c='black')
plt.xlabel('X')
plt.ylabel('Y')
# Get rid of colored axes planes
# First remove fill
ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False

# Now set color to white (or whatever is "invisible")
ax.xaxis.pane.set_edgecolor('w')
ax.yaxis.pane.set_edgecolor('w')
ax.zaxis.pane.set_edgecolor('w')

# Bonus: To get rid of the grid as well:
ax.grid(False)

plt.show()

