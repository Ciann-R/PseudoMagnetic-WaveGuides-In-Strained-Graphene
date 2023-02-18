#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#---------------------------Preamble/Constants/Strainparameters--------------------------#

area = 100E-9
sigma_x = 15E-9
sigma_y = 10E-9
Height = 1E-9
FWHM_y = 2 * np.sqrt(2* np.log(2)) * sigma_y
FWHM_x = 2 * np.sqrt(2* np.log(2)) * sigma_x
x0 = 0
x00 = 0
y0 = 0
y00 = 0

#--------------------------------Gridspace/2Dplotspace-----------------------------------#

xl = np.linspace(-area/2, area/2, 250)
yl = np.linspace(-area/2, area/2, 250)
dx = xl[1] - xl[0]
dy = yl[1] - yl[0]

x, y = np.meshgrid(xl, yl)
xx = np.arange(-area/2, area/2 , dx)

PseudoScalarConst = 3
PseudoGaugeConst = -(constants.hbar * 3.37) / (2 * 0.142E-9 * constants.e)

R_x = (x -x0)**2
R_y = (y - y0)**2
R_rs = (y - y00)**2
R_xx = (xx - x0)**2
R_rsxx = (xx - x00)**2

#---------------------------------DisplacementVectors------------------------------------#

def GaussianFold(A, sigma_y, R_y):
    return A * np.exp(-R_y / (2 * sigma_y**2))

def GaussianBubble(A, sigma_x, sigma_y, R_x, R_y):
    return A * np.exp(-(R_y) / (2*sigma_y**2) -(R_x) / (2*sigma_x**2))

#-------------------------------------StrainFields---------------------------------------#

def OutOfPlaneStrain(z, dx, dy):
    gradient = np.gradient(z, dx, dy)
    epsxx = 0.5*gradient[1]*gradient[1]
    epsyy = 0.5*gradient[0]*gradient[0]
    epsxy = 0.5*gradient[0]*gradient[1]
    return epsxx, epsyy, epsxy

def PseudoScalarPotential(const, epsxx, epsyy):
    phi = const * (epsxx + epsyy)
    return phi

def PseudoGaugeField(const, epsxx, epsyy, epsxy):
    Ax = const * (epsxx - epsyy)
    Ay = const * (-2*epsxy)
    return Ax, Ay

def PseudoMagneticField(Ay, Ax, dx, dy):
    B = np.gradient(Ay, dx, dy)[1] - np.gradient(Ax, dx, dy)[0]
    return B

def TwoDPlot(Z, sigma_y, const, dx):
    gradient = np.gradient(Z, dx)
    eps = 0.5*gradient*gradient
    PSP = const*(eps)
    B = np.gradient(PSP, dx)
    return PSP, B

#-----------------------------------CallingFunctions-------------------------------------#

#Strain profiles
GF = GaussianFold(Height, sigma_y, R_y)
GB = GaussianBubble(Height, sigma_x, sigma_y, R_x, R_y)
RS = GaussianFold(Height, sigma_y, R_y) + GaussianFold(Height, sigma_y, R_rs)

#----------------------------------MeshStrainProfiles------------------------------------#

fig = plt.figure(figsize=plt.figaspect(0.35))

ax = fig.add_subplot(1,2,1,projection='3d')
plt.title("{TITLE OF STRAIN}")
ax.scatter(x, y, {FUNCTION}, s=0.15, cmap=cm.viridis)

ax = fig.add_subplot(1,2,2)
plt.title("Contour Of Strain Profile",fontsize=20, **csfont, fontweight='bold', y=1.1)
strain_ = ax.pcolormesh(x, y, {FUNCTION}, cmap = cm.inferno, shading='auto')
ax.set_aspect('equal')
plt.xlabel("x",fontsize=18, **csfont)
plt.ylabel("y",fontsize=18, **csfont)
cbar = plt.colorbar(strain_)
cbar.set_label('nm', rotation=0,fontsize=18, **csfont)

fig.tight_layout()
plt.show()

#---------------------------------------PGFs/PSFs----------------------------------------#

OPS = OutOfPlaneStrain({FUNCTION}, dx, dy)
PSP = PseudoScalarPotential(PseudoScalarConst, OPS[0], OPS[1])
PGF = PseudoGaugeField(PseudoGaugeConst, OPS[0], OPS[1], OPS[2])


Q = 5


fig = plt.figure(figsize=plt.figaspect(0.35))


ax = fig.add_subplot(1,2,1)
plt.title("{TITLE OF FUNCTION}",fontsize=20, **csfont, fontweight='bold', y=1.1)
arrows = ax.quiver (x[::Q, ::Q], y[::Q, ::Q], PGF[0][::Q, ::Q], PGF[1][::Q, ::Q],
color='blue')


ax.set_aspect('equal')
plt.xlabel("x",fontsize=18, **csfont)
plt.ylabel("y",fontsize=18, **csfont)
ax = fig.add_subplot(1,2,2)
plt.title("Pseudo Scalar Field",fontsize=20, **csfont, fontweight='bold', y=1.1)
ScalerF = ax.pcolormesh(x, y, PSP, cmap = cm.plasma, shading='auto')
ax.set_aspect('equal')


plt.xlabel("x",fontsize=18, **csfont)
plt.ylabel("y",fontsize=18, **csfont)


cbar = plt.colorbar(ScalerF)
cbar.set_label('$\phi$(meV)', rotation=90,fontsize=18, **csfont)
fig.tight_layout()
plt.show()
#-----------------------------------------PMFs-------------------------------------------#

{FUNCTION} = {FUNCTION}(Height, sigma_y, R_xx, ...)
PMF = PseudoMagneticField(PGF[1], PGF[0], dx, dy)
TDP = TwoDPlot(GF, sigma_x, PseudoScalarConst, dx)
maxfield = max (abs(PMF.max() ), abs(PMF.min()) )


fig = plt.figure(figsize=plt.figaspect(0.35))


ax = fig.add_subplot(1,2,1)
plt.title("{TITLE OF FUCNTION}",fontsize=20, **csfont, fontweight='bold', y=1.1)
plt.plot(xx, TDP[0] / sigma_y, color='purple', label='(y)/_y')
plt.plot(xx, -TDP[1], color='red', label='B')
plt.xlabel("y",fontsize=18, **csfont)
plt.ylabel("B(T)",fontsize=18, **csfont)
leg = plt.legend(loc='upper right')


ax = fig.add_subplot(1,2,2)
plt.title("{TITLE OF FUCNTION}",fontsize=20, **csfont, fontweight='bold', y=1.1)
PseudoMagF = ax.pcolormesh (x, y, PMF, cmap = cm.seismic, vmin= -maxfield, vmax= maxfield)
ax.set_aspect('equal')
plt.xlabel("x",fontsize=18, **csfont)
plt.ylabel("y",fontsize=18, **csfont)
cbar = plt.colorbar(PseudoMagF)
cbar.set_label('B(T)', rotation=90,fontsize=18, **csfont)


fig.tight_layout()

plt.show()

