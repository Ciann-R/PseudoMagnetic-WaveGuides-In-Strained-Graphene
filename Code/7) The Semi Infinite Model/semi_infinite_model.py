#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#--------------------------------------Nano-Ribbon---------------------------------------#

a = 1.42*10E-10 #Lattice constant
acc = a/np.sqrt(3)

Unitcell_x = np.array([0, a/(2*np.sqrt(3)), a*np.sqrt(3)/2, 2*a/np.sqrt(3) ])
Unitcell_y = np.array([0.5*a, 0, 0, 0.5*a])
Unitcell_s = np.array([1,-1,1,-1])

a1 = np.array([np.sqrt(3)*a, 0])
a2 = np.array([0, a])

#-----------------------------------------------------------------------------------------#

coords_left = x - np.sqrt(3)*a
coords_right = x + np.sqrt(3)*a


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


#Sheet dimensions

n_x = 15
n_y = 15
x = np.array([])
y = np.array([])

sublattice = np.array([])

for i in np.arange(n_x):
    for j in np.arange(n_y):
        x = np.concatenate( (x, Unitcell_x + i * a1[0] + j* a2[0]) )
        y = np.concatenate( (y, Unitcell_y + i * a1[1] + j *a2[1]))
        sublattice = np.concatenate ( (sublattice, Unitcell_s))
        
        
fig, ax = plt.subplots()
ax.set_aspect('equal')
ax.set_title("Nano-Ribbon Graphene Model", fontsize=20, **csfont, x= 0.55, y=1.1)
plt.xlabel('Armchair Direction', fontsize=17.5, **csfont)
plt.ylabel('Zig Zag Direction', fontsize=17.5, **csfont)


ax.plot(x[np.where(sublattice==1)], y[np.where(sublattice==1)], 'o',
                                                            markersize=4,c='purple')

ax.plot(x[np.where(sublattice==-1)], y[np.where(sublattice==-1)], 'o',
                                                            markersize=4,c='black')

ax.plot(coords_right[np.where(sublattice==1)], ycoords[np.where(sublattice==1)], 'x', \par
                                                            markersize=4, c = 'purple')

ax.plot(coords_right[np.where(sublattice==-1)], ycoords[np.where(sublattice==-1)], 'x', \par
                                                            markersize=4, c = 'black')

ax.plot(coords_left[np.where(sublattice==1)], ycoords[np.where(sublattice==1)], 'x', \par
                                                            markersize=4, c = 'purple')

ax.plot(coords_left[np.where(sublattice==-1)], ycoords[np.where(sublattice==-1)], 'x', \par
                                                            markersize=4, c = 'black')

centre = np.array((xcoords, ycoords)).T
left = np.array((coords_left, ycoords)).T
right = np.array((coords_right, ycoords)).T


def DistanceMatrix(coords_1, coords_2, acc, t0):
    dist = distance_matrix(coords_1, coords_2)
    hops=np.isclose(dist, acc)
    H = -t0 * hops
    return dist, hops, H


def EigenSoln(H):
    vals, vecs = np.linalg.eig(H)
    el = vals.argsort()
    vals = vals[el]
    vecs = vecs[:,el]
    return vals, vecs


CC = DistanceMatrix(centre, centre, acc, t0)
CL = DistanceMatrix(centre, left, acc, t0)
CR = DistanceMatrix(centre, right, acc, t0)
ES = EigenSoln(CC[2])


fix, ax = plt.subplots(figsize=(3, 8))
[ax.plot( [0, 1], [val, val], 'k-' ) for val in ES[0]]
plt.show()


kx = np.linspace (-3, 3, 2000)

def Hamiltonian(kx, t0, hops, hopsl, hopsr, d):
    H = t0 * hops + t0 * hopsl *np.exp(-1j * kx * d ) + t0 * hopsr *np.exp(1j * kx * d )
    evals=np.linalg.eigvalsh(H)
    return evals

full_evals=[]

for k in kx:
    evals = Hamiltonian(k, t0, CC[1], CL[1], CR[1], np.sqrt(3)*a)
    full_evals.append(evals)
full_evals = np.array(full_evals)

csfont = {'fontname':'Times New Roman'}

plt.title("Electronic Bandstructure ", fontsize=20, **csfont, fontweight='bold', y=1.1)
[plt.plot(kx, full_evals[:, i], 'k') for i in np.arange(len(centre))]
plt.xlabel('Kx',fontsize=18, **csfont, fontweight='bold')
plt.ylabel('E(k)', fontsize=18, **csfont, fontweight='bold')
plt.show()


#------------------------------------StrainParameters--------------------------------------#

kx = np.linspace (-3, 3, 2000)
ky = np.linspace (-3, 3, 2000)
A = 3
B = 3.37
a = 1

sigma_x = 10
sigma_y = 15

#Fold parameters
R_xc = (ycoords - y0)**2
R_xcc = (ycoords - y00)**2
R_xl = (coords_left - x0)**2
R_xr = (coords_right - x0)**2

#--------------------------------AcquiringStrainCoords-------------------------------------#

def GaussianFold(A, sigma_x, R_x):
    return A * np.exp(-R_x / (2 * sigma_x**2))


RS = GaussianFold(A, sigma_x, R_xc) + GaussianFold(A, sigma_x, R_xcc)


strain_centre = np.array([xcoords, ycoords, GaussianFold(A, sigma_x, R_xc)]).T
strain_right = np.array([coords_right, ycoords, GaussianFold(A, sigma_x, R_xc) ]).T
strain_left = np.array([coords_left, ycoords, GaussianFold(A, sigma_x, R_xc) ]).T

strain_centre1 = np.array([xcoords, ycoords, RS]).T
strain_right1 = np.array([coords_right, ycoords, RS ]).T
strain_left1 = np.array([coords_left, ycoords, RS ]).T


#---------------------------------DistanceMatrices/NNHopping-------------------------------#

def StrainDistanceMatrix(coords_1, coords_2, strain_1, strain_2, acc):
    dist = distance_matrix(coords_1, coords_2)
    zdist = distance_matrix(strain_1, strain_2)
    hops=np.isclose(dist, acc)
    return dist, zdist, hops


def StrainHopping(zdist, B, t0, acc):
    return t0 * np.exp(-B *((zdist - acc) / acc))


#----------------------------------------Hamiltonian---------------------------------------#
def StrainHamiltonian(kx, tc, tl, tr, hops, hopsl, hopsr, d):
    H = tc * hops + tl * hopsl *np.exp(-1j * kx * d ) + tr * hopsr *np.exp(1j * kx * d )
    evals = np.linalg.eigvalsh(H)
    return evals

def StrainHamiltonian2(kx, tc, tl, tr, hops, hopsl, hopsr, d):
    H = tc * hops + tl * hopsl *np.exp(-1j * kx * d ) + tr * hopsr *np.exp(1j * kx * d )
    evals, evecs = np.linalg.eigh(H)
    return evals, evecs



#-----------------------------CallingStrainGeometryFunctions------------------------------#

#Repeat for other strain geometries outlined

CS = StrainDistanceMatrix(strain_centre[:, 0:2], strain_centre[:, 0:2],
                                                    strain_centre, strain_centre, acc)

CLS = StrainDistanceMatrix(strain_centre[:, 0:2], strain_left[:, 0:2],
                                                    strain_centre, strain_left, acc)

CRS = StrainDistanceMatrix(strain_centre[:, 0:2], strain_right[:, 0:2],
                                                    strain_centre, strain_right, acc)

tc = StrainHopping(CS[1], B, t0, acc)
tl = StrainHopping(CLS[1], B, t0, acc)
tr = StrainHopping(CRS[1], B, t0, acc)

strain_evals=[]

for k in kx:
    s_evals = StrainHamiltonian(k, tc, tl, tr, CS[2], CLS[2], CRS[2], np.sqrt(3)*a)
    strain_evals.append(s_evals)
    
strain_evals = np.array(strain_evals)



H_s = CS[2]*tc
H_sval = np.linalg.eigvalsh(H_s)


csfont = {'fontname':'Times New Roman'}
fix, ax = plt.subplots(figsize=(3, 8))
plt.ylabel('Energy (eV)', fontsize=17.5, **csfont)
[ax.plot( [0, 1], [val, val], 'k-', c='r' ) for val in H_sval]
[ax.plot( [0, 1], [val, val], 'k-' ) for val in ES[0]]
plt.show()


[plt.plot(kx, strain_evals[:, i], 'k', c='r') for i in np.arange(len(centre))]
[plt.plot(kx, full_evals[:, i], 'k') for i in np.arange(len(centre))]


plt.title("Gaussian Fold", fontsize=20, **csfont, fontweight='bold', y=1.1)
plt.xlabel('Kx',fontsize=18, **csfont, fontweight='bold')
plt.ylabel('E(k)', fontsize=18, **csfont, fontweight='bold')
plt.show()



band = int(len(centre)/2)


kspec = 0.15
kspec1 = -0.15
[plt.plot(kx, strain_evals[:, i], 'k') for i in np.arange(len(centre))]
plt.plot(kx, strain_evals[:, band], 'r')

plt.plot(kspec, StrainHamiltonian(kspec, tc, tl, tr, CS[2], CLS[2], CRS[2], np.sqrt(3)*a)
                                                                        [band], 'ro')

plt.plot(kspec1, StrainHamiltonian(kspec1, tc, tl, tr, CS[2], CLS[2], CRS[2], np.sqrt(3)*a)
                                                                        [band], 'bo')
plt.show()


R_xc1 = (avy-y0)**2

def GaussianFold1(A, sigma_x, R_x):
    return A * np.exp(-R_x / (2 * sigma_x**2))


GF = GaussianFold1(A, sigma_x, R_xc1)

avy = np.convolve(ycoords[shortlist], np.ones((3,))/3, mode='valid')
avprob = np.convolve(projection[shortlist], np.ones((3,))/3, mode='valid')

fig = plt.figure()

ax1= fig.add_subplot(211)
pl = ax1.scatter(avy1, avprob1, c=avprob1, cmap="Reds", s=30, edgecolors="black")
ax2=fig.add_subplot(212)
ax2.plot(avy1, GF, c='black')

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(pl, cax=cbar_ax)

# plt.colorbar(pl)
# pl.set_label('X+Y')
plt.title("Gaussian Fold", fontsize=20, **csfont, fontweight='bold',x = -7.7, y=1.1)

plt.show()

