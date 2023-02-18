#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#-----------------Real space Hamiltonian generator for the finite sheet.-----------------#

t0 = -2.8
epsilon = 0.02
positions = np.array((x, y)).T
dist = distance_matrix(positions, positions)
hops=np.isclose(dist, acc )
H = np.where(hops==True, t0, 0.0)
onsite_E = epsilon* np.ones((n))
H = H + np.diag(onsiteE)

#------------------------------------Energy Spectrum-------------------------------------#

#Solving for the eigenvalues, eigenvectors.
eigenvalues, eigenvectors = np.linalg.eig(H)
elements = eigenvalues.argsort()
eigenvalues = eigenvalues[elements]
eigenvectors = eigenvectors[:,elements]

fix, ax = plt.subplots(figsize=(3, 8))
[ax.plot( [0, 1], [val, val], 'k-' ) for val in eigenvalues]
plt.show()

#---------------------------------------Eigenmodes---------------------------------------#

state = 50
normal = 1000 ## factor between the probability and the dot radius in the plot

fig2, ax =plt.subplots()
ax.set_title("Eigenmode Spectrum for 5 x 5 ", fontsize=20, **csfont)
plt.xlabel('Armchair Direction', fontsize=17.5, **csfont)
plt.ylabel('Zig Zag Direction', fontsize=17.5, **csfont)
ax.scatter(positions[:,0], positions[:,1],s=normal*abs(eigenvectors[:,state])**2)
plt.plot()

#-----------------------------------------DOS--------------------------------------------#
n = 100
En = np.linspace(-3, 3, n)
b=10E-2


def Lorentzian(En,vals,b):
    a = (En - vals)/(b/2)
    return (2 / (np.pi * b)) / (1+a**2)


def DOS(En, vals, b, n):
    dos = np.zeros([n])
    for eig in vals:
        dos += Lorentzian(En, eig, b)
    return dos/len(vals)


csfont = {'fontname':'Times New Roman'}
plt.plot(En, DOS(En, vals, b, n),label="$D(E)$", c='black')
plt.xlabel("E",fontsize=18, **csfont)
plt.ylabel("D(E)",fontsize=18, **csfont)
plt.title("Density of States (DOS) - Finite Model", fontsize=20, **csfont)
plt.axvline(0, c='gray',ls='--', lw=1)
plt.axhline(0.04, c='gray',ls='--', lw=1)
plt.legend(loc='best', fontsize=10)
plt.show()

