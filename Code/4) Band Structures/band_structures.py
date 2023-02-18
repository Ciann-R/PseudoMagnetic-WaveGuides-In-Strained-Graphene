#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#--------------------Band Structure in the First Brilloiun Zone--------------------------#
epsilon = 0
t = -1
a = 1
S = 50

#reciprocal space
b1 = np.array([np.pi / a*np.sqrt(3), 0])
b2 = np.array([0, np.pi / a])
bx = b1[0]
by = b2[1]

#Nearest neighbour vectors
d1 = a * np.array((-1, 0))
d2 = a/2 * np.array((1, np.sqrt(3)))
d3 = a/2 * np.array((1, -np.sqrt(3)))


def Hopping(k, r):
    Hops = t * np.exp(1j * np.dot(k, r))
    return Hops


def HoppingConj(k, r):
    Hops_Conj = t * np.exp(-1j * np.dot(k, r))
    return Hops_Conj


def Hamiltonian(k):
    H = np.matrix([[epsilon, Hopping(k, d1) + Hopping(k, d2) + Hopping(k, d3)],
                   [HoppingConj(k, d1) + HoppingConj(k, d2) + HoppingConj(k, d3),
                    epsilon]])
    Eig = np.real(np.linalg.eigvals(H))
    return Eig

k = [(kx, kx/np.sqrt(3)) for kx in np.linspace(0, bx, 50)]



csfont = {'fontname':'Times New Roman'}
fig = plt.figure(figsize=(10,10))
axis = fig.add_subplot(111, projection='3d')
axis.set_title("3D Bandstructure of Graphene", fontsize=26, **csfont, fontweight='bold')
plt.xlabel('kx', fontsize=20, **csfont,fontweight='bold')
plt.ylabel('ky', fontsize=20, **csfont, fontweight='bold')
axis.set_zlabel('E(k)',fontsize=20, **csfont, fontweight='bold')


kx = np.linspace(-2*bx, 2*bx, S)
ky = np.linspace(-2*bx/np.sqrt(3), 2*bx/np.sqrt(3), S)
X, Y = np.meshgrid(kx, ky)


Z = np.array(list(map(Hamiltonian, zip(X.flatten(), Y.flatten()))))
axis.plot_surface(X, Y, Z[:,0].reshape((S,S)), cmap=cm.Greys)
axis.plot_surface(X, Y, Z[:,1].reshape((S,S)), cmap=cm.Purples);
axis.xaxis.pane.fill = False
axis.yaxis.pane.fill = False
axis.zaxis.pane.fill = False
axis.xaxis.pane.set_edgecolor('w')
axis.yaxis.pane.set_edgecolor('w')
axis.zaxis.pane.set_edgecolor('w')
axis.grid(False)

