#!/usr/bin/env python
# coding: utf-8

# In[ ]:


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
ax.set_title("5 x 5 Monolayer Graphene", fontsize=20, **csfont, x= 0.55, y=1.1)
plt.xlabel('Armchair Direction', fontsize=17.5, **csfont)
plt.ylabel('Zig Zag Direction', fontsize=17.5, **csfont)
ax.plot(x[np.where(sublattice==1)], y[np.where(sublattice==1)], 'o',
markersize=4,c='purple')
ax.plot(x[np.where(sublattice==-1)], y[np.where(sublattice==-1)], 'o',
markersize=4,c='black')

