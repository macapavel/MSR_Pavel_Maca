# In[73]:


#        

from os import name
from textwrap import fill
#from turtle import width
import numpy as np
import matplotlib.pyplot as plt 
import openmc
import os

openmc.config['cross_sections']='/home/pmaca/git/openmc/libraries/endfb-vii.1-hdf5/cross_sections.xml'

## Remove old files
#os.system('rm geometry.xml')
#os.system('rm materials.xml')
#os.system('rm settings.xml')
#os.system('rm tallies.xml')
#os.system('rm *.h5')

## import opnmc plt
#import openmc.plot

# OpenMC simulations paramaters
#batches = 20
#inactive = 10
#particles = 10000

## Materials definitions 

# salt (FLiBeUF)

FLiBeUF = openmc.Material(1, name="FLiBeUF")
FLiBeUF.add_element('F', 2.0)
FLiBeUF.add_nuclide('Li6', 0.01)
FLiBeUF.add_nuclide('Li7', 99.99)
FLiBeUF.add_element('Be', 1.0)
FLiBeUF.add_nuclide('U235', 9.75)  
FLiBeUF.add_nuclide('U238', 79.95)
FLiBeUF.add_nuclide('U234', 0.20)
FLiBeUF.add_nuclide('U236', 0.1)
FLiBeUF.set_density('g/cm3', 2.645)


# graphite
graphite = openmc.Material(2, name='graphite')
graphite.add_element('C', 99.99980, percent_type='wo')
graphite.add_nuclide('B10', 0.00004, percent_type='wo')
graphite.add_nuclide('B11', 0.00016, percent_type='wo')
graphite.set_density('g/cm3', 1.74)
#graphite.add_s_alpha_beta('c_graphite')


# stainless steel 
steel = openmc.Material(3, name='steel')
steel.add_element('C', 0.1900, percent_type='wo')
steel.add_element('Si', 1.0048, percent_type='wo')
steel.add_element('Mn', 1.0274, percent_type='wo')
steel.add_element('P', 0.0413, percent_type='wo')
steel.add_element('S', 0.0260, percent_type='wo')
steel.add_element('Cr', 18.1986, percent_type='wo')
steel.add_element('Ni', 11.3803, percent_type='wo')
steel.add_element('Fe', 66.6811, percent_type='wo')
steel.add_element('Mo', 1.4504, percent_type='wo')
steel.set_density('g/cm3', 8.00)

# helium     
helium=openmc.Material(name='helium')
helium.add_element('He', 1.0)
helium.set_density('g/cm3', 0.000055)

# instantiate a materials collection and export to xml 
materials_file = openmc.Materials([FLiBeUF, graphite, steel, helium])
materials_file.export_to_xml()

# hexagonal lattice  


# In[3]:


r_pin = openmc.ZCylinder(r=1.5, name='salt')
FLiBeUF_cell = openmc.Cell(fill=FLiBeUF, region=-r_pin)
graphite_cell = openmc.Cell(fill=graphite, region=+r_pin)
he_cell = openmc.Cell(fill=helium, region=-r_pin)
# 
pin_universe = openmc.Universe(cells=(FLiBeUF_cell, graphite_cell))
he_universe = openmc.Universe(cells=(he_cell, graphite_cell))


all_graphite_cell = openmc.Cell(fill=graphite,)
outer_universe = openmc.Universe(cells=(all_graphite_cell,))


# In[4]:

'''
################################################
Filling a hexagon lattice
'''
hex_lattice = openmc.HexLattice()
hex_lattice.center = (-5.0, 5.0)
hex_lattice.pitch = (10.0,)  
hex_lattice.outer = outer_universe


#

# print(hex_lattice.show_indices(num_rings=10))


# In[6]:


pin_pos = np.array([
                                       [4,28],      [5, 0],      [4, 2],             
                                [4,27],      [5,23],      [5, 1],      [4, 3],             
                                       [5,22],      [6, 0],      [5, 2],      [4, 4],             
                                [5,21],      [6,17],      [6, 1],      [5, 3],              
                          [5,20],      [6,16],      [7, 0],      [6, 2],      [5, 4],             
                                [6,15],      [7,11],      [7, 1],      [6, 3],      [4, 6],      
                          [5,19],      [7,10],      [8, 0],      [7, 2],      [5, 5],            
                   [4,23],      [6,14],      [8, 5],      [8, 1],      [6, 4],      [4, 7],       
                          [5,18],      [7, 9],      [9, 0],      [7, 3],      [5, 6],      [3, 9],      
                   [4,22],      [6,13],      [8, 4],      [8, 2],      [6, 5],      [4, 8],       
                          [5,17],      [7, 8],      [8, 3],      [7, 4],      [5, 7],      [3,10],       
                   [4,21],      [6,12],      [7, 7],      [7, 5],      [6, 6],      [4, 9],       
                          [5,16],      [6,11],      [7, 6],      [6, 7],      [5, 8],              
                                [5,15],      [6,10],      [6, 8],      [5, 9],      [4,10],       
                          [4,19],      [5,14],      [6, 9],      [5,10],      [4,11],             
                                [4,18],      [5,13],      [5,11],      [4,12],             
                                       [4,17],      [5,12],      [4,13],      [3,14],       
                                [3,21],      [4,16],      [4,14],      [3,15],       
                                      [3,20],      [4,15],      [3,16],      
]) 
# space between ij is probably not necessary 

# In[7]:


########## INITIATE A LATTICE FILL WITH GRAPHITE ############# 
n = 10 # number of ring
lattice_list = []
for i in (range(n)):
   # this loop create list for each ring that filled with 'outer_universe'
   dummy_list = [outer_universe]*6*i
   lattice_list.insert(0,dummy_list) 
lattice_list[n-1] = [outer_universe]


# In[8]:

################# FILL FUEL TO CORRESPONDING POS ###############
for pos in pin_pos:
   i = int(pos[0])  
   j = int(pos[1])  
   lattice_list[i][j] = pin_universe  


# In[9]:


################# FILL He TO CORRESPONDING POS ###############

he_pos = np.array([
                                                      [0, 0],
                                                [0,53],      [0, 1],
                                          [0,52],      [1, 0],      [0, 2],
                                    [0,51],      [1,47],      [1, 1],      [0, 3],
                              [0,50],      [1,46],      [2, 0],      [1, 2],      [0, 4],
                        [0,49],      [1,45],      [2,41],      [2, 1],      [1, 3],      [0, 5],
                  [0,48],      [1,44],      [2,40],      [3, 0],      [2, 2],      [1, 4],      [0, 6],
            [0,47],      [1,43],      [2,39],      [3,35],      [3, 1],      [2, 3],      [1, 5],      [0, 7],
      [0,46],      [1,42],      [2,38],      [3,34],      [4, 0],      [3, 2],      [2, 4],      [1, 6],      [0, 8],
[0,45],      [1,41],      [2,37],      [3,33],      [4,29],      [4, 1],      [3, 3],      [2, 5],      [1, 7],      [0, 9],
      [1,40],      [2,36],      [3,32],                                             [3, 4],      [2, 6],      [1, 8],
[0,44],      [2,35],      [3,31],                                                          [3, 5],      [2, 7],      [0,10],
      [1,39],      [3,30],      [4,26],                                                          [3, 6],      [1, 9],
[0,43],      [2,34],      [4,25],                                                          [4, 5],      [2, 8],      [0,11],
      [1,38],      [3,29],                                                                          [3, 7],      [1,10],
[0,42],      [2,33],      [4,24],                                                                       [2, 9],      [0,12],
      [1,37],      [3,28],                                                                          [3, 8],      [1,11],
[0,41],      [2,32],                                                                                    [2,10],      [0,13],
      [1,36],      [3,27],                                                                                       [1,12],
[0,40],      [2,31],                                                                                    [2,11],      [0,14],
      [1,35],      [3,26],                                                                                       [1,13],
[0,39],      [2,30],                                                                                    [2,12],      [0,15],
      [1,34],      [3,25],                                                                          [3,11],      [1,14],
[0,38],      [2,29],      [4,20],                                                                       [2,13],      [0,16],
      [1,33],      [3,24],                                                                          [3,12],      [1,15],
[0,37],      [2,28],      [3,23],                                                          [3,13],      [2,14],      [0,17],
      [1,32],      [2,27],      [3,22],                                                             [2,15],      [1,16],
[0,36],      [1,31],      [2,26],                                                          [2,16],      [1,17],      [0,18],
      [0,35],      [1,30],      [2,25],                                           [2,17],      [1,18],      [0,19],
            [0,34],      [1,29],      [2,24],      [3,19],      [3,17],      [2,18],      [1,19],      [0,20],
                  [0,33],      [1,28],      [2,23],      [3,18],      [2,19],      [1,20],      [0,21],
                        [0,32],      [1,27],      [2,22],      [2,20],      [1,21],      [0,22],
                              [0,31],      [1,26],      [2,21],      [1,22],      [0,23],
                                    [0,30],      [1,25],      [1,23],      [0,24],
                                          [0,29],      [1,24],      [0,25],
                                                [0,28],      [0,26],
                                                      [0,27],                                                           
])
#for pos in he_pos:
#    i = int(pos[0]) # First indice is i
#    j = int(pos[1]) # Second indice is j
#    lattice_list[i][j] = he_universe # assign material (fuel pin) to position


# In[10]:
#print(lattice_list)

hex_lattice.universes = lattice_list[:]

# geometry plotting

# def geom; export to xml

# tally 
#tally = openmc.Tally(3)
#tally.socres = ['flux']
#tally. filters = [openmc.CellFilter(FLiBeUF_cell)]
#tallies = openmc.Tallies([tally])


# In[11]:


# print(hex_lattice)


# In[74]:

#        bottom is on -70, then add 140, so top is 70  radius 66 
c_cont = openmc.model.RightCircularCylinder([0,0,-70],140,66)

inner_region = -c_cont
inner_cell = openmc.Cell(region = inner_region)
inner_cell.fill = hex_lattice

# In[75]:


flibe_cont = openmc.model.RightCircularCylinder([0,0,-72],144,68)

flibe_region = -flibe_cont & +c_cont
flibe_cell = openmc.Cell(region = flibe_region)
flibe_cell.fill = FLiBeUF

# In[76]:

# steel radius just (should be) 71 
steel_cont = openmc.model.RightCircularCylinder([0,0,-75],150,71)
steel_tube = openmc.model.RightCircularCylinder([0,0,-150],300,6)

flibe_tube = openmc.model.RightCircularCylinder([0,0,-150],300,5)

he_cont = openmc.model.RightCircularCylinder([0,0,-150],300,90)

steel_region = (+flibe_cont & -steel_cont & +steel_tube) | (-steel_tube & +flibe_tube & + flibe_cont & -he_cont)
steel_cell = openmc.Cell(region = steel_region)
steel_cell.fill = steel
######  
flibe_tube_region = -flibe_tube & +flibe_cont 
flibe_tube_cell = openmc.Cell(region=flibe_tube_region)
flibe_tube_cell.fill = FLiBeUF

he_region = -he_cont & +steel_cont & +steel_tube
he_cell =  openmc.Cell(region=he_region)
he_cell.fill = helium


# In[77]:

# shoud be the radius 100 not 110???  
outer_surface = openmc.ZCylinder(r=100, boundary_type='vacuum')
outer_cell = openmc.Cell(region= +he_cont & -outer_surface )
# add the geometry 
universe = openmc.Universe(cells=[inner_cell,flibe_cell,flibe_tube_cell,steel_cell,he_cell,outer_cell])
geometry = openmc.Geometry(universe)
geometry.export_to_xml()


# In[78]:

#plt 
plot = openmc.Plot.from_geometry(geometry)
plot.color_by = 'material'
plot.basis = 'xy'
plot.colors = colors = {
    graphite: 'aquamarine',
    FLiBeUF: 'blue',
    steel: 'green',
    helium: 'yellow'
}

# display the picture
plot.to_ipython_image()
plt.show()
#plot.image_format = 'png'

#plt xz
#plt.show(universe.plot(basis='yz',origin = (0,0,0)))

#plt.show(universe.plot(width=(120, 320), origin = [0,0,0], basis='xz'))



### tally like ,, Tallies vsechny tally,tally- prvni tally, k_effective_tally - druhe tally
tallies = openmc.Tallies()
# crate mesh which will be used for tally 
mesh = openmc.RegularMesh()
mesh.dimension = [100, 100]
mesh.lower_left = [-110, -110]
mesh.upper_right = [110, 110]

# create mesh filter for tally 
mesh_filter = openmc.MeshFilter(mesh)

# create mesh tally to score flux and fission rate 
tally = openmc.Tally(name='flux')
tally.filters = [mesh_filter]
tally.scores = ['flux', 'fission']
tallies.append(tally)



###
k_effective_tally = openmc.Tally(name='k-effective')
k_effective_tally.scores = ['k-effective']
# add tally to model 
#tallies.append(k_effective_tally)  
tallies.export_to_xml()


##### setting like in BWR example - in 27 

# openmc simulation parameters
point = openmc.stats.Point((0, 0, 0))
src = openmc.Source(space=point)
settings = openmc.Settings()
settings.source = src
settings.batches = 100
settings.inactive = 10
settings.particles = 10000
settings.output = {'tallies': True}
settings.export_to_xml()  


# connecting all parts of model
# the order the setting and tally must be first setting then tally! 
model = openmc.model.Model(geometry, materials_file, settings,tallies)

# run simulation
output = openmc.run()
#openmc.run()

## obtaining results 
##k_effective_results = output.tallies['k-effective'].mean 
#print(f'Kriticky faktor (k-effective): {k_effective_tally}')



################################################




