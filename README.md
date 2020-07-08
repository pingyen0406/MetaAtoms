# Simulating metalens as the superposition of the point source
This method takes every meta-atom as an independent point source(electric field of spherical wave)  and adds them together (superposition of scalar field).
## Input data 
The input data should have radii-transmission and radii-phase relations. 

## Processing flow
1. Read the input json file
2. Creat a 2xN matrix which contains the position of every meta-atom
3. Do massage to the input data and create the library radius, transmission and phase
4. Generate the designed phase profile and use interpolation to put corresponding meta-atom to each position
5. Calculate the field by the transmission and initial phase data and plot it out
One_dim.m and Two_dim.m are 1-D and 2-D examples
