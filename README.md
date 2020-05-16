# Simulating metalens as the superposition of the point source
This method takes every meta-atom as an independent point source and superpose them together. 
## Input data 
The input data should have radius-transmission and radius-phase relations. 
## General usage
1. Read the input file
2. Creat a 2xN matrix which contains the position of every meta-atom
3. Do massage to the input data and create the library
4. Generate the designed phase profile and use interpolation to put corresponding meta-atom to each position
5. Calculate the field by the transmission and initial phase data and plot it out

  One_dim.m and Two_dim.m correspond to 1-D and 2-D examples
