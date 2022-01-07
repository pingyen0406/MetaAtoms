# Generating radius matrix of desired metalens phase profile
Choose desired phase profile and find corresponding radius of the meta-atom.
## Prerequisites 
The input data should have radii-transmission and radii-phase relations. You can use RCWA, FDTD,...etc to calculate the relation.

## Processing flow
1. Setting parameters in the input.json(transmission & phase data, lens
size, radius & height range and step)
2. For the 1st run, set a breakpoint at line 60. You need to choose the
range of the radii.
3. After determining the radii range, just follow the pop out dialog.
Note that if you choose "custom", the input desired phase matrix should
have value between 0 to 2pi.

