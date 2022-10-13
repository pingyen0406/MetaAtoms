# Generating radius matrix of desired metalens phase profile
Choose desired phase profile and find corresponding radius of the meta-atom.
## Prerequisites 
The input data should have radii-transmission and radii-phase relations. You can use RCWA, FDTD,...etc to calculate the relation.

## Standard using flow
1. Setting parameters in the input.json(transmission & phase data, lens
size, radius & height range and step)
2. For the 1st run, set a breakpoint at line 60. You need to customly determine the
range of the radii.
3. After determining the radii range, just follow the pop out dialog.
Note that if you choose "custom", the input desired phase matrix should
have values between 0 to 2pi.

## Relevent Publications
[Integrated Metasurfaces on Silicon Photonics for Emission Shaping and Holographic Projection](https://arxiv.org/abs/2205.10537)  
Ping-Yen Hsieh, Shun-Lin Fang, Yu-Siang Lin, Wen-Hsien Huang, Jia-Min Shieh, Peichen Yu, and You-Chia Chang  
Nanophotonics (Accepted, 2022)  

[Shaping Free-space Emission with Monolithically Integrated Metalenses on Silicon Photonic Waveguides](https://doi.org/10.1364/CLEO_SI.2022.SM4P.7)  
Ping-Yen Hsieh, Yu-Siang Lin, Shun-Lin Fang, and You-Chia Chang  
CLEO 2022: Science and Innovations, SM4P.7  
