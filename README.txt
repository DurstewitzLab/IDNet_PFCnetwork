***********************************************************************
 Simulation code from Hass, Hertaeg and Durstewitz (2016), "A detailed,
 data data-driven network model of prefrontal cortex reproduces key 
 features of in vivo activity", PLoS Comput Biol
***********************************************************************

This package should provide all files needed to simulate the network 
model introduced in the paper above. Run the simulation with 1000 
neurons for 1000 ms using the script 'RunIDNet.m'. These parameters can
be changed in the script, as well as a number of others which are often
varied. All other parameters are defined in 'ConfigIDNet.m' and stored
in a common structure 'SimPar'. In 'IDNetSim.m', parameters specified 
by distributions are randomly drawn and the actual simulation program 
'IDNet.c' is run as a MEX file. The C code is compiled at the beginning 
of 'RunIDNet.m', this line can be commented out after the first run. The
program generates a large number of temporary files ending with '.dat',
which are deleted after the results are stored in the output file.

The simulation produces a spike train STMtx for all simulated neurons 
as well as more detailed variables such as the membrane potential V for
those neurons specified in the parameter 'ViewList'. Finally, it uses 
STMtx to draw a raster plot of the spike times. 

If model parameters are changed that affect the postsynaptic potential 
(PSP) that is elicited in any neuron type by a given synaptic peak 
conductance (gmax), the script 'Run_update_inv_con_PSP.m' needs to be 
run after making a number of changes to the code that are described in 
that script. The results are used to change the parameters of the 
function 'inv_con_PSP.m' in the way that is also described in the script. 


For questions or comments please contact:
Dr Joachim Hass (joachim.hass@zi-mannheim.de)
or 
Prof Daniel Durstewitz (daniel.durstewitz@zi-mannheim.de)

Dept. Theoretical Neuroscience, Bernstein-Center for Computational
Neuroscience, Central Institute of Mental Health, Medical Faculty
Mannheim of Heidelberg University, Mannheim, Germany


Copyright (C) 2016 for all software in this package by the Authors of
the above cited paper. This software may be used under the terms of the 
General Public License (www.gnu.org/copyleft/gpl.txt).
