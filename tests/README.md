#saha_hydrogen_eos_astrobear_test.f90 and saha_hydrogen_eos_astrobear_test_declarations.f90
_Overview_     
This module is used to test the temperature interpolation techniques used in saha_hydrogen_eos_astrobear.f90. It first creates the table in the same manner as saha_hydrogen_eos_astrobear.f90. Next, for each density in a range of evenly distributed densities and each internal energy in a range of even distributed internal energies, a temperature is interpolated. With that temperature, the Saha equation is used to iterate until convergence, and the internal energy is recalculated. The percent error between the initial internal energy and the interpolated internal energy is calculated and stored in the table. The maximum and average values can then be found from the outputted file.    
_To compile_: gfortran -o saha_hydrogen_eos_astrobear_test saha_hydrogen_eos_astrobear_test_declarations.f90 saha_hydrogen_eos_astrobear_test.f90 -freal-4-real-16     
_Future Goals_     
The only thing that could be done with this file more fine-grained testing, with perhaps smaller density or internal energy increments. The techniques developed in this file should be replicated to test temperature interpolation from the other quantities (pressure, specific heat, etc.). This file can also be altered to test the accuracy of different interpolation techniques.     

