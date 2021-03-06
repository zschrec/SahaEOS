﻿# SahaEOS
This is an implementation of the Saha equation of state (EOS). These modules form the foundation of implementing the Saha EOS in AstroBEAR (https://astrobear.pas.rochester.edu/trac/). The paper on the algorithms and test of AstroBEAR can be found at http://iopscience.iop.org/0067-0049/182/2/519/pdf/0067-0049_182_2_519.pdf. Each module consists of the main program and a declarations file containing everything the main program needs. The three main modules are described below, as well as what must be done in the future. The initial files folder contains the “pre-cleaned” versions of the files that are included only for the sake of completeness. The outfiles folder contains the output files of the modules, and the tests folder contains the module used to test the table generation and temperature interpolation techniques. All units are in CGS unless the variable contains “_eV” at the end, where those variables are in electron-volts. It is highly recommended that Math_Capstone_Paper.pdf is read before trying to understand these modules.   

#saha_solver.f90 with saha_solver_declarations.f90
_Overview_    
This module uses the Saha equation to iterate until electron density converges. It does so over a range of temperatures with specified compositions of hydrogen, helium, carbon, nitrogen, and oxygen. It outputs a file containing the various ionization fractions, and the ionization energy, for each temperature. Due to the lack of partition function values for carbon, nitrogen, and oxygen, we only consider the zeroth, first, and second ionization levels.    
_To compile_: gfortran -o saha_solver saha_solver_declarations.f90 saha_solver.f90    
_Configuration_        
To run an iteration with different concentrations of elements, the values in the “x” array can be changed. The index-element pairs are 1-hydrogen, 2-helium, 3-carbon, 4-nitrogen, and 5-oxygen.  
The upper_T, lower_T, and increment_T values can be changed to iterate over different temperature ranges.   
The density of the mixture can be changed by changing the value of the variable rho.  
The initial electron density to start the iterations and the tolerance of the calculation can be changed, however their current values have worked well in our tests.    
_Future Goals_    
The main thing that should be added here is determining the partition function values for carbon, nitrogen, and oxygen. Once these values are determined, these elements can be fully added to the iterations. Refer to saha_math.pdf to understand how to calculate the ionization degree using only the ratios that the Saha equation gives. 

#saha_hydrogen_eos.f90 with saha_hydrogen_eos_declarations.f90
_Overview_    
While only considering hydrogen, this module uses the Saha equation to iterate until electron density converges. Once convergence is reached, the pressure, internal energy, and specific heat capacity are calculated. This is done over a range of temperatures. It outputs a file containing the ionization fraction, the pressure, the internal energy, and the specific heat capacity for each temperature.    
_To compile_: gfortran -o saha_hydrogen_eos saha_hydrogen_eos_declarations.f90 saha_hydrogen_eos.f90     
_Configuration_        
The upper_T, lower_T, and increment_T values can be changed to iterate over different temperature ranges.   
The density of the mixture can be changed by changing the variable rho.  
The initial electron density to start the iterations and the tolerance of the calculation can be changed, however their current values have worked well in our tests.    
_Future Goals_    
More elements should be included in the calculations of temperature, internal energy, and specific heat. To fully add more than just helium, the future goals of saha_solver.f90 must be completed first. The equations to include additional elements also must be developed. See eos_math.pdf for the initial formulations of how to do this.

#saha_hydrogen_eos_astrobear.f90 with saha_hydrogen_eos_astrobear_declarations.f90
_Overview_    
This is the module that generates the tables that AstroBEAR will use. Considering only hydrogen,  the current program generates a table for a range of temperatures and densities. For each temperature-density pair, the Saha equation is used to iterate until electron density converges. Once convergence is reached, internal energy is calculated and stored in the table. Once the table is created, the user can enter a density and internal energy and the temperature is interpolated from the table. The table is generated such that the values around the edges of the table are outside our “region of interest”.        
_To compile_: gfortran -o saha_hydrogen_eos_astrobear saha_hydrogen_eos_astrobear_declarations.f90 saha_hydrogen_eos_astrobear.f90 -freal-4-real-16     
_Configuration_        
The upper_T, lower_T, and increment_T values can be changed to iterate over different temperature ranges.   
The upper_rho, lower_rho, and increment_rho values can be changed to iterate over different density ranges.   
The initial electron density to start the iterations and the tolerance of the calculation can be changed, however their current values have worked well in our tests.    
_Future Goals_     
There are several tasks that must be accomplished in order to integrate this into AstroBEAR. The first is to implement the methods to interpolate temperature using the other variables (pressure, specific heat, total energy). This requires the generation of tables similar to that of internal energy. Next, these tables must be added to AstroBEAR such that they are generated at the beginning of the simulation. Due to the amount of time it takes to generate these tables, it may be beneficial to look into doing this using parallel techniques. The functions in the saha.f90 module in AstroBEAR must then be changed so that they interpolate temperature from the table, given a thermodynamic variable based on the form of the “q” array. From the temperature, the Saha equation must be used to iterate until convergence, and then the desired quantity can be calculated. As an example, suppose you have internal energy and density and you need pressure: 1) interpolate temperature from the table 2) use that temperature and density in the iterative method with the Saha equation 3) calculate pressure from the ionization fraction that follows from convergence. The final step in integrating this with AstroBEAR is to develop a test problem and run simulations with different EOS modules to test. 
