# SahaEOS
This is an implementation of the Saha equation of state (EOS). These modules form the foundation of implementing the Saha EOS in AstroBEAR (https://astrobear.pas.rochester.edu/trac/). Each module consists of the main program and a declarations files containing everything the main program needs. The three main modules are be described below, as well as what must be done in the future. The initial files folder contains the “pre-cleaned” versions of the files that are included only for the sake of completeness. All units are in CGS unless the variable contains “_eV” at the end, where those variables are in electron-volts. It is highly recommended that the Math_Capstone_Paper.pdf is read before trying to understand these modules.   

#saha_solver.f90 with saha_solver_declarations.f90
This module uses the Saha equation to iterate until electron density converges. It does so over a range of temperatures with specified compositions of hydrogen, helium, carbon, nitrogen, and oxygen. It outputs a file containing the various ionization fractions and the ionization energy for each temperature. Due to the lack of partition function values for carbon, nitrogen, and oxygen, we only consider the zeroth, first, and second ionization levels.    
To run the iteration with different concentrations of elements, the values in the “x” array can changed. The index-element pairs are 1-hydrogen, 2-helium, 3-carbon, 4-nitrogen, and 5-oxygen.  
The upper_T, lower_T, and increment_T values can be changed to iterate over different temperature ranges.   
The density of the mixture can be changed by changing the value of the variable rho.  
The initial electron density to start the iteration and the tolerance of the calculation can be changed, however their current values have worked well in our tests.    
_Future Goals_    
The main thing that should be added here is determining the partition function values for carbon, nitrogen, and oxygen. Once these values are determined, they can be fully added to the iterations. Refer to the saha_math.doc file to understand how to calculate the ionization degree using only the ratios that the Saha equation gives. 

#saha_hydrogen_eos.f90 with saha_hydrogen_eos_declarations.f90
While only considering hydrogen, this module uses the Saha equation to iterate until electron density converges. Once convergence is reached, the pressure, internal energy, and specific heat capacity are calculated and it does this over a range of temperatures. It outputs a file containing the ionization fraction, the pressure, the internal energy, and the specific heat capacity for each temperature.  
The upper_T, lower_T, and increment_T values can be changed to iterate over different temperature ranges.   
The density of the mixture can be changed by changing the variable rho.  
The initial electron density to start the iteration and the tolerance of the calculation can be changed, however their current values have worked well in our tests.    
_Future Goals_    
More elements should be included in the calculations of temperature, internal energy, and specific heat. To fully add more than just helium, the future goals of saha_solver.f90 must be completed first. The equations to include additional elements also must be developed. 


#saha_hydrogen_eos_astrobear.f90 with saha_hydrogen_eos_astrobear_declarations.f90
This is the module that generates the tables that AstroBEAR will use. Considering only hydrogen,  the current program generates a table for a range of temperatures and densities. For each temperature-density pair, the Saha equation is used to iterate until electron density converges. Once convergence is reached, internal energy is calculated and stored in the table. Once the table is created, the user can enter a density and internal energy and the temperature is interpolated from the table. The table is generated such that the values around the edges of the table are outside our “region of interest”. 
