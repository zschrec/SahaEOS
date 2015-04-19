! /home/schrecengost/Documents/AstroBEAR/openmpi-1.8.4/build/bin/mpif90 -o saha_hydrogen saha_hydrogen_declarations.f90 saha_hydrogen.f90


PROGRAM saha_hydrogen

    USE SahaHydrogenDeclarations
    IMPLICIT NONE
    
    PRINT *,'Enter initial electron density'
    READ *, N_e

    N_e1 = 0
    N_e_initial = N_e
    N_e_last = N_e
 
    OPEN(10, FILE='hydrogen.dat')

    WRITE(10,*) 'T',' ', 'ionization_fraction', ' ', 'Pressure', & 
                    ' ', 'Internal_Energy', ' ', 'Specific_Heat'
        
    !K
    DO T = lower_T, upper_T, increment_T
           
        N_e = N_e_initial
        N_e_last = N_e_initial
        N_e1 = 0
          
        DO 

            ratio_h(1) = const * (T)**(3.0/2.0) * & 
              (2.0 * B_h(2) * exp(-X_h(1)/(k_eV*T))) / (B_h(1) * N_e)

            y_h(1) = ratio_h(1) / (ratio_h(1) + 1.0)

            v_e(1) = 1.0*y_h(1)

            N_e1 =  ((rho*N_o) * (x(1)/A(1)) * v_e(1))

            IF (abs((N_e - N_e1)/N_e) < tol) exit

            IF(abs((N_e_last - N_e1)/N_e_last) > 1e-3) then 
              N_e_last = N_e
              N_e = N_e1
            ELSE 
              N_e_last = N_e
              N_e = (N_e + N_e1) / 2.0
            END IF
        
        END DO  
 
        WRITE(10,*) T, ' ', ratio_h(1) / (ratio_h(1) + 1.0), &
                        ' ', Pressure(ratio_h(1) / (ratio_h(1) + 1.0),T), &
                        ' ', InternalEnergy(ratio_h(1) / (ratio_h(1) + 1.0),T), &
                        ' ', SpecificHeatConstantVolume(ratio_h(1) / (ratio_h(1) + 1.0),T)

    END DO

END PROGRAM saha_hydrogen


