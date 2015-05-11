
PROGRAM saha_hydrogen_eos_astrobear

    USE SahaHydrogenEOSAstroBEARDeclarations
    IMPLICIT NONE
    
    N_e1 = 0
    N_e_initial = N_e
    N_e_last = N_e
 
    OPEN(10, FILE='hydrogen_internal_energy.dat')

    WRITE(10, "(A, A, $)"), ' ', 'T'
    rho_iter = 1.0*lower_rho
    DO rho_index=1,num_rho
      WRITE(10, "(A, F100.30,$)"), ' ', 10.0**rho_iter
      rho_iter=rho_iter+increment_rho
    END DO
    WRITE(10,*), ' '

    rho_index = 1
    rho_iter = 1.0*lower_rho
    DO WHILE (rho_iter < (upper_rho+increment_rho))
        
        T_index = 1

        !K
        DO T = lower_T, upper_T, increment_T

            N_e = N_e_initial
            N_e_last = N_e_initial
            N_e1 = 0
            T_real = 1.0*T

            DO 
                ratio_h(1) = const * (T)**(3.0/2.0) * & 
                    (2.0 * B_h(2) * exp(-X_h_eV(1)/(k_eV*T))) / (B_h(1) * N_e)
        
                y_h(1) = ratio_h(1) / (ratio_h(1) + 1.0)
        
                v_e(1) = 1.0*y_h(1)
        
                N_e1 =  (10.0**rho_iter*N_o) * (x(1)/A(1)) * v_e(1)
        
                IF (abs((N_e - N_e1)/N_e) < tol) exit
        
                IF(abs((N_e_last - N_e1)/N_e_last) > 1e-3) then 
                    N_e_last = N_e
                    N_e = N_e1
                ELSE 
                    N_e_last = N_e
                    N_e = (N_e + N_e1) / 2.0
                END IF
                
            END DO

            
            internal_energies(T_index,rho_index) = &
                InternalEnergy_RHO(ratio_h(1) / (ratio_h(1) + 1.0), T_real,10.0**rho_iter)

            T_index = T_index+1

        END DO
            
        rho_index = rho_index+1
        rho_iter = rho_iter+increment_rho
    END DO

    T_index = 1
    
    DO T = lower_T, upper_T, increment_T
      WRITE(10, "(I6, $)"), T
      DO rho_index=1,num_rho,1
        WRITE(10, "(A, F100.30,$)"), ' ', internal_energies(T_index,rho_index)
      END DO
      WRITE(10,*), ' '
      FLUSH(10)
      T_index=T_index+1
    END DO

    PRINT *,'Enter density '
    READ *, rho_enter
    PRINT *,'Enter internal energy'
    READ *, internal_energy_enter
  
    DO WHILE (rho_enter > 0)
  
      T_real = GetTemperatureInternalEnergy(internal_energies,rho_enter,internal_energy_enter)
  
      PRINT *, T_real, Pressure_RHO(ratio_h(1) / (ratio_h(1) + 1.0), T_real, rho_enter), &
            SpecificHeatConstantVolume(ratio_h(1) / (ratio_h(1) + 1.0), T_real)

      N_e = N_e_initial
      N_e_last = N_e_initial
      N_e1 = 0
          
      DO 
        ratio_h(1) = const * (T_real)**(3.0/2.0) * & 
          (2.0 * B_h(2) * exp(-X_h_eV(1)/(k_eV*T_real))) / (B_h(1) * N_e)
        
        y_h(1) = ratio_h(1) / (ratio_h(1) + 1.0)
        
        v_e(1) = 1.0*y_h(1)
        
        N_e1 =  (rho_enter*N_o) * (x(1)/A(1)) * v_e(1)
        
        IF (abs((N_e - N_e1)/N_e) < tol) exit
        
          IF(abs((N_e_last - N_e1)/N_e_last) > 1e-3) then 
            N_e_last = N_e
            N_e = N_e1
          ELSE 
            N_e_last = N_e
            N_e = (N_e + N_e1) / 2.0
          END IF
                
      END DO

      PRINT *, InternalEnergy_RHO(ratio_h(1) / (ratio_h(1) + 1.0), T_real,rho_enter), ' ', &
          ABS(InternalEnergy_RHO(ratio_h(1) / (ratio_h(1) + 1.0), T_real,rho_enter) - &
              internal_energy_enter)/(internal_energy_enter*1.0) * 100

    END DO

    CLOSE(10)

END PROGRAM saha_hydrogen_eos_astrobear


