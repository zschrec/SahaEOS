! mpif90 -o saha_hydrogen_internal_energy saha_hydrogen_declarations.f90 saha_hydrogen_internal_energy.f90


PROGRAM saha_hydrogen_internal_energy

    USE SahaHydrogenDeclarations
    IMPLICIT NONE
    
    !PRINT *,'Enter initial electron density'
    !READ *, N_e

    N_e1 = 0
    N_e_initial = N_e
    N_e_last = N_e
 
    OPEN(10, FILE='hydrogen_internal_energy.dat')

    WRITE(10,*) 'T', ' ', (rhos(rho_iter)*1.0e-3, rho_iter=1,num_rho)
    

    DO rho_iter = 1, num_rho 
        
        T_iter = 1

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
        
                N_e1 =  ((rhos(rho_iter)*N_o) * (x(1)/A(1)) * v_e(1))
        
                IF (abs((N_e - N_e1)/N_e) < tol) exit
        
                IF(abs((N_e_last - N_e1)/N_e_last) > 1e-3) then 
                    N_e_last = N_e
                    N_e = N_e1
                ELSE 
                    N_e_last = N_e
                    N_e = (N_e + N_e1) / 2.0
                END IF
                
            END DO

            !PRINT *, T_iter, ' ', rho_iter

            !internal_energies(T_iter,rho_iter) = ratio_h(1) / (ratio_h(1) + 1.0)
            internal_energies(T_iter,rho_iter) = &
                InternalEnergy_RHO(ratio_h(1) / (ratio_h(1) + 1.0), T*1.0d20,rhos(rho_iter)*1.0e-3)


            !IF (T == upper_T) THEN 
            !  PRINT *, T, ' ',internal_energies(T_iter,rho_iter)
            !END IF

            T_iter = T_iter+1

        END DO
    
    END DO

!    PRINT *, T, ' ',internal_energies(200,1)
    
!    T_iter = 1
!    DO T = lower_T, upper_T, increment_T
        !IF (T == upper_T) THEN 
        !  PRINT *, T, ' ',internal_energies(T_iter,rho_iter)
        !END IF
        !WRITE(10,*) T, ' ', (internal_energies(T_iter,rho_iter), rho_iter=1,num_rho) 
        !T_iter=T_iter+1
!      WRITE(10, "(I5, $)") T
!      DO rho_iter=1,num_rho
      
        !PRINT *, T_iter, ' ', rho_iter
!        IF (T == upper_T) THEN 
!          PRINT *, T, ' ',internal_energies(T_iter,rho_iter)
!        END IF
!        WRITE(10, "(A, F100.30,$)"), ' ', internal_energies(T_iter,rho_iter)
!      END DO
!      WRITE(10,*) ' '
            
      
!      T_iter=T_iter+1
!    END DO

  
    PRINT *,'Enter density '
    READ *, rho_enter
    PRINT *,'Enter internal energy'
    READ *, internal_energy_enter

    T_real = GetTemperatureInternalEnergy(internal_energies,rho_enter*1.0e-3,internal_energy_enter)

    PRINT *, T_real, Pressure_RHO(ratio_h(1) / (ratio_h(1) + 1.0), T_real, rho_enter*1.0e-3), &
          SpecificHeatConstantVolume(ratio_h(1) / (ratio_h(1) + 1.0), T_real)

!    WRITE(10,*) T, ' ', ratio_h(1) / (ratio_h(1) + 1.0) !, &
!                   ' ', Pressure(ratio_h(1) / (ratio_h(1) + 1.0),T), &
!                   ' ', InternalEnergy(ratio_h(1) / (ratio_h(1) + 1.0),T), &
!                   ' ', SpecificHeatConstantVolume(ratio_h(1) / (ratio_h(1) + 1.0),T)

END PROGRAM saha_hydrogen_internal_energy


