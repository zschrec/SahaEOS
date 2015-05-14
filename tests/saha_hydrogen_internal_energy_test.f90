! mpif90 -o saha_hydrogen_internal_energy saha_hydrogen_declarations.f90 saha_hydrogen_internal_energy.f90


PROGRAM saha_hydrogen_internal_energy

    USE SahaHydrogenDeclarations
    IMPLICIT NONE
   
    REAL (KIND=ikind), DIMENSION(500, 39) :: ie_test
    REAL (KIND=ikind), DIMENSION(500) :: T_test
    REAL (KIND=ikind) :: ie_iter
    REAL (KIND=ikind) :: ie_inc
    INTEGER :: ie_iter_index
    INTEGER :: ie_iter_l_index = 9
    INTEGER :: ie_iter_u_index = 193
    INTEGER :: ie_iter_num = 500
    REAL :: rho_iter_test
    INTEGER :: rho_index_test

!    PRINT *,(upper_rho-lower_rho)/increment_rho
!    rho_iter = 10.0**lower_rho + (10.0**upper_rho-10.0**lower_rho)/((upper_rho-lower_rho)/increment_rho)
!    rho_index = 2
!    DO WHILE (rho_iter < 10.0**upper_rho) 
!      PRINT *, rho_iter, rho_index
!      rho_iter = rho_iter + (10.0**upper_rho-10.0**lower_rho)/((upper_rho-lower_rho)/increment_rho)
!      rho_index = rho_index + 1
!    
!    END DO
!
!    RETURN



    N_e1 = 0
    N_e_initial = N_e
    N_e_last = N_e
 
    OPEN(10, FILE='hydrogen_internal_energy_test.dat')

    rho_iter_test = 10.0**lower_rho + (10.0**upper_rho-10.0**lower_rho)/((upper_rho-lower_rho)/increment_rho)
    DO WHILE (rho_iter_test < 10.0**(upper_rho-increment_rho)) 
      WRITE(10, "(A, F100.30,$)"), ' ', rho_iter_test
      rho_iter_test = rho_iter_test + (10.0**upper_rho-10.0**lower_rho)/((upper_rho-lower_rho)/increment_rho)
    END DO
    WRITE(10,*), ' '

    rho_index = 1
    rho_iter = 1.0*lower_rho
    DO WHILE (rho_iter < (upper_rho+increment_rho))
        
        T_iter = 1

        !K
        DO T = lower_T, upper_T, increment_T

            N_e = N_e_initial
            N_e_last = N_e_initial
            N_e1 = 0
          
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

            internal_energies(T_iter,rho_index) = &
                InternalEnergy_RHO(ratio_h(1) / (ratio_h(1) + 1.0), T*1.0,10.0**rho_iter)

            T_iter = T_iter+1

        END DO
            
        rho_index = rho_index+1
        rho_iter = rho_iter+increment_rho
    END DO

    PRINT *, 'Table Done'


 !   rho_index = 1
 !   rho_iter = 1.0*lower_rho
 !   DO WHILE (rho_iter < (upper_rho+increment_rho))
        

 !       rho_index = rho_index+1
 !       rho_iter = rho_iter+increment_rho
 !   END DO

    PRINT *,(upper_rho-lower_rho)/increment_rho
    rho_iter_test = 10.0**lower_rho + (10.0**upper_rho-10.0**lower_rho)/((upper_rho-lower_rho)/increment_rho)
    rho_index_test = 2
    DO WHILE (rho_iter_test < 10.0**(upper_rho-increment_rho)) 
      PRINT *, rho_iter_test, rho_index_test
      PRINT *, internal_energies(ie_iter_l_index, rho_index_test), internal_energies(ie_iter_u_index, rho_index_test)

      ie_inc = (internal_energies(ie_iter_u_index, rho_index_test) - internal_energies(ie_iter_l_index, rho_index_test)) /&
          (ie_iter_num*1.0)
      ie_iter = internal_energies(ie_iter_l_index, rho_index_test)
      
      N_e1 = 0
      N_e_initial = N_e
      N_e_last = N_e
      !DO ie_iter_index = ie_iter_l_index, ie_iter_u_index
      DO ie_iter_index = 1, 500
 
        !PRINT *, ie_iter_index

        T_real = GetTemperatureInternalEnergy(internal_energies,rho_iter_test,ie_iter)
        
        !PRINT *, T_real

        N_e = N_e_initial
        !PRINT *, N_e
        N_e_last = N_e_initial
        !PRINT *, N_e_last
        N_e1 = 0
        !PRINT *, N_e1
          
        DO 

          
          ratio_h(1) = const * (T_real)**(3.0/2.0) * & 
            (2.0 * B_h(2) * exp(-X_h_eV(1)/(k_eV*T_real))) / (B_h(1) * N_e)
        
          y_h(1) = ratio_h(1) / (ratio_h(1) + 1.0)
        
          v_e(1) = 1.0*y_h(1)
        
          N_e1 =  (rho_iter_test*N_o) * (x(1)/A(1)) * v_e(1)
        
          IF (abs((N_e - N_e1)/N_e) < tol) exit
        
          IF(abs((N_e_last - N_e1)/N_e_last) > 1e-3) then 
            N_e_last = N_e
            N_e = N_e1
          ELSE 
            N_e_last = N_e
            N_e = (N_e + N_e1) / 2.0
          END IF

          !PRINT *, N_e


                
        END DO

        !ie_test(ie_iter_index, rho_index_test-1) = ie_iter

        ie_test(ie_iter_index, rho_index_test-1) = &
          ABS(InternalEnergy_RHO(ratio_h(1) / (ratio_h(1) + 1.0), T_real,rho_iter_test)-ie_iter)/ie_iter * 100.0

        ie_iter = ie_iter + ie_inc
      END DO

      rho_iter_test = rho_iter_test + (10.0**upper_rho-10.0**lower_rho)/((upper_rho-lower_rho)/increment_rho)
      rho_index_test = rho_index_test + 1
   
    END DO

    DO T = 1, 500, 1

      DO rho_index=1,39,1
        WRITE(10, "(A, F100.30,$)"), ' ', ie_test(T,rho_index)
      
      END DO
      
      WRITE(10,*), ' '
      
      FLUSH(10)
      
    END DO

    CLOSE(10)

END PROGRAM saha_hydrogen_internal_energy


