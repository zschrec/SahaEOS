

PROGRAM saha_solver
  USE SahaSolverDeclarations
  IMPLICIT NONE
  
  N_e_initial = N_e
  N_e_last = N_e
 
  OPEN(10, FILE='outfiles/saha_solver.dat')
  WRITE(10,*) 'T', ' ', 'n_h1/n_h', ' ', 'n_he1/n_he', ' ', 'n_he2/n_he', ' ', &
        'n_c1/n_c', ' ', 'n_c2/n_c', ' ', 'n_n1/n_n', ' ', 'n_n2/n_n', ' ', &
        'n_o1/n_o', ' ', 'n_o2/n_o', ' ', 'ionization_energy' 
       
  !K
  DO T = lower_T, upper_T, increment_T
    
    N_e = N_e_initial
    N_e_last = N_e_initial
    N_e1 = 0
          
    DO

      ! Calculate the Saha Equation for each element
      ratio_h(1) = const * (T)**(3.0/2.0) * & 
        (2.0 * B_h(2) * EXP(-X_h_eV(1)/(k_eV*T))) / (B_h(1) * N_e)

      ratio_he(1) = const * (T)**(3.0/2.0) * & 
        (2.0 * B_he(2) * EXP(-X_he_eV(1)/(k_eV*T))) / (B_he(1) * N_e)
      ratio_he(2) = const * (T)**(3.0/2.0) * & 
        (2.0 * B_he(3) * EXP(-X_he_eV(2)/(k_eV*T))) / (B_he(2) * N_e)
      
      ratio_c(1) = const * (T)**(3.0/2.0) * & 
        (2.0 * B_c(2) * EXP(-X_c_eV(1)/(k_eV*T))) / (B_c(1) * N_e)
      ratio_c(2) = const * (T)**(3.0/2.0) * & 
        (2.0 * B_c(3) * EXP(-X_c_eV(2)/(k_eV*T))) / (B_c(2) * N_e)
      
      ratio_n(1) = const * (T)**(3.0/2.0) * & 
        (2.0 * B_n(2) * EXP(-X_n_eV(1)/(k_eV*T))) / (B_n(1) * N_e)
      ratio_n(2) = const * (T)**(3.0/2.0) * & 
        (2.0 * B_n(3) * EXP(-X_n_eV(2)/(k_eV*T))) / (B_n(2) * N_e)
      
      ratio_o(1) = const * (T)**(3.0/2.0) * & 
        (2.0 * B_o(2) * EXP(-X_o_eV(1)/(k_eV*T))) / (B_o(1) * N_e)
      ratio_o(2) = const * (T)**(3.0/2.0) * & 
        (2.0 * B_o(3) * EXP(-X_o_eV(2)/(k_eV*T))) / (B_o(2) * N_e)
      
      ! Calculate the degree of ionization (aka ionization fraction)
      y_h(1) = 1.0 / ((1.0/ratio_h(1)) + 1.0)

      y_he(1) = 1.0 / ((1.0/ratio_he(1)) + 1 + ratio_he(2))
      y_he(2) = 1.0 / ( (1.0 / ( ratio_he(1)*ratio_he(2) ))  + (1.0 / ratio_he(2) ) + 1.0)
      
      y_c(1) = 1.0 / ((1.0/ratio_c(1)) + 1 + ratio_c(2))
      y_c(2) = 1.0 / ( (1.0 / ( ratio_c(1)*ratio_c(2) ))  + (1.0 / ratio_c(2) ) + 1.0)
      
      y_n(1) = 1.0 / ((1.0/ratio_n(1)) + 1 + ratio_n(2))
      y_n(2) = 1.0 / ( (1.0 / ( ratio_n(1)*ratio_n(2) ))  + (1.0 / ratio_n(2) ) + 1.0)
      
      y_o(1) = 1.0 / ((1.0/ratio_o(1)) + 1 + ratio_o(2))
      y_o(2) = 1.0 / ( (1.0 / ( ratio_o(1)*ratio_o(2) ))  + (1.0 / ratio_o(2) ) + 1.0)
      
      ! Calculate average electron contribution due to each element
      v_e(1) = 1.0*y_h(1)
      v_e(2) = 1.0*y_he(1) + 2.0*y_he(2)
      v_e(3) = 1.0*y_c(1) + 2.0*y_c(2)
      v_e(4) = 1.0*y_n(1) + 2.0*y_n(2)
      v_e(5) = 1.0*y_o(1) + 2.0*y_o(2)
            
      ! Recompute electron density
      N_e1 =  (rho*N_o) * &
        ((x(1)/A(1)) * v_e(1) + & 
        (x(2)/A(2)) * v_e(2) + &
        (x(3)/A(3)) * v_e(3) + &
        (x(4)/A(4)) * v_e(4) + &
        (x(5)/A(5)) * v_e(5))
     
      IF (ABS((N_e - N_e1)/N_e) < tol) EXIT
        
      IF(ABS((N_e_last - N_e1)/N_e_last) > 1e-3) THEN 
        N_e_last = N_e
        N_e = N_e1
      ELSE 
        N_e_last = N_e
        N_e = (N_e + N_e1) / 2.0
      END IF
        
    END DO  

    ! Log temperature, ionization fractions, and ionization energy
    WRITE(10,*) T,' ', y_h(1), ' ', y_he(1), ' ', y_he(2), ' ', y_c(1), ' ', y_c(2), ' ', &
      y_n(1), ' ', y_n(2), ' ', y_o(1), ' ', y_o(2), ' ', &
      (y_h(1))* (rho*N_o*x(1)/A(1))*X_h_eV(1) + &
      (y_he(1))*(rho*N_o*x(2)/A(2))*X_he_eV(1) + &
      (y_he(2))*(rho*N_o*x(2)/A(2))*X_he_eV(2) + &
      (y_c(1))* (rho*N_o*x(3)/A(3))*X_c_eV(1) + &
      (y_c(2))* (rho*N_o*x(3)/A(3))*X_c_eV(2) + &
      (y_n(1))* (rho*N_o*x(4)/A(4))*X_n_eV(1) + &
      (y_n(2))* (rho*N_o*x(4)/A(4))*X_n_eV(2) + &
      (y_o(1))* (rho*N_o*x(5)/A(5))*X_o_eV(1) + &
      (y_o(2))* (rho*N_o*x(5)/A(5))*X_o_eV(2)

          
  END DO

END PROGRAM saha_solver

