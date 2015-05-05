!


PROGRAM saha_helium
  USE SahaSolverDeclarations
  IMPLICIT NONE
  
  N_e_initial = N_e
  N_e_last = N_e
 
  OPEN(10, FILE='outfiles/saha_helium.dat')
  WRITE(10,*) 'T', ' ', 'n_h1/n_h', ' ', 'n_he1/n_he', ' ', 'n_he2/n_he', ' ', 'ionization_energy'
       
  !K
  DO T = lower_T, upper_T, increment_T
    
    N_e = N_e_initial
    N_e_last = N_e_initial
    N_e1 = 0
          
    DO
      ratio_h(1) = const * (T)**(3.0/2.0) * & 
        (2.0 * B_h(2) * EXP(-X_h_eV(1)/(k_eV*T))) / (B_h(1) * N_e)
      ratio_he(1) = const * (T)**(3.0/2.0) * & 
        (2.0 * B_he(2) * EXP(-X_he_eV(1)/(k_eV*T))) / (B_he(1) * N_e)
      ratio_he(2) = const * (T)**(3.0/2.0) * & 
        (2.0 * B_he(3) * EXP(-X_he_eV(2)/(k_eV*T))) / (B_he(2) * N_e)
      
      y_h(1) = 1.0 / ((1.0/ratio_h(1)) + 1.0)
      y_he(1) = 1.0 / ((1.0/ratio_he(1)) + 1 + ratio_he(2))
      y_he(2) = 1.0 / ( (1.0 / ( ratio_he(1)*ratio_he(2) ))  + (1.0 / ratio_he(2) ) + 1.0)
      
      v_e(1) = 1.0*y_h(1)
      v_e(2) = 1.0*y_he(1) + 2.0*y_he(2)
            
      N_e1 =  ((rho*N_o) * (x(1)/A(1)) * v_e(1)) + & 
        ((rho*N_o) * (x(2)/A(2)) * v_e(2)) 
     
      IF (ABS((N_e - N_e1)/N_e) < tol) EXIT
        
      IF(ABS((N_e_last - N_e1)/N_e_last) > 1e-3) THEN 
        N_e_last = N_e
        N_e = N_e1
      ELSE 
        N_e_last = N_e
        N_e = (N_e + N_e1) / 2.0
      END IF
        
    END DO  
 
    WRITE(10,*) T,' ', y_h(1), ' ', y_he(1), ' ', y_he(2), ' ', &
      (y_h(1))* (rho*N_o*x(1)/A(1))*X_h_eV(1) + &
      (y_he(1))*(rho*N_o*x(2)/A(2))*X_he_eV(1) + &
      (y_he(2))*(rho*N_o*x(2)/A(2))*X_he_eV(2)

          
  END DO

END PROGRAM saha_helium

