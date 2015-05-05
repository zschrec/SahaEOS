!


PROGRAM saha_hydrogen
  USE SahaSolverDeclarations
  IMPLICIT NONE
     
  N_e_initial = N_e
  N_e_last = N_e
 
  OPEN(10, FILE='outfiles/saha_hydrogen.dat')
  WRITE(10,*) 'T',' ', 'n1/n', ' ', 'ionization_energy'
       
  !K
  DO T = lower_T, upper_T, increment_T
    
    N_e = N_e_initial
    N_e_last = N_e_initial
    N_e1 = 0
          
    DO
      ratio_h(1) = const * (T)**(3.0/2.0) * & 
        (2.0 * B_h(2) * EXP(-X_h_eV(1)/(k_eV*T))) / (B_h(1) * N_e)
      
      y_h(1) = 1.0 / ((1.0/ratio_h(1)) + 1.0)
      v_e(1) = 1.0*y_h(1)
      N_e1 = (rho*N_o) * (x(1)/A(1)) * v_e(1)
      
      IF (ABS((N_e - N_e1)/N_e) < tol) EXIT
        
      IF(ABS((N_e_last - N_e1)/N_e_last) > 1e-3) THEN 
        N_e_last = N_e
        N_e = N_e1
      ELSE 
        N_e_last = N_e
        N_e = (N_e + N_e1) / 2.0
      END IF
        
    END DO  
 
    WRITE(10,*) T,' ', y_h(1), ' ', (y_h(1))* (rho*N_o*x(1)/A(1))*X_h_eV(1)

          
  END DO

END PROGRAM saha_hydrogen

