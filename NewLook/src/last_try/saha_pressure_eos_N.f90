

PROGRAM saha_pressure_eos_N
  USE SahaPressureEOSDeclarationsN
  IMPLICIT NONE
  
  N_e_initial = N_e
  N_e_last = N_e

  OPEN(10, FILE='outfiles/saha_pressure_eos_N.dat')
  WRITE(10,*) 'T', ' ', 'n_h1/n_h', ' ', 'Pressure_Initial', ' ', 'Pressure_Calc'

  ! CGS
  DO P = lower_P, upper_P, increment_P
  
    N_e = N_e_initial
!    N_e_last = N_e_initial
    N_e1 = 0.0

    DO

      PRINT *, 'P ', P
      PRINT *, 'N_e ', N_e

      ! Calculate the Saha Equation for hydrogen 
      ratio_h(1) = const * (2.0 * B_h(2)/B_h(1)) * &
          ((1.0*P) / (n(1) + N_e))**(3.0/2.0) * &
          EXP((-X_h(1) * (n(1) + N_e) )/(1.0*P)) * (1.0/N_e) 
        
      ! Calculate the degree of ionization (aka ionization fraction)
      y_h(1) = 1.0 / ((1.0/ratio_h(1)) + 1.0)

      ! Calculate average electron contribution 
      v_e(1) = 1.0*y_h(1)
 
      G_N_e = n(1) * v_e(1) - N_e

      dRatio = const * (2.0 * B_h(2)/B_h(1)) * (1.0*P)**(3.0/2.0) * (n(1) + N_e)**(-5.0/2.0) * &
          EXP((-X_h(1) * (n(1) + N_e) )/(1.0*P)) * (1.0/N_e) * &
          ((-3.0/2.0) - (n(1) + N_e)*(X_h(1)/(1.0*P) + 1.0/N_e)) 
       
      dG_N_e = n(1) * ( 1.0 / (ratio_h(1) + 1))**2.0 * dRatio - 1.0

      PRINT *, '*****'
      PRINT *, G_N_e
      PRINT *, dRatio
      PRINT *, dG_N_e

!      ! Recompute electron density
!      N_e1 = n(1) * v_e(1) 
     
      ! Recompute electron density
      N_e1 = N_e - G_N_e/dG_N_e 
     
      IF (ABS((N_e - N_e1)/N_e) < tol) EXIT
      
      IF(ABS((N_e_last - N_e1)/N_e_last) > 1e-3) THEN
        PRINT *, 'HERE'
        N_e_last = N_e
        N_e = N_e1
      ELSE
        PRINT *, 'THERE'
        N_e_last = N_e
        N_e = (N_e + N_e1) / 2.0
      END IF
 
    END DO  

    ! Calculate the Saha Equation for hydrogen 
    ratio_h(1) = const * (2.0 * B_h(2)/B_h(1)) * &
        ((1.0*P) / (n(1) + N_e1))**(3.0/2.0) * &
        EXP((-X_h(1) * (n(1) + N_e1) )/(1.0*P)) * (1.0/N_e1) 

    ! Calculate the degree of ionization (aka ionization fraction)
    y_h(1) = 1.0 / ((1.0/ratio_h(1)) + 1.0)

    T = (1.0*P) / (k*(n(1) + N_e1)) 

    ! Log temperature, ionization fraction, pressure, internal energy, and specific heat
    WRITE(10,*) T,' ', y_h(1), ' ', P, ' ', Pressure(y_h(1), T), ' ', N_e1
          
  END DO

END PROGRAM saha_pressure_eos_N

