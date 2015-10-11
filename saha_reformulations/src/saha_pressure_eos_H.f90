

PROGRAM saha_pressure_eos_H
  USE SahaPressureEOSDeclarationsH
  IMPLICIT NONE
  
  N_e_initial = N_e
  N_e_last = N_e

  OPEN(10, FILE='outfiles/saha_pressure_eos_H.dat')
  WRITE(10,*) 'T', ' ', 'n_h1/n_h'

  ! CGS
  DO P = lower_P, upper_P, increment_P
  
    N_e = N_e_initial
    N_e_last = N_e_initial
    N_e1 = 0.0

    here = 0
    there = 0

    DO

      PRINT *, 'P ', P
      PRINT *, 'N_e ', N_e

      A_h(1) = const * (2.0 * B_h(2) / B_h(1)) * (1.0*P)**(3.0/2.0)

      B_Ne = (SUM(n) + N_e)**(-3.0/2.0)
      
      C_h_Ne(1) = EXP((-X_h(1)*(SUM(n) + N_e))/(1.0*P))
      
      D_Ne = 1.0/N_e

      ratio_h(1) = A_h(1)*B_Ne*C_h_Ne(1)*D_Ne

      dB_Ne = (-3.0/2.0) * (SUM(n) + N_e)**(-5.0/2.0)
      dC_h_Ne(1) = C_h_Ne(1) * (-X_h(1)/(1.0*P))
      dD_Ne = -1.0/(N_e**2.0)

      dratio_h(1) = A_h(1) * &
          (dB_Ne*C_h_Ne(1)*D_Ne + B_Ne*(dC_h_Ne(1)*D_Ne + C_h_Ne(1)*dD_Ne))

      ! Calculate the degree of ionization (aka ionization fraction)
      y_h(1) = 1.0 / ((1.0/ratio_h(1)) + 1.0)
      
      ! Calculate average electron contribution 
      v_e(1) = 1.0*y_h(1)
            
      ! Recompute electron density
      N_e_h = n(1) * v_e(1) 

      dN_e_h = n(1) * (ratio_h(1) + 1.0)**(-2.0) * dratio_h(1)

      G_N_e = N_e_h - N_e
      dG_N_e = dN_e_h - 1.0

      N_e1 = N_e - G_N_e/dG_N_e 

      IF (ABS((N_e - N_e1)/N_e) < tol) EXIT
      
      !IF(ABS((N_e_last - N_e1)/N_e_last) > 1e-3) THEN
      IF(ABS((N_e_last - N_e1)/N_e_last) > 10**(0 - here + there)) THEN
        PRINT *, 'HERE'
        N_e_last = N_e
        N_e = N_e1
        there = 0
        here = here + 1
      ELSE
        PRINT *, 'THERE'
        N_e_last = N_e
        N_e = (N_e + N_e1) / 2.0
        here = 0
        there = there + 1
      END IF

      !IF (here == 1000) here = 0
      !IF (there == 1000) there = 0

    END DO  

    N_e = N_e1

    A_h(1) = const * (2.0 * B_h(2) / B_h(1)) * (1.0*P)**(3.0/2.0)
    B_Ne = (SUM(n) + N_e)**(-3.0/2.0)
    C_h_Ne(1) = EXP((-X_h(1)*(SUM(n) + N_e))/(1.0*P))
    D_Ne = 1.0/N_e

    ratio_h(1) = A_h(1)*B_Ne*C_h_Ne(1)*D_Ne

    ! Calculate the degree of ionization (aka ionization fraction)
    y_h(1) = 1.0 / ((1.0/ratio_h(1)) + 1.0)

    T = (1.0*P) / (k*(SUM(n) + N_e)) 

    ! Log temperature, ionization fraction, pressure, internal energy, and specific heat
    WRITE(10,*) T,' ', y_h(1)
          
  END DO

END PROGRAM saha_pressure_eos_H

