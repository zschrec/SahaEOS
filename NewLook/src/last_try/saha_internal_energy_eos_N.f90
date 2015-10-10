

PROGRAM saha_internal_energy_eos_N
  USE SahaInternalEnergyEOSDeclarationsN
  IMPLICIT NONE
  
  N_e_initial = N_e
  N_e_last = N_e

  OPEN(10, FILE='outfiles/saha_internal_energy_eos_N.dat')
  WRITE(10,*) 'T', ' ', 'n_h1/n_h', ' ', 'Pressure_Initial', ' ', 'Pressure_Calc'

  ! CGS
  DO U = lower_U, upper_U, increment_U
  
    N_e = N_e_initial
    N_e_last = N_e_initial
    N_e1 = 0.0

    tol_count = 0

    DO

      PRINT *, 'U ', U
      PRINT *, 'N_e ', N_e

      IF (N_e == 0) THEN
        EXIT
      END IF
     
!      IF (tol_count == 2000) THEN
!        tol_count = 0
!      END IF

!      IF(1.0*U - N_e*X_h(1) < 0) THEN
!        PRINT *, '******************************** ZIP ', tol_count
!      !  N_e = N_e/(10.0**(LOG10(N_e)/2.0))
!        N_e = N_e/(10.0**(1.0+tol_count/1000))
!      END IF

      A_n = const * (2.0*B_h(2)/B_h(1)) * (1.5)**(-3.0/2.0)
      PRINT *, 'A_n ', A_n

      f_factor = 1.0*U - N_e*X_h(1)
      IF (f_factor < 0) THEN
!        !f_factor = 1.0*U - (N_e/(10.0**(LOG10(N_e)/2.0)))*X_h(1)
        T = (1.0*U - N_e_last*X_h(1)) / (k*(n(1) + N_e_last)) 
        f_factor = 1.5*N_e_last*k*T
      END IF
      B_n = ((f_factor)/(n(1)+N_e))**(3.0/2.0)

!!      B_n = ((1.0*U - N_e*X_h(1))/(n(1)+N_e))**(3.0/2.0)
!      PRINT *, 'B_n ', B_n

!      C_n = EXP( ( -1.5 * X_h(1) * (n(1) + N_e) ) / ( 1.0*U - N_e * X_h(1) ) )
      C_n = EXP( ( -1.5 * X_h(1) * (n(1) + N_e) ) / ( f_factor) )
      PRINT *, 'C_n ', C_n

      D_n = (1.0/N_e)
      PRINT *, 'D_n ', D_n

      dB_n = (-3.0/2.0) * ( (f_factor ) / (n(1) + N_e) )**(1.0/2.0) * &
          ( (1.0*U + N_e * X_h(1)) / (n(1) + N_e)**2.0 )
      
!!      dB_n = (-3.0/2.0) * ( (1.0*U - N_e*X_h(1) ) / (n(1) + N_e) )**(1.0/2.0) * &
!!          ( (1.0*U + N_e * X_h(1)) / (n(1) + N_e)**2.0 )
!      PRINT *, 'dB_n ', dB_n

!      dC_n = -1.5*X_h(1)*C_n* ( (1.0*U + n(1)*X_h(1)) / (1.0*U - N_e*X_h(1))**2.0 )
      dC_n = -1.5*X_h(1)*C_n* ( (1.0*U + n(1)*X_h(1)) / (f_factor)**2.0 )
      PRINT *, 'dC_n ', dC_n

      dD_n = -1.0/(N_e**2.0)
      PRINT *, 'dD_n ', dD_n

      ! Calculate the Saha Equation for hydrogen 
      ratio_h(1) = A_n*B_n*C_n*D_n
      PRINT *, 'ratio', ratio_h(1)

!      ratio_h(1) = const * (2.0 * B_h(2)/B_h(1)) * (1.5)**(-3.0/2.0) * &
!          ((1.0*U - N_e*X_h(1)) / (n(1) + N_e))**(3.0/2.0) * &
!          EXP((-1.5*X_h(1) * (n(1) + N_e) )/(1.0*U - N_e*X_h(1))) * (1.0/N_e) 
        
      ! Calculate the degree of ionization (aka ionization fraction)
      y_h(1) = 1.0 / ((1.0/ratio_h(1)) + 1.0)

      ! Calculate average electron contribution 
      v_e(1) = 1.0*y_h(1)
 
      G_N_e = n(1) * v_e(1) - N_e

      !dRatio = const * (2.0 * B_h(2)/B_h(1)) * (1.0*P)**(3.0/2.0) * (n(1) + N_e)**(-5.0/2.0) * &
      !    EXP((-X_h(1) * (n(1) + N_e) )/(1.0*P)) * (1.0/N_e) * &
      !    ((-3.0/2.0) - (n(1) + N_e)*(X_h(1)/(1.0*P) + 1.0/N_e)) 
      
      dRatio = A_n * ( dB_n*C_n*D_n + B_n*(dC_n*D_n + C_n*dD_n) )

!      dRatio = ratio_h(1) * ( (-3.0/2.0)*((1.0*U + X_h(1)*n(1))/((1.0*U-N_e*X_h(1))*(n(1)+N_e)))-&
!          1.5*X_h(1)*((1.0*U+n(1)*X_h(1))/(1.0*U-N_e*X_h(1))**2.0) - (1.0/N_e))

!      dRatio = const * (2.0 * B_h(2)/B_h(1)) * (1.5)**(-3.0/2.0) * &
!          ((1.0*U - N_e*X_h(1))/(n(1) - N_e))**(1.0/2.0) * &
!          EXP((-1.5*X_h(1) * (n(1) + N_e) )/(1.0*U - N_e*X_h(1))) * (1.0/N_e) * &
!          ((-3.0/2.0)*((1.0*U+n(1)*X_h(1))/((n(1)+N_e)**2.0)) + &
!          ((1.0*U-N_e*X_h(1))/(n(1)+N_e))* &
!          ( (-1.5*X_h(1)*(U+n(1)*X_h(1)))/((U-N_e*X_h(1))**2.0) - (1.0/N_e) ))
       
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
      

      IF(tol_count > 100 .AND. tolA == tol_m3) THEN
        tolA = tol_p3
      ELSE IF(tol_count > 100 .AND. tolA == tol_p3) THEN
        tolA = tol_m3
      END IF

      !IF(ABS((N_e_last - N_e1)/N_e_last) > 1e-3) THEN
      IF(ABS((N_e_last - N_e1)/N_e_last) > tolA) THEN
        PRINT *, 'HERE'
        N_e_last = N_e
        N_e = N_e1
      ELSE
        PRINT *, 'THERE'
        N_e_last = N_e
        N_e = (N_e + N_e1) / 2.0
      END IF
 
      tol_count = tol_count + 1

    END DO  

!    IF(N_e1 == 0) THEN
!      WRITE(10,*) 'A',' ', 'B', ' ', 'C', ' ', 'D', ' ', 'E'
!      CONTINUE
!    END IF

    N_e = N_e1

!    A_n = const * (2.0*B_h(2)/B_h(1)) * (1.5)**(-3.0/2.0)
      
!    B_n = ((f_factor)/(n(1)+N_e))**(3.0/2.0)

!!    B_n = ((1.0*U - N_e*X_h(1))/(n(1)+N_e))**(3.0/2.0)
!    C_n = EXP( ( -1.5 * X_h(1) * (n(1) + N_e) ) / ( 1.0*U - N_e * X_h(1) ) )
!    D_n = (1.0/N_e)

!    ratio_h(1) = A_n*B_n*C_n*D_n

    ratio_h(1) = const * (2.0 * B_h(2)/B_h(1)) * (1.5)**(-3.0/2.0) * &
        ((1.0*U - N_e*X_h(1)) / (n(1) + N_e))**(3.0/2.0) * &
        EXP((-1.5*X_h(1) * (n(1) + N_e) )/(1.0*U - N_e*X_h(1))) * (1.0/N_e) 

    ! Calculate the Saha Equation for hydrogen 
!    ratio_h(1) = const * (2.0 * B_h(2)/B_h(1)) * (1.5)**(-3.0/2.0) * &
!         ((1.0*U - N_e*X_h(1)) / (n(1) + N_e))**(3.0/2.0) * &
!         EXP((-1.5*X_h(1) * (n(1) + N_e) )/(1.0*U - N_e*X_h(1))) * (1.0/N_e) 

    ! Calculate the degree of ionization (aka ionization fraction)
    y_h(1) = 1.0 / ((1.0/ratio_h(1)) + 1.0)

    T = (1.0*U - N_e*X_h(1)) / (k*(n(1) + N_e)) 

    ! Log temperature, ionization fraction, pressure, internal energy, and specific heat
    WRITE(10,*) T,' ', y_h(1), ' ', U, ' ', InternalEnergy(y_h(1), T), ' ', N_e1
          
  END DO

END PROGRAM saha_internal_energy_eos_N

