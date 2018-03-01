      SUBROUTINE USER_RHS_PER_FOURIER

      include 'header'
! Here, you can add terms to the right hand side of the momentum
! and scalar equations.
! The right hand side foring arrays, CF1, CF2, CF3, CFTH are in Fourier space.
! The velocity and scalars are available in Fourier space.
! CS1 and CSTH1 are available as working variables.

      integer i,j,k,n

      real*8 alpha, K0, F_0

! **** LOW WAVENUMBER (FURUE, 2003) FORCING ****


      K0=7.
      F_0=1.0e3           ! Amplitude of forcing

      DO J=0,TNKY
        DO K=0,TNKZ_S
          do I=0,NKX_S
            IF ( (KX2_S(I)+KZ2_S(K)+KY2(J).LE.100.) .AND.
     &            (KY(J).NE.0) .AND. (KX2_S(I)+KZ2_S(K).NE.0) ) THEN
              call RANDOM_NUMBER(alpha)
              alpha=2.*pi*alpha ! Random phase of forcing
              CS1(I,K,J)=F_0*cexp(cmplx(0,alpha))*(4.*pi)**(-0.5)*
     &                      (KX2_S(I)+KZ2_S(K)+KY2(J))**(1./4.)*
     &                      (KX2_S(I)+KZ2_S(K)+KY2(J)+K0**2)**(-3)
              CF1(I,K,J)=CF1(I,K,J)+CS1(i,k,j)*KY(j)*KX_S(i)/sqrt(
     &          (KX2_S(i)+KZ2_S(k)+KY2(j))*(KX2_S(i)+KZ2_S(k)))
              CF2(I,K,J)=CF2(I,K,J)+CS1(i,k,j)*sqrt(KX2_S(i)+KZ2_S(k))
     &          /sqrt(KX2_S(i)+KZ2_S(k)+KY2(j))
              CF3(I,K,J)=CF3(I,K,J)+CS1(i,k,j)*KY(j)*KZ_S(k)/sqrt(
     &          (KX2_S(i)+KZ2_S(k)+KY2(j))*(KX2_S(i)+KZ2_S(k)))
              CFTH(I,K,J,1)=CFTH(I,K,J,1)+CS1(i,k,j)*CI/sqrt(RI_TAU(1))
            END IF
          end do
        END DO
      END DO



! For example, to add a linear damping term (e.g. -alpha*U) to the RHS:
!       alpha=-0.1d0
!       DO J=JSTART,JEND
!         DO K=0,TNKZ
!           DO I=0,NXP-1
!             CF1(I,K,J)=CF1(I,K,J)-alpha*CU1(I,K,J)
!           END DO
!         END DO
!       END DO

! For U2 do this...
! Note that the only thing that changes are the bounds of the J index
!       DO J=2,NY
!         DO K=0,TNKZ
!           DO I=0,NXP-1
!             CF2(I,K,J)=CF2(I,K,J)-alpha*CU2(I,K,J)
!           END DO
!         END DO
!       END DO

! For scalars, do this...
!       DO J=JSTART,JEND
!         DO K=0,TNKZ
!           DO I=0,NXP-1
!             CFTH(I,K,J,N)=CFTH(I,K,J,N)-alpha*CTH(I,K,J,N)
!           END DO
!         END DO
!       END DO
!      END DO

      RETURN
      END
