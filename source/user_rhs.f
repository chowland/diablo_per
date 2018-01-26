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

      call RANDOM_NUMBER(alpha)
      alpha=2.*pi*alpha ! Random phase of forcing
      K0=7.
      F_0=1.0           ! Amplitude of forcing

      DO J=0,TNKY
        DO K=0,TNKZ_S
          do I=0,NKX_S
            IF ( (KX2_S(I)+KZ2_S(K)+KY2(J).LE.100.) .AND.
     &            (KY(J).NE.0) .AND. (KX2_S(I)+KZ2_S(K).NE.0) ) THEN
              CS1(I,K,J)=cexp(cmplx(0,alpha))*(8.*pi)**(-0.5)*
     &                      (KX2_S(I)+KZ2_S(K)+KY2(J)+K0**2)**(-3)
              CF1(I,K,J)=CF1(I,K,J)+CS1(I,K,J)*KY(J)*KX_S(I)/
     &       (KX2_S(I)+KZ2_S(K)+KY2(J))**(1./4.)/SQRT(KX2_S(I)+KZ2_S(K))
              CF2(I,K,J)=CF2(I,K,J)+CS1(I,K,J)*KY(J)*KZ_S(K)/
     &       (KX2_S(I)+KZ2_S(K)+KY2(J))**(1./4.)/SQRT(KX2_S(I)+KZ2_S(K))
              CF3(I,K,J)=CF3(I,K,J)+CS1(I,K,J)*SQRT(KX2_S(I)+KZ2_S(K))/
     &           (KX2_S(I)+KZ2_S(K)+KY2(J))**(1./4.)
            END IF
          end do
        END DO
      END DO

      DO J=0,TNKY_TH
        DO K=0,TNKZ_S_TH
          do I=0,NKX_S_TH
          IF ((KX2_S_TH(I)+KZ2_S_TH(K)+KY2_TH(J).LE.100.).AND.
     &        (KY_TH(J).NE.0).AND.(KX2_S_TH(I)+KZ2_S_TH(K).NE.0)) THEN
              CSTH1(I,K,J)=cexp(cmplx(0,alpha))*(8.*pi)**(-0.5)*
     &                  (KX2_S_TH(I)+KZ2_S_TH(K)+KY2_TH(J)+K0**2)**(-3)
              CFTH(I,K,J,1)=CFTH(I,K,J,1)+CSTH1(I,K,J)*CI*
     &       (KX2_S_TH(I)+KZ2_S_TH(K)+KY2_TH(J))**(1./4.)
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
