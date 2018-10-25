      SUBROUTINE USER_RHS_PER_FOURIER

      include 'header'
! Here, you can add terms to the right hand side of the momentum
! and scalar equations.
! The right hand side foring arrays, CF1, CF2, CF3, CFTH are in Fourier space.
! CF1, CF2, CF3 are available as working variables, but also define the output.
! CFTH already contains the background stratification, so is not a working variable.
! The velocity and scalars are available in Fourier space.
! CS1 and CSTH1 are available as working variables.

      integer i,j,k,n

      real*8 alpha, K0, F0, B2, puf, pff, puf_sum, PF, PF_sum
      real*8 P_AIM, F01, F02, kappa2, K2, PS, PS_sum, FS1, FS2

      P_AIM=RI_TAU(1)*target_Reb*NU
      puf=0.d0
      PS=0.d0
      PF=0.d0
      if (FIRST_TIME) pff=0.d0

      if (F_TYPE.eq.1) then
! **** LOW WAVENUMBER (FURUE, 2003) FORCING ****

      K0=7.

      DO J=0,TNKY
        DO K=0,TNKZ_S
          do I=0,NKX_S
            kappa2=KX2_S(I)+KZ2_S(K)
            K2=kappa2+KY2(J)
            IF ((K2.LE.100.) .AND. (KY(J).NE.0)) THEN
              IF (kappa2.EQ.0) THEN
                IF (FORCE_SHEAR) THEN
                CS1(I,K,J)=(4.*pi)**-0.5*K2**0.25*(K2+K0**2)**-3
                call RANDOM_NUMBER(alpha)
                alpha=2.d0*pi*alpha
                call RANDOM_NUMBER(F0)
                CF1(I,K,J)=sqrt(F0)*CS1(I,K,J)*cexp(cmplx(0,alpha))
                call RANDOM_NUMBER(alpha)
                alpha=2.d0*pi*alpha
                CF3(I,K,J)=sqrt(1-F0)*CS1(I,K,J)*cexp(cmplx(0,alpha))
                PS=PS+conjg(CU1(I,K,J))*CF1(I,K,J)+
     &              conjg(CU3(I,K,J))*CF3(I,K,J)
                PF=PF+0.5*abs(CS1(I,K,J))**2
                END IF
              ELSE
              call RANDOM_NUMBER(alpha)
              alpha=2.d0*pi*alpha ! Random phase of forcing
              CS1(I,K,J)=cexp(cmplx(0,alpha))*(4.*pi)**(-0.5)*
     &                      K2**(1./4.)*(K2+K0**2)**(-3)
              CF1(I,K,J)=CS1(i,k,j)*KY(j)*KX_S(i)/sqrt(K2*kappa2)
              CF2(I,K,J)=CS1(i,k,j)*sqrt(kappa2)/sqrt(K2)
              CF3(I,K,J)=CS1(i,k,j)*KY(j)*KZ_S(k)/sqrt(K2*kappa2)
              CSTH1(I,K,J)=CS1(i,k,j)*CI/sqrt(RI_TAU(1))
              puf=puf+conjg(CU1(I,K,J))*CF1(I,K,J)+conjg(CU2(I,K,J))
     &            *CF2(I,K,J)+conjg(CU3(I,K,J))*CF3(I,K,J)
     &            +RI_TAU(1)*conjg(CTH(I,K,J,1))*CSTH1(I,K,J)
              if (FIRST_TIME) then
              pff=pff+abs(CS1(I,K,J))**2
              end if
              END IF
            END IF
          end do
        END DO
      END DO

      CALL MPI_ALLREDUCE(puf,puf_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     &              MPI_COMM_WORLD,IERROR)
      if (FIRST_TIME) then
      CALL MPI_ALLREDUCE(pff,pff_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     &              MPI_COMM_WORLD,IERROR)
      end if
      if (FORCE_SHEAR) then
      CALL MPI_ALLREDUCE(PS,PS_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     &              MPI_COMM_WORLD,IERROR)
      CALL MPI_ALLREDUCE(PF,PF_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     &              MPI_COMM_WORLD,IERROR)
      end if

      F01=(-puf_sum+SQRT(puf_sum**2+4*pff_sum*DELTA_T*P_AIM))
     &      /(2*DELTA_T*pff_sum)
      F02=(-puf_sum-SQRT(puf_sum**2+4*pff_sum*DELTA_T*P_AIM))
     &      /(2*DELTA_T*pff_sum)
      if (FORCE_SHEAR) then
        FS1=(-PS_sum+SQRT(PS_sum**2+4*PF_sum*DELTA_T*P_AIM))
     &      /(2*DELTA_T*PF_sum)
        FS2=(-PS_sum-SQRT(PS_sum**2+4*PF_sum*DELTA_T*P_AIM))
     &      /(2*DELTA_T*PF_sum)
      end if
      if (abs(F01)<abs(F02)) THEN
        DO J=0,TNKY
          DO K=0,TNKZ_S
            DO I=0,NKX_S
            kappa2=KX2_S(I)+KZ2_S(K)
            if (kappa2.eq.0) then
              if (FORCE_SHEAR) then
                if (abs(FS1)<abs(FS2)) then
                  CF1(I,K,J)=FS1*CF1(I,K,J)
                  CF3(I,K,J)=FS1*CF3(I,K,J)
                else
                CF1(I,K,J)=FS2*CF1(I,K,J)
                CF3(I,K,J)=FS2*CF3(I,K,J)
                end if
              end if
            else
              CF1(I,K,J)=F01*CF1(I,K,J)
              CF2(I,K,J)=F01*CF2(I,K,J)
              CF3(I,K,J)=F01*CF3(I,K,J)
              CFTH(I,K,J,1)=CFTH(I,K,J,1)+F01*CSTH1(I,K,J)
            end if
            END DO
          END DO
        END DO
      else
        DO J=0,TNKY
          DO K=0,TNKZ_S
            DO I=0,NKX_S
            kappa2=KX2_S(I)+KZ2_S(K)
            if (kappa2.eq.0) then
              if (FORCE_SHEAR) then
                if (abs(FS1)<abs(FS2)) then
                CF1(I,K,J)=FS1*CF1(I,K,J)
                CF3(I,K,J)=FS1*CF3(I,K,J)
                else
                CF1(I,K,J)=FS2*CF1(I,K,J)
                CF3(I,K,J)=FS2*CF3(I,K,J)
                end if
              end if
            else
              CF1(I,K,J)=F02*CF1(I,K,J)
              CF2(I,K,J)=F02*CF2(I,K,J)
              CF3(I,K,J)=F02*CF3(I,K,J)
              CFTH(I,K,J,1)=CFTH(I,K,J,1)+F02*CSTH1(I,K,J)
            end if
            END DO
          END DO
        END DO
      end if


      else if (F_TYPE.eq.2) then
! **** MAFFIOLI (2017) FORCING WITH CONSTANT POWER INJECTION ****

      DO J=0,TNKY
        DO K=0,TNKZ_S
          DO I=0,NKX_S
            kappa2=KX2_S(I)+KZ2_S(K)
            IF ((kappa2.gt.6.25) .and. (KY(J).EQ.0)
     &        .and. (kappa2.lt.12.25)) THEN
              call RANDOM_NUMBER(alpha)
              CF1(I,K,J)=KZ_S(K)*cexp(cmplx(0,2.d0*pi*alpha))
     &                  /sqrt(pi)/3**1.5
              CF3(I,K,J)=-KX_S(I)*cexp(cmplx(0,2.d0*pi*alpha))
     &                  /sqrt(pi)/3**1.5
              puf=puf+conjg(CU1(I,K,J))*CF1(I,K,J)
     &              +conjg(CU3(I,K,J))*CF3(I,K,J)
              if (FIRST_TIME) then
              pff=pff+0.5*(abs(CF1(I,K,J))**2+abs(CF3(I,K,J))**2)
              end if
            END IF
          END DO
        END DO
      END DO

      CALL MPI_ALLREDUCE(puf,puf_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     &              MPI_COMM_WORLD,IERROR)
      if (FIRST_TIME) then
      CALL MPI_ALLREDUCE(pff,pff_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     &              MPI_COMM_WORLD,IERROR)
      end if

      F01=(-puf_sum+SQRT(puf_sum**2+4*pff_sum*DELTA_T*P_AIM))
     &      /(2*DELTA_T*pff_sum)
      F02=(-puf_sum-SQRT(puf_sum**2+4*pff_sum*DELTA_T*P_AIM))
     &      /(2*DELTA_T*pff_sum)
      if (abs(F01)<abs(F02)) THEN
        DO J=0,TNKY
          DO K=0,TNKZ_S
            DO I=0,NKX_S
              CF1(I,K,J)=F01*CF1(I,K,J)
              CF3(I,K,J)=F01*CF3(I,K,J)
            END DO
          END DO
        END DO
      else
        DO J=0,TNKY
          DO K=0,TNKZ_S
            DO I=0,NKX_S
              CF1(I,K,J)=F02*CF1(I,K,J)
              CF3(I,K,J)=F02*CF3(I,K,J)
            END DO
          END DO
        END DO
      end if

      else if (F_TYPE.EQ.3) then
! **** PROPAGATING WAVE FORCING ****
      
      end if

      RETURN
      END
