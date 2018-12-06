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
      call RANDOM_SEED

      if (F_TYPE.eq.1) then
! **** LOW WAVENUMBER (FURUE, 2003) FORCING ****

      K0=7.

      DO J=0,TNKY
        DO K=0,TNKZ_S
          do I=0,NKX_S
            if ((KX_S(I).EQ.0) .AND. ((KZ_S(K).LT.0) .OR. 
     &            ((KZ_S(K).EQ.0) .AND. (KY(J).LT.0)))) then
            else
              if ((i+RANKZ*(NKX_S+1).LE.NKX) .AND.
     &            (k+RANKY*(TNKZ_S+1).LE.TNKZ)) then
                kappa2=KX2_S(I)+KZ2_S(K)
                K2=kappa2+KY2(J)
                IF ((K2.LE.100.) .AND. (KY(J).NE.0)
     &            .AND. (kappa2.NE.0)) THEN
                  call RANDOM_NUMBER(alpha)
                  alpha=2.d0*pi*alpha ! Random phase of forcing
                  CS1(I,K,J)=cexp(cmplx(0,alpha))*
     &                        K2**(1./4.)*(K2+K0**2)**(-3)
                  CF1(I,K,J)=CS1(i,k,j)*KY(j)*KX_S(i)/sqrt(K2*kappa2)
                  CF2(I,K,J)=-CS1(i,k,j)*sqrt(kappa2)/sqrt(K2)
                  CF3(I,K,J)=CS1(i,k,j)*KY(j)*KZ_S(k)/sqrt(K2*kappa2)
                  CSTH1(I,K,J)=CS1(i,k,j)*CI/sqrt(RI_TAU(1))
              puf=puf+conjg(CU1(I,K,J))*CF1(I,K,J)+conjg(CU2(I,K,J))
     &              *CF2(I,K,J)+conjg(CU3(I,K,J))*CF3(I,K,J)
     &              +RI_TAU(1)*conjg(CTH(I,K,J,1))*CSTH1(I,K,J)
                  if (FIRST_TIME) then
                  pff=pff+abs(CS1(I,K,J))**2
                  end if
                END IF
              end if
            end if
          end do
        END DO
      END DO

      CALL MPI_ALLREDUCE(puf,puf_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     &              MPI_COMM_WORLD,IERROR)
      if (FIRST_TIME) then
      CALL MPI_ALLREDUCE(pff,pff_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     &              MPI_COMM_WORLD,IERROR)
      end if
      puf_sum=2.d0*puf_sum
      if (FIRST_TIME) pff_sum=2.d0*pff_sum

      F01=(-puf_sum+SQRT(puf_sum**2+4*pff_sum*DELTA_T*P_AIM))
     &      /(2*DELTA_T*pff_sum)
      F02=(-puf_sum-SQRT(puf_sum**2+4*pff_sum*DELTA_T*P_AIM))
     &      /(2*DELTA_T*pff_sum)
      if (abs(F01)<abs(F02)) THEN
        DO J=0,TNKY
          DO K=0,TNKZ_S
            DO I=0,NKX_S
              CF1(I,K,J)=F01*CF1(I,K,J)
              CF2(I,K,J)=F01*CF2(I,K,J)
              CF3(I,K,J)=F01*CF3(I,K,J)
              CFTH(I,K,J,1)=CFTH(I,K,J,1)+F01*CSTH1(I,K,J)
            END DO
          END DO
        END DO
      else
        DO J=0,TNKY
          DO K=0,TNKZ_S
            DO I=0,NKX_S
              CF1(I,K,J)=F02*CF1(I,K,J)
              CF2(I,K,J)=F02*CF2(I,K,J)
              CF3(I,K,J)=F02*CF3(I,K,J)
              CFTH(I,K,J,1)=CFTH(I,K,J,1)+F02*CSTH1(I,K,J)
            END DO
          END DO
        END DO
      end if


      else if (F_TYPE.eq.2) then
! **** MAFFIOLI (2017) FORCING WITH CONSTANT POWER INJECTION ****

      DO J=0,TNKY
        DO K=0,TNKZ_S
          DO I=0,NKX_S
          if ((KX_S(I).EQ.0) .AND. ((KZ_S(K).LT.0) .OR. 
     &            ((KZ_S(K).EQ.0) .AND. (KY(J).LT.0)))) then
            else
              if ((i+RANKZ*(NKX_S+1).LE.NKX) .AND.
     &            (k+RANKY*(TNKZ_S+1).LE.TNKZ)) then
                kappa2=KX2_S(I)+KZ2_S(K)
                IF ((kappa2.gt.6.25) .and. (KY(J).EQ.0)
     &                .and. (kappa2.lt.12.25)) THEN
                  call RANDOM_NUMBER(alpha)
                  CF1(I,K,J)=KZ_S(K)*cexp(cmplx(0,2.d0*pi*alpha))
     &                          /sqrt(pi)/3**1.5
                  CF3(I,K,J)=-KX_S(I)*cexp(cmplx(0,2.d0*pi*alpha))
     &                          /sqrt(pi)/3**1.5
                  puf=puf+conjg(CU1(I,K,J))*CF1(I,K,J)
     &                      +conjg(CU3(I,K,J))*CF3(I,K,J)
                  if (FIRST_TIME) then
                pff=pff+0.5*(abs(CF1(I,K,J))**2+abs(CF3(I,K,J))**2)
                  end if
                END IF
              end if
            end if
          END DO
        END DO
      END DO

      CALL MPI_ALLREDUCE(puf,puf_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     &              MPI_COMM_WORLD,IERROR)
      if (FIRST_TIME) then
      CALL MPI_ALLREDUCE(pff,pff_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     &              MPI_COMM_WORLD,IERROR)
      end if
      puf_sum=2.d0*puf_sum
      if (FIRST_TIME) pff_sum=2.d0*pff_sum

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
      do j=0,TNKY
        do k=0,TNKZ_S
          do i=0,NKX_S
            if ((KX_S(I).EQ.0) .AND. ((KZ_S(K).LT.0) .OR. 
     &            ((KZ_S(K).EQ.0) .AND. (KY(J).LT.0)))) then
            else
              if ((i+RANKZ*(NKX_S+1).LE.NKX) .AND.
     &            (k+RANKY*(TNKZ_S+1).LE.TNKZ)) then
                kappa2=KX2_S(I)+KZ2_S(K)
                K2=kappa2+KY2(J)
                if ((K2.le.100.) .and. (KY(j).ne.0)
     &                .and. (kappa2.ne.0)) then
                  CS1(i,k,j)=cexp(cmplx(0,f_phase(i,k,j)-
     &                    sqrt(RI_TAU(1)*kappa2/K2)*TIME))*
     &                    (kappa2**3*K2)**(-0.5)
                  CF1(I,K,J)=CS1(i,k,j)*KY(j)*KX_S(i)/sqrt(K2*kappa2)
                  CF2(I,K,J)=-CS1(i,k,j)*sqrt(kappa2)/sqrt(K2)
                  CF3(I,K,J)=CS1(i,k,j)*KY(j)*KZ_S(k)/sqrt(K2*kappa2)
                  CSTH1(I,K,J)=CS1(i,k,j)*CI/sqrt(RI_TAU(1))
                  puf=puf+conjg(CU1(i,k,j))*CF1(i,k,j)
     &                +conjg(CU2(i,k,j))*CF2(i,k,j)
     &                +conjg(CU3(i,k,j))*CF3(i,k,j)
     &                +RI_TAU(1)*conjg(CTH(i,k,j,1))*CSTH1(i,k,j)
                  if (FIRST_TIME) then
                    pff=pff+abs(CS1(i,k,j))**2
                  end if
                end if
              end if
            end if
          end do
        end do
      end do

      call MPI_ALLREDUCE(puf,puf_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     &          MPI_COMM_WORLD,IERROR)
      if (FIRST_TIME) then
        call MPI_ALLREDUCE(pff,pff_sum,1,MPI_DOUBLE_PRECISION,
     &          MPI_SUM,MPI_COMM_WORLD,IERROR)
      end if
      puf_sum=2.d0*puf_sum
      if (FIRST_TIME) pff_sum=2.d0*pff_sum

      F01=(-puf_sum+SQRT(puf_sum**2+4*pff_sum*DELTA_T*P_AIM))
     &      /(2*DELTA_T*pff_sum)
      F02=(-puf_sum-SQRT(puf_sum**2+4*pff_sum*DELTA_T*P_AIM))
     &      /(2*DELTA_T*pff_sum)

      if (abs(F01)<abs(F02)) then
        do j=0,TNKY
          do k=0,TNKZ_S
            do i=0,NKX_S
              CF1(i,k,j)=F01*CF1(i,k,j)
              CF2(i,k,j)=F01*CF2(i,k,j)
              CF3(i,k,j)=F01*CF3(i,k,j)
              CFTH(i,k,j,1)=CFTH(i,k,j,1)+F01*CSTH1(i,k,j)
            end do
          end do
        end do
      else
        do j=0,TNKY
          do k=0,TNKZ_S
            do i=0,NKX_S
              CF1(i,k,j)=F02*CF1(i,k,j)
              CF2(i,k,j)=F02*CF2(i,k,j)
              CF3(i,k,j)=F02*CF3(i,k,j)
              CFTH(i,k,j,1)=CFTH(i,k,j,1)+F02*CSTH1(i,k,j)
            end do
          end do
        end do
      end if

      end if

      RETURN
      END
