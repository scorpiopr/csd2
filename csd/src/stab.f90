MODULE STAB_CLASS
USE FORWARD_CLASS
USE HOVER_CLASS
IMPLICIT NONE
!->+++++++++++++++++++++++++++++++++++++++++++++
TYPE,PUBLIC :: STABSIS

    INTEGER :: NDOFC,NDOFC2,NMN,NMN2
    real(rdt),ALLOCATABLE :: at(:,:)
    real(rdt),ALLOCATABLE :: VEGRS(:),VEGIS(:)

    real(rdt),ALLOCATABLE :: atm(:,:)
    real(rdt),ALLOCATABLE :: VEGRSm(:),VEGISm(:)
    
    CONTAINS

    PROCEDURE,PUBLIC :: SOLVS
    PROCEDURE,PUBLIC :: SOLVSM
    PROCEDURE,PUBLIC :: CONSTRUCT_HOVSTAB
    PROCEDURE,PUBLIC :: getmatrixa
    PROCEDURE,PUBLIC :: getmatrixaM

    PROCEDURE,PUBLIC :: OUTPUT_STAB
    PROCEDURE,PRIVATE :: eigensolvs
    PROCEDURE,PRIVATE :: FVNLIN2
    PROCEDURE,PRIVATE :: fvnlin2m

END TYPE STABSIS

CONTAINS
SUBROUTINE CONSTRUCT_HOVSTAB(THIS,HOV,INPU) 
IMPLICIT NONE
CLASS(STABSIS) :: THIS
CLASS(HOVER) :: HOV
TYPE(INPUT) :: INPU

    THIS%NMN=HOV%BLADE%NMN
    THIS%NMN2=HOV%BLADE%NMN*2
    THIS%NDOFC=HOV%BLADE%NDOFC
    THIS%NDOFC2=HOV%BLADE%NDOFC*2

    ALLOCATE(THIS%at(THIS%NDOFC2,THIS%NDOFC2))
    ALLOCATE(THIS%VEGRS(THIS%NDOFC))
    ALLOCATE(THIS%VEGIS(THIS%NDOFC))

    ALLOCATE(THIS%atm(THIS%NMN2,THIS%NMN2))
    ALLOCATE(THIS%VEGRSm(THIS%NMN))
    ALLOCATE(THIS%VEGISm(THIS%NMN))


END SUBROUTINE


!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
subroutine getmatrixa(THIS,HOV) 
implicit none
CLASS(STABSIS) :: THIS
CLASS(HOVER) :: HOV
integer err
integer i,j,cho
real(rdt)   dlxi
real(rdt), allocatable:: xv(:), fxv(:) ,jacb(:, :), xvtp(:), fvtp(:)
real(rdt), allocatable:: UMMCE(:,:),UKMCE(:,:),UCMCE(:,:)
REAL(RDT),ALLOCATABLE :: MMTINV(:,:),TP(:,:),EYES(:,:)
real(rdt) :: rlxv
integer :: ndofc
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> main body:
    
    rlxv=0.0001
    ndofc=THIS%NDOFC

      ALLOCATE(xv(ndofc), fxv(ndofc), jacb(ndofc, ndofc), xvtp(ndofc), fvtp(ndofc), &
            UMMCE(ndofc,ndofc),UKMCE(ndofc,ndofc),UCMCE(ndofc,ndofc),&
            EYES(ndofc,ndofc),MMTINV(ndofc,ndofc),TP(ndofc,ndofc),&
            STAT=ERR)
                  
	IF(ERR.NE.0) THEN
		WRITE(*, 998)
		STOP
	END IF
	
    call THIS%fvnlin2(HOV,fxv) 
	
    do j=1,3
        cho=j
         if( CHO .eq. 1 ) then
              xv=HOV%BLADE%UQS
          else if( CHO .eq. 2 ) then
              xv=HOV%BLADE%UQS1T
          else if( CHO .eq. 3 ) then
               xv=HOV%BLADE%UQS2T
          end if
          
          do 5 i=1, NDOFC
                xvtp=xv 
                if(ABS(xv(i)) .lt.1.0d-8) then
                    dlxi=1.0d-9   
                else
                    dlxi=rlxv*xv(i) 
                end if
                if( ABS(dlxi) .lt. 1.0d-15 ) then
                      jacb(:,i)=0.0D0
                      pause
                else
                    xvtp(i)=xvtp(i)+dlxi
                          if( CHO .eq. 1 ) then
                                HOV%BLADE%UQS=xvtp
                          else if( CHO .eq. 2 ) then
                              HOV%BLADE%UQS1T=xvtp
                          else if( CHO .eq. 3 ) then
                               HOV%BLADE%UQS2T=xvtp
                          end if
                    call THIS%fvnlin2(HOV,fvtp) 
                     
                    jacb(:, i)=(fvtp-fxv)/dlxi 
                end if
5       continue

         if( CHO .eq. 1 ) then
            UKMCE=JACB
          else if( CHO .eq. 2 ) then
            UCMCE=JACB
          else if( CHO .eq. 3 ) then
            UMMCE=JACB
          end if
    end do

        EYES=0.0
        DO I=1,NDOFC
            DO J=I,I
                EYES(I,J)=1.0
            END DO
        END DO    
    
        THIS%AT=0.0
        THIS%AT(1:NDOFC,NDOFC+1:2*NDOFC)=EYES
           
        call mtxinv2(MMTINV, 1, NDOFC, UMMCE)  
        TP=MATMUL(MMTINV,UKMCE)
        THIS%AT(NDOFC+1:2*NDOFC,1:NDOFC)=-TP
        
        TP=MATMUL(MMTINV,UCMCE)
        THIS%AT(NDOFC+1:2*NDOFC,NDOFC+1:2*NDOFC)=-TP


      DEALLOCATE(xv,fxv, jacb,xvtp,fvtp,UMMCE,UKMCE,UCMCE,&
                  EYES,MMTINV,&
            STAT=ERR)
                  
	IF(ERR.NE.0) THEN
		WRITE(*, 999)
		STOP
	END IF

998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
end subroutine getmatrixa
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
subroutine getmatrixam(THIS,HOV) 
implicit none
CLASS(STABSIS) :: THIS
CLASS(HOVER) :: HOV
integer err
integer i,j,cho
real(rdt)   dlxi
real(rdt), allocatable:: xv(:), fxv(:) ,jacb(:, :), xvtp(:), fvtp(:)
real(rdt), allocatable:: UMMCE(:,:),UKMCE(:,:),UCMCE(:,:)
REAL(RDT),ALLOCATABLE :: MMTINV(:,:),TP(:,:),EYES(:,:)
REAL(RDT),ALLOCATABLE ::QSP(:),Q1TSP(:),Q2TSP(:),&
            q(:),q1t(:),q2t(:)
real(rdt) :: rlxv
integer :: ndofc,nmn
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> main body:
    
    rlxv=0.001
    ndofc=THIS%NDOFC
    nmn=THIS%nmn

      ALLOCATE(xv(NMN), fxv(NMN), jacb(NMN, NMN), xvtp(NMN), fvtp(NMN), &
            UMMCE(NMN,NMN),UKMCE(NMN,NMN),UCMCE(NMN,NMN),&
            EYES(NMN,NMN),MMTINV(NMN,NMN),TP(NMN,NMN),&
            QSP(NMN),Q1TSP(NMN),Q2TSP(NMN),&
            q(NMN),q1t(NMN),q2t(NMN),&
            STAT=ERR)
                  
	IF(ERR.NE.0) THEN
		WRITE(*, 998)
		STOP
	END IF

!    CALL MVMU2(HOV%BLADE%UQS, THIS%NDOFC, THIS%NMN, THIS%VTEG, HOV%BLADE%q)
!    CALL MVMU2(HOV%BLADE%UQS1T, THIS%NDOFC, THIS%NMN, THIS%VTEG, HOV%BLADE%q1t)
!    CALL MVMU2(HOV%BLADE%UQS2T, THIS%NDOFC, THIS%NMN, THIS%VTEG, HOV%BLADE%q2t)

    QSP=HOV%BLADE%Q
    Q1TSP=HOV%BLADE%Q1T
    Q2TSP=HOV%BLADE%Q2T

!    QSP=0.0D0
!    Q1TSP=0.0D0
!    Q2TSP=0.0D0

    call THIS%fvnlin2m(HOV,fxv) 

    do j=1,3
        cho=j
         if( CHO .eq. 1 ) then
!    HOV%aerof%bsc=0.07853981635D0
              xv=QSP
          else if( CHO .eq. 2 ) then
!HOV%aerof%bsc=0.07853981635D0
              xv=Q1TSP
          else if( CHO .eq. 3 ) then
!HOV%aerof%bsc=0.25
               xv=Q2TSP
          end if

          do 5 i=1, NMN
                xvtp=xv 
                if(ABS(xv(i)) .lt.1.0d-8) then
                    dlxi=1.0d-7   
                else
                    dlxi=rlxv*xv(i) 
                end if
                if( ABS(dlxi) .lt. 1.0d-15 ) then
                      jacb(:,i)=0.0D0
                      pause
                else
                    xvtp(i)=xvtp(i)+dlxi
                          if( CHO .eq. 1 ) then
                                HOV%BLADE%q=xvtp
                                HOV%BLADE%q1t=Q1TSP
                                HOV%BLADE%q2t=Q2TSP
                          else if( CHO .eq. 2 ) then
                               HOV%BLADE%q=QSP
                                HOV%BLADE%q1t=xvtp
                                HOV%BLADE%q2t=Q2TSP
                          else if( CHO .eq. 3 ) then
                               HOV%BLADE%q=QSP
                                HOV%BLADE%q1t=Q1TSP
                                HOV%BLADE%q2t=xvtp
                          end if

                    call THIS%fvnlin2m(HOV,fvtp) 
                    jacb(:, i)=(fvtp)/dlxi 
                end if
5       continue

         if( CHO .eq. 1 ) then
            UKMCE=JACB
          else if( CHO .eq. 2 ) then
            UCMCE=JACB
          else if( CHO .eq. 3 ) then
            UMMCE=JACB
          end if
    end do
!    WRITE(*,*)         UKMCE
!    PAUSE
!    WRITE(*,*)         UCMCE
!    PAUSE
!    WRITE(*,*)         UMMCE
!    PAUSE

        EYES=0.0
        DO I=1,NMN
            DO J=I,I
                EYES(I,J)=1.0
            END DO
        END DO    
    
        THIS%ATM=0.0
        THIS%ATM(1:NMN,NMN+1:2*NMN)=EYES

        call mtxinv2(MMTINV, 1, NMN, UMMCE)  
        TP=MATMUL(MMTINV,UKMCE)
        THIS%ATM(NMN+1:2*NMN,1:NMN)=-TP
        
        TP=MATMUL(MMTINV,UCMCE)
        THIS%ATM(NMN+1:2*NMN,NMN+1:2*NMN)=-TP

      DEALLOCATE(xv,fxv, jacb,xvtp,fvtp,UMMCE,UKMCE,UCMCE,&
                  EYES,MMTINV,&
            STAT=ERR)
                  
	IF(ERR.NE.0) THEN
		WRITE(*, 999)
		STOP
	END IF

998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
end subroutine getmatrixam
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE OUTPUT_STAB(THIS,HOV)
IMPLICIT NONE
CLASS(STABSIS) :: THIS
CLASS(HOVER) :: HOV

INTEGER :: I,J
    open(111, file='hovereigi.dat',position='append')
    WRITE(111,997)  HOV%th0,THIS%vegiSM(1:7)/HOV%OMG      
    CLOSE(111)
    OPEN(111,FILE='hovereigr.dat',position='append')
    WRITE(111,997)  HOV%th0,THIS%VEGRSM(1:7)/HOV%OMG         
    CLOSE(111)

!    open(111, file='hovereigi.dat',position='append')
!    WRITE(111,997)  HOV%th0,THIS%vegiS(1:7)/HOV%OMG      
!    CLOSE(111)
!    OPEN(111,FILE='hovereigr.dat',position='append')
!    WRITE(111,997)  HOV%th0,THIS%VEGRS(1:7)/HOV%OMG         
!    CLOSE(111)

!    open(111, file='hovereigi.dat',position='append')
!    WRITE(111,997)  -HOV%BLADE%LDTIP(3)/PI*180.0D0,THIS%vegiS(1:8)/HOV%OMG      
!    CLOSE(111)
!    OPEN(111,FILE='hovereigr.dat',position='append')
!    WRITE(111,997)  -HOV%BLADE%LDTIP(3)/PI*180.0D0,THIS%VEGRS(1:8)/HOV%OMG         
!    CLOSE(111)

    OPEN(111,FILE='stabMATRIX.dat')
    DO I=1,THIS%NMN2
        WRITE(111,997)  (THIS%ATM(I,J),J=1,THIS%NMN2)
    END DO
    CLOSE(111)

997 format(1000f40.12)
END SUBROUTINE
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
subroutine fvnlin2(THIS,HOV,fvtp)  
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
implicit none
CLASS(STABSIS) :: THIS
CLASS(HOVER) :: HOV
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
real(rdt), intent(out):: fvtp(THIS%NDOFC)

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->        
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
integer SOL
REAL(RDT) :: PPSI,TH1C,TH1S
REAL(RDT) :: TP1(THIS%NDOFC),TP2(THIS%NDOFC),TP3(THIS%NDOFC)

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

    SOL=1
    PPSI=0.0D0
    TH1C=0.0D0
    TH1S=0.0D0
        
    CALL HOV%BLADE%UPDATE()
    CALL HOV%AEROF%MLEGGAUSS_AERO(HOV%OMG,HOV%BLADE,HOV%TH0,TH1C,TH1S,PPSI,SOL)
    CALL HOV%BLADE%MLEGGAUSS_ST_BEAM(HOV%TH0,TH1C,TH1S,PPSI,HOV%OMG,SOL)
    CALL HOV%BLADE%FORMULATE_ST()  

    TP1=MATMUL(HOV%BLADE%UMMC,HOV%BLADE%UQS2T)
    TP2=MATMUL(HOV%BLADE%UCMC,HOV%BLADE%UQS1T)
    TP3=MATMUL(HOV%BLADE%UKMC,HOV%BLADE%UQS) 
    
    
      fvtp(1:THIS%NDOFC)=TP1+TP2+TP3+HOV%BLADE%UFVC
      
end subroutine fvnlin2
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
subroutine fvnlin2m(THIS,HOV,fvtp)  
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
implicit none
CLASS(STABSIS) :: THIS
CLASS(HOVER) :: HOV
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
real(rdt), intent(out):: fvtp(THIS%NMN)
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->        
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
integer SOL
REAL(RDT) :: PPSI,TH1C,TH1S
REAL(RDT) :: TP1(THIS%NMN),TP2(THIS%NMN),TP3(THIS%NMN)
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

    CALL MVMU2(HOV%BLADE%UQS, HOV%BLADE%NDOFC, HOV%BLADE%NMN, HOV%BLADE%VTEG, HOV%BLADE%q)
    CALL MVMU2(HOV%BLADE%UQS1T, HOV%BLADE%NDOFC, HOV%BLADE%NMN, HOV%BLADE%VTEG, HOV%BLADE%q1t)
    CALL MVMU2(HOV%BLADE%UQS2T, HOV%BLADE%NDOFC, HOV%BLADE%NMN, HOV%BLADE%VTEG, HOV%BLADE%q2t)
    
    SOL=1
    PPSI=0.0D0
    TH1C=0.0D0
    TH1S=0.0D0
        
    CALL HOV%BLADE%UPDATE()
    CALL HOV%AEROF%MLEGGAUSS_AERO(HOV%OMG,HOV%BLADE,HOV%TH0,TH1C,TH1S,PPSI,SOL)
    CALL HOV%BLADE%MLEGGAUSS_ST_BEAM(HOV%TH0,TH1C,TH1S,PPSI,HOV%OMG,SOL)
    CALL HOV%BLADE%FORMULATE_ST()  

    CALL HOV%BLADE%GETMKCFT()

    TP1=MATMUL(HOV%BLADE%mmt,HOV%BLADE%q2t)
    TP2=MATMUL(HOV%BLADE%cmt,HOV%BLADE%q1t)
    TP3=MATMUL(HOV%BLADE%kmt,HOV%BLADE%q) 
    
      fvtp(1:THIS%nmn)=TP1+TP2+TP3-HOV%BLADE%fv
!      write(*,*) fvtp
!      pause
end subroutine fvnlin2m
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE SOLVS(THIS)
IMPLICIT NONE
CLASS(STABSIS) :: THIS
INTEGER :: IMK

        CALL THIS%eigensolvs(THIS%NDOFC2, THIS%AT, THIS%NDOFC, THIS%VEGRS, THIS%VEGIS, IMK)

END SUBROUTINE
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE SOLVSM(THIS)
IMPLICIT NONE
CLASS(STABSIS) :: THIS
INTEGER :: IMK

        CALL THIS%eigensolvs(THIS%NMN2, THIS%ATM, THIS%NMN, THIS%VEGRSM, THIS%VEGISM, IMK)

END SUBROUTINE
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine eigensolvs(THIS,NN, AK,  NMD, VEGR, VEGI, IMK)
implicit none
CLASS(STABSIS) :: THIS
INTEGER, INTENT(IN):: NN, NMD
REAL(RDT), INTENT(IN):: AK(NN, NN)
REAL(RDT), INTENT(OUT):: VEGR(NMD), VEGI(NMD)
INTEGER, INTENT(OUT):: IMK


CHARACTER BALANC, JOBVL, JOBVR, SENSE
INTEGER N, LDA
INTEGER LDVL, LDVR, ILO, IHI
INTEGER LWORK,MARK
INTEGER INFO
DOUBLE PRECISION ABNRM, BBNRM
DOUBLE PRECISION , ALLOCATABLE:: A(:, :)
DOUBLE PRECISION ,  ALLOCATABLE:: WI(:), WR(:)
DOUBLE PRECISION , ALLOCATABLE:: VL(:, :), VR(:, :),SCALE(:)
DOUBLE PRECISION , ALLOCATABLE:: RCONDE(:), RCONDV(:)
DOUBLE PRECISION , ALLOCATABLE:: WORK(:)
INTEGER, ALLOCATABLE:: IWORK(:)
LOGICAL, ALLOCATABLE:: BWORK(:)

INTEGER I,J,COUNTER
DOUBLE PRECISION :: TP1
DOUBLE PRECISION ,ALLOCATABLE :: VEGRT(:),VEGIT(:),TP2(:)

INTEGER ERR

!*  DGGEV, DGGEVX computes for a pair of N-by-N real nonsymmetric matrices (A,B)
!*  the generalized eigenvalues, and optionally, the left and/or right
!*  generalized eigenvectors.           
!>   CALL DGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI, &
!>           BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
                      
!*  DGEGV, This routine is deprecated and has been replaced by routine DGGEV.


      BALANC='N' 
      JOBVL='N'
      JOBVR='N'
      SENSE='N'
      
      N=NN
      LDA=NN
      LDVL=NN
      LDVR=NN  
      LWORK=2*NN*NN+8*NN+16
     

      allocate(A(LDA, N), &
                  WR(N), WI(N), &
                  VL(N, N), VR(N, N), SCALE(N), &
                  RCONDE(N), RCONDV(N), &
                  WORK(LWORK), &
                  IWORK(N+6), &
                  VEGRT(N),VEGIT(N), &
                  TP2(N),&
                  stat=err)
	if(err.ne.0) then
		write(*, 999)
		stop
      end if
      
      A=AK 
     
      CALL  DGEEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA,&
                WR, WI, VL, LDVL, VR, LDVR, ILO,&
                IHI, SCALE, ABNRM,RCONDE, &
                RCONDV, WORK, LWORK, IWORK, INFO )


      IF( INFO .EQ. 0 ) THEN
          IMK=0 
          VEGIT(:)=WI(:)
          VEGRT(:)=WR(:)


          DO I=1,NN  
            DO J=I,NN
!               IF(VEGRT(J).LT.-1E10) THEN
!                        VEGRT(J)=-VEGRT(J)
!               END IF
               IF(VEGIT(J).LT.VEGIT(I)) THEN
                    TP1=VEGRT(I)
                    VEGRT(I)=VEGRT(J)
                    VEGRT(J)=TP1
                    
                    TP1=VEGIT(I)
                    VEGIT(I)=VEGIT(J)
                    VEGIT(J)=TP1    
                      
                END IF
            END DO
          END DO 
          
          COUNTER=0
          DO I=1,N
            IF(ABS(VEGIT(I)-0.0D0).LE.1.0D-12) THEN
                COUNTER=COUNTER+1
            END IF
        END DO
        IF(MOD(COUNTER,2).EQ.1) THEN
            WRITE(*,*) COUNTER,'ODD'
            PAUSE
        END IF
        
          DO I=1,N
            IF(VEGIT(I).GT.0.0D0) THEN
                MARK=I
                EXIT
            END IF
        END DO

!       OPEN(111,FILE='../dd.dat')
!        do i=1,n
!            WRITE(111,997)  VEGIT(i),VEGRT(i)
!        end do         
!        close(111) 
!        stop

          VEGR=VEGRT(MARK-COUNTER/2:MARK-COUNTER/2+NMD-1)
          VEGI=VEGIT(MARK-COUNTER/2:MARK-COUNTER/2+NMD-1) 
      ELSE
            IMK=1
            RETURN          
      END IF 

      deallocate(A,  &
                  WI, WR,  &
                  VL, VR, SCALE, &
                  RCONDE, RCONDV, &
                  WORK, &
                  IWORK, &
                  VEGRT,VEGIT, &
                  TP2,&
                  stat=err)
	if(err.ne.0) then
		write(*, 998)
		stop
      end if
997 format(4000f36.12)
998 FORMAT(3X,'DEALLOCATING ARRAY ERROR, PROGRAM STOPPED!')      
999 FORMAT(3X,'ALLOCATING ARRAY ERROR, PROGRAM STOPPED!')
end SUBROUTINE      
END MODULE
