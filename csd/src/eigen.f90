MODULE EIGEN_CLASS
USE BEAMH_CLASS
IMPLICIT NONE
!->+++++++++++++++++++++++++++++++++++++++++++++
TYPE,PUBLIC,EXTENDS(BEAMH) :: EIGEN

    REAL(RDT) :: OMGST

    INTEGER :: soc
    REAL(RDT) :: TFLCTRL

    INTEGER :: NMNF
    INTEGER :: NMN,NMNV,NMNW,NMNT,NMNA
    REAL(RDT), ALLOCATABLE:: VEG(:),VTEG(:,:),VTEGT(:,:)

    REAL(RDT),ALLOCATABLE :: VEGF(:), VTEGF(:,:)

    INTEGER :: NC,NQ
    CHARACTER,ALLOCATABLE :: CHARA(:) 
    INTEGER,ALLOCATABLE :: CHN(:) 
    REAL(RDT),ALLOCATABLE :: VM(:,:),WM(:,:),UM(:,:),FM(:,:),UQ(:,:)
    REAL(RDT),ALLOCATABLE :: FRE(:),FREH(:)

    CONTAINS

    PROCEDURE,PUBLIC :: CONSTRUCT_EIGEN
    PROCEDURE,PUBLIC :: MLEGGAUSS_EIGEN
    PROCEDURE,PUBLIC ::  eigensolv
    PROCEDURE,PUBLIC ::  RPOC
    PROCEDURE,PUBLIC :: DECONSTRUCT_EIGEN
    PROCEDURE,PUBLIC :: GETVTEG
    
    PROCEDURE,PRIVATE :: eigensolvn
    PROCEDURE,PRIVATE :: eigensolvg
    PROCEDURE,PUBLIC :: DECOMUQ
    PROCEDURE,PUBLIC :: COMUQ
    PROCEDURE,PRIVATE :: OUTPUTFRE

    PROCEDURE,PUBLIC :: CONSTRUCT_EIGENH
    PROCEDURE,PUBLIC :: eigensolvh
    PROCEDURE,PUBLIC ::  PROCH
    PROCEDURE,PUBLIC :: OUTPUTORF
    PROCEDURE,PUBLIC :: OUTPUTORX
    
END TYPE EIGEN
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++PUBLIC++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
CONTAINS
!->+++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE CONSTRUCT_EIGEN(THIS,INPU)
USE INPUT_CLASS
CLASS(EIGEN) :: THIS
CLASS(INPUT) :: INPU

    CALL THIS%CONSTRUCT_BEAM(INPU)

    THIS%OMGST=INPU%MAT%OMG
    THIS%NMNF=INPU%IPT%NMNF
    THIS%NMN=INPU%IPT%NMN
    THIS%NMNV=INPU%IPT%NMNV
    THIS%NMNW=INPU%IPT%NMNW
    THIS%NMNT=INPU%IPT%NMNT
    THIS%NMNA=INPU%IPT%NMNA

    THIS%NC=INPU%MAT%NBPL
    THIS%NQ=INPU%MAT%NBPL*2-1
    THIS%SOC=INPU%EGN%SOC
    THIS%TFLCTRL=INPU%EGN%TFLCTRL

    ALLOCATE(THIS%VEGF(THIS%NMNF),THIS%VTEGF(THIS%NDOFC,THIS%NMNF))

    ALLOCATE(THIS%CHARA(THIS%NMNF), THIS%CHN(THIS%NMNF),&
                    THIS%VM(THIS%NC,THIS%NMNF*2),THIS%WM(THIS%NC,THIS%NMNF*2),&
                    THIS%UM(THIS%NQ,THIS%NMNF),THIS%FM(THIS%NQ,THIS%NMNF),&
                    THIS%UQ(THIS%NDOFC,THIS%NMNF),&
                    THIS%FRE(THIS%NMNF),THIS%FREH(THIS%NMNF)) 

    ALLOCATE(THIS%VEG(THIS%NMN), THIS%VTEG(THIS%NDOFC,THIS%NMN), &
                        THIS%VTEGT(THIS%NMN,THIS%NDOFC))
END SUBROUTINE
!->+++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE CONSTRUCT_EIGENH(THIS,INPU)
USE INPUT_CLASS
CLASS(EIGEN) :: THIS
CLASS(INPUT) :: INPU

    CALL THIS%CONSTRUCT_BEAMH(INPU)

    THIS%OMGST=INPU%MAT%OMG
    THIS%NMNF=INPU%MAT%NBPL
    THIS%NC=INPU%MAT%NBPL
    THIS%NQ=INPU%MAT%NBPL*2-1
    THIS%SOC=INPU%EGN%SOC
    THIS%TFLCTRL=INPU%EGN%TFLCTRL

    ALLOCATE(THIS%VEGF(THIS%NMNF),THIS%VTEGF(THIS%NDOFCC,THIS%NMNF))

    ALLOCATE(THIS%CHARA(THIS%NMNF), THIS%CHN(THIS%NMNF),&
                    THIS%VM(THIS%NC,THIS%NMNF*2),THIS%WM(THIS%NC,THIS%NMNF*2),&
                    THIS%UM(THIS%NQ,THIS%NMNF),THIS%FM(THIS%NQ,THIS%NMNF),&
                    THIS%UQ(THIS%NDOFCC,THIS%NMNF),&
                    THIS%FRE(THIS%NMNF),THIS%FREH(THIS%NMNF)) 

    ALLOCATE(THIS%VEG(THIS%NMN), THIS%VTEG(THIS%NDOFC,THIS%NMN), &
                        THIS%VTEGT(THIS%NMN,THIS%NDOFC))    
END SUBROUTINE
!->++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE MLEGGAUSS_EIGEN(THIS,TH0,TH1C,TH1S,OMG)
CLASS(EIGEN) :: THIS
REAL(RDT),INTENT(IN) :: TH0,TH1C,TH1S,OMG
INTEGER :: SOL
REAL(RDT) :: PSI

    PSI=0.0D0
    SOL=0

    CALL THIS%UPDATE()
    CALL THIS%MLEGGAUSS_ST_BEAM(TH0,TH1C,TH1S,PSI,OMG,SOL)

END SUBROUTINE
!->+++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE RPOC(THIS) 
IMPLICIT NONE
CLASS(EIGEN) :: THIS

INTEGER :: I,J,K
REAL(RDT) :: V,W,U,F,BIG
INTEGER :: FLAG
REAL(RDT) :: VMT(THIS%NC,2*THIS%NMNF),WMT(THIS%NC,2*THIS%NMNF)
REAL(RDT) :: UMT(THIS%NQ,THIS%NMNF),FMT(THIS%NQ,THIS%NMNF)
INTEGER :: VN,WN,UN,FN

    THIS%FRE=THIS%VEGF/THIS%OMGST
    THIS%FREH=THIS%VEGF/2/PI

    THIS%UQ=THIS%VTEGF

    VN=0
    WN=0
    FN=0
    UN=0

    CALL THIS%DECOMUQ(THIS%UQ,VMT,WMT,FMT,UMT,THIS%NMNF)
    
    FMT=FMT/THIS%TFLCTRL
    
    DO I=1,THIS%NMNF
        V=ABS(VMT(1,I*2-1))
        W=ABS(WMT(1,I*2-1))
        F=ABS(FMT(1,I))
        U=ABS(UMT(1,I))
        
        BIG=-0.1
         IF(V.GT.BIG) THEN
            BIG=V
            FLAG=1
            THIS%CHARA(I)='L'
            VN=VN+1
        END IF
        IF(W.GT.BIG) THEN
            BIG=W
            FLAG=2
            THIS%CHARA(I)='F'
            WN=WN+1
            VN=VN-1
        END IF
        IF(F.GT.BIG) THEN
            BIG=F
            FN=FN+1
            IF(FLAG.EQ.1) THEN
                VN=VN-1
            ELSE IF(FLAG.EQ.2) THEN
                WN=WN-1
            END IF
            FLAG=3
            THIS%CHARA(I)='T'
        END IF   
        IF(U.GT.BIG) THEN
            BIG=U
            UN=UN+1
            IF(FLAG.EQ.1) THEN
                VN=VN-1
            ELSE IF(FLAG.EQ.2) THEN
                WN=WN-1
            ELSE IF(FLAG.EQ.3) THEN
                FN=FN-1
            END IF
            FLAG=4
            THIS%CHARA(I)='A'
        END IF
        
        IF(FLAG.EQ.1) THEN
            VMT(:,I*2-1:I*2)=VMT(:,I*2-1:I*2)/VMT(1,I*2-1)
            WMT(:,I*2-1:I*2)=0.0
            FMT(:,I)=0.0
            UMT(:,I)=0.0
            THIS%CHN(I)=VN
        ELSE IF(FLAG.EQ.2) THEN
            VMT(:,I*2-1:I*2)=0.0
            WMT(:,I*2-1:I*2)=WMT(:,I*2-1:I*2)/WMT(1,I*2-1)
            FMT(:,I)=0.0
            UMT(:,I)=0.0
            THIS%CHN(I)=WN
        ELSE IF(FLAG.EQ.3) THEN
            VMT(:,I*2-1:I*2)=0.0
            WMT(:,I*2-1:I*2)=0.0
            !FMT(:,I)=(FMT(:,I)-FMT(THIS%NQ,I))/(FMT(1,I)-FMT(THIS%NQ,I))
            FMT(:,I)=FMT(:,I)/FMT(1,I)!按诸自由度中最大幅值归一化
            UMT(:,I)=0.0
            THIS%CHN(I)=FN
        ELSE IF(FLAG.EQ.4) THEN
            VMT(:,I*2-1:I*2)=0.0
            WMT(:,I*2-1:I*2)=0.0
            FMT(:,I)=0.0
            UMT(:,I)=UMT(:,I)/UMT(1,I)
            THIS%CHN(I)=UN
        END IF
    END DO
    
        THIS%VM=VMT
        THIS%WM=WMT
        THIS%FM=FMT
        THIS%UM=UMT
    
END SUBROUTINE 
!->+++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE PROCH(THIS)
IMPLICIT NONE
CLASS(EIGEN) :: THIS

INTEGER :: I

    THIS%FRE=THIS%VEGF/THIS%OMGST
    THIS%FREH=THIS%VEGF/2/PI

    THIS%CHARA='F'
    do i=1,THIS%NMNF
        THIS%CHN(i)=i
    end do

      do 580 i=1, THIS%NMNF
            IF(THIS%CHARA(I).EQ.'F') THEN
                    THIS%WM(1:THIS%NC-1, i*2-1)=THIS%VTEGF(1:THIS%NDOFCC-1:2, i)/THIS%VTEGF(1,i)
                    THIS%WM(THIS%NC, i*2-1)=0.0D0

                    THIS%WM(1:THIS%NC-1, i*2)=THIS%VTEGF(2:THIS%NDOFCC-1:2, i)/THIS%VTEGF(1,i)
                    THIS%WM(THIS%NC, i*2)=THIS%VTEGF(THIS%NDOFCC, i)/THIS%VTEGF(1,i)
            END IF
580 continue

END SUBROUTINE
!->+++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE OUTPUTORF(THIS)
IMPLICIT NONE
CLASS(EIGEN) :: THIS
INTEGER :: I,J,K,ACC,FLAG(THIS%NMNF),TT

        OPEN(77, FILE='ORIF.DAT')
        
        FLAG=0 
        ACC=0

        DO I=1,THIS%NMNF
            IF(THIS%CHARA(I).EQ.'F') THEN
                ACC=ACC+1
                FLAG(ACC)=I
            END IF
        END DO
        
        DO J=1,THIS%NELE
            WRITE(77,997) THIS%ELEMEN(J)%NOD(1)%RLUD,&
            (THIS%WM(THIS%NC+1-J,FLAG(I)*2-1),I=1,ACC)
        END DO
        WRITE(77,997) THIS%ELEMEN(THIS%NELE)%NOD(2)%RLUD,&
            (THIS%WM(1,FLAG(I)*2-1),I=1,ACC)
 
    
      close(77)

997 format(400f26.6) 
END SUBROUTINE OUTPUTORF

!    DO J=1, 10
!        WRITE(18,997) (FRE(I,J),I=1,ROTAN),((FRE(I,J)-FREINPUT(J))/FREINPUT(J)*100,I=1,ROTAN)
!    END DO 
!->+++++++++++++++++++++++++++++++++++++++++++++
subroutine eigensolv(THIS)
implicit none
CLASS(EIGEN) :: THIS
integer :: imk

    IF(THIS%SOC.EQ.1) THEN
        CALL THIS%eigensolvn(IMK,THIS%NDOFC,THIS%UKMC,THIS%UMMC)
    ELSE IF(THIS%SOC.EQ.2) THEN
        CALL THIS%eigensolvg(IMK,THIS%NDOFC,THIS%UKMC,THIS%UMMC)
    END IF

    IF(IMK.NE.0) THEN
        CALL THIS%OUTPUTMATRIX()
        WRITE(*,*) 'eigensolv ERROR'
        STOP
    END IF

end subroutine
!->+++++++++++++++++++++++++++++++++++++++++++++
subroutine eigensolvh(THIS)
implicit none
CLASS(EIGEN) :: THIS
integer :: imk

    IF(THIS%SOC.EQ.1) THEN
        CALL THIS%eigensolvn(IMK,THIS%NDOFCC,THIS%UKMCC,THIS%UMMCC)
    ELSE IF(THIS%SOC.EQ.2) THEN
        CALL THIS%eigensolvg(IMK,THIS%NDOFCC,THIS%UKMCC,THIS%UMMCC)
    END IF

    IF(IMK.NE.0) THEN
        WRITE(*,*) 'eigensolv ERROR'
        STOP
    END IF

end subroutine
!->+++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE DECONSTRUCT_EIGEN(THIS)
IMPLICIT NONE
CLASS(EIGEN) :: THIS

    CALL THIS%DECONSTRUCT_BEAM()

    DEALLOCATE(THIS%VEGF,THIS%VTEGF)

    DEALLOCATE(THIS%CHARA, THIS%CHN,&
                    THIS%VM,THIS%WM,&
                    THIS%UM,THIS%FM,&
                    THIS%UQ,&
                    THIS%FRE,THIS%FREH) 
END SUBROUTINE
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++PRIVATE++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE DECOMUQ(THIS,UQ,VMT,WMT,FMT,UMT,NMNF)
IMPLICIT NONE
CLASS(EIGEN) :: THIS
INTEGER,INTENT(IN) :: NMNF
REAL(RDT),INTENT(IN) :: UQ(THIS%NDOFC,NMNF)
REAL(RDT),INTENT(OUT) :: VMT(THIS%NC,2*NMNF),WMT(THIS%NC,2*NMNF)
REAL(RDT),INTENT(OUT) :: UMT(THIS%NQ,NMNF),FMT(THIS%NQ,NMNF)

INTEGER :: NEXT
INTEGER :: I,J,K,ACC,BCC,II
INTEGER :: ITP,IPP
INTEGER :: COUNTER
INTEGER,ALLOCATABLE :: MARK(:)
REAL(RDT),ALLOCATABLE :: UQTP(:,:)

    ALLOCATE(MARK(THIS%NDOF),UQTP(THIS%NDOF,NMNF))
    MARK=0
    MARK(THIS%CNTS(1:THIS%CNT))=1
    COUNTER=0
    DO I=1,THIS%NDOF
        IF(MARK(I).EQ.0) THEN
            COUNTER=COUNTER+1
            UQTP(I,:)=UQ(COUNTER,:)
        ELSE
            UQTP(I,:)=0.0D0
        END IF
    END DO

    VMT=0
    WMT=0
    FMT=0
    UMT=0
    ACC=1
    BCC=1
    DO I=1,2
        DO J=1,NMNF
            VMT(ACC,I+(J-1)*2)=UQTP(I,J)
        END DO
    END DO
    
    DO I=3,4
        DO J=1,NMNF
            WMT(ACC,I-2+(J-1)*2)=UQTP(I,J)
        END DO
    END DO 
        
    DO J=1,NMNF
         FMT(BCC,J)=UQTP(5,J)
    END DO
    
    DO J=1,NMNF
         UMT(BCC,J)=UQTP(6,J)
    END DO
	
    NEXT=THIS%UVC0
    DO K=THIS%NELE,1,-1

        ACC=ACC+1
        BCC=BCC+1
  
            DO J=1,NMNF
                 FMT(BCC,J)=UQTP(NEXT+1,J)
            END DO
            DO J=1,NMNF
                 UMT(BCC,J)=UQTP(NEXT+2,J)
            END DO
            BCC=BCC+1
               
            DO I=THIS%UVC1+1,THIS%UVC1+2
                DO J=1,NMNF
                    VMT(ACC,I-THIS%UVC1+(J-1)*2)=UQTP(NEXT+I,J)
                END DO
            END DO
            DO I=THIS%UVC1+3,THIS%UVC1+4
                DO J=1,NMNF
                    WMT(ACC,I-THIS%UVC1-2+(J-1)*2)=UQTP(NEXT+I,J)
                END DO
            END DO 
            
            DO J=1,NMNF
                 FMT(BCC,J)=UQTP(NEXT+THIS%UVC1+5,J)
              END DO
            DO J=1,NMNF
                 UMT(BCC,J)=UQTP(NEXT+THIS%UVC1+6,J)
            END DO   
   
        NEXT=NEXT+THIS%UVC2    
    END DO

END SUBROUTINE
!->+++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE COMUQ(THIS,UQ,VMT,WMT,FMT,UMT,NMNF)
IMPLICIT NONE
CLASS(EIGEN) :: THIS
INTEGER,INTENT(IN) :: NMNF
REAL(RDT),INTENT(OUT) :: UQ(THIS%NDOFC,NMNF)
REAL(RDT),INTENT(IN) :: VMT(THIS%NC,2*NMNF),WMT(THIS%NC,2*NMNF)
REAL(RDT),INTENT(IN) :: UMT(THIS%NQ,NMNF),FMT(THIS%NQ,NMNF)

INTEGER :: NEXT
INTEGER :: I,J,K,ACC,BCC,II
INTEGER :: ITP,IPP
INTEGER :: COUNTER
INTEGER,ALLOCATABLE :: MARK(:)
REAL(RDT),ALLOCATABLE :: UQTP(:,:)

    ALLOCATE(MARK(THIS%NDOF),UQTP(THIS%NDOF,NMNF))

    ACC=1
    BCC=1
    DO I=1,2
        DO J=1,NMNF
            UQTP(I,J)=VMT(ACC,I+(J-1)*2)
        END DO
    END DO
    
    DO I=3,4
        DO J=1,NMNF
            UQTP(I,J)=WMT(ACC,I-2+(J-1)*2)
        END DO
    END DO 
        
    DO J=1,NMNF
         UQTP(5,J)=FMT(BCC,J)
    END DO
    
    DO J=1,NMNF
         UQTP(6,J)=UMT(BCC,J)
    END DO

    NEXT=THIS%UVC0
    DO K=THIS%NELE,1,-1
        ACC=ACC+1
        BCC=BCC+1
  
            DO J=1,NMNF
                 UQTP(NEXT+1,J)=FMT(BCC,J)
            END DO
            DO J=1,NMNF
                 UQTP(NEXT+2,J)=UMT(BCC,J)
            END DO
            BCC=BCC+1
               
            DO I=THIS%UVC1+1,THIS%UVC1+2
                DO J=1,NMNF
                    UQTP(NEXT+I,J)=VMT(ACC,I-THIS%UVC1+(J-1)*2)
                END DO
            END DO
            DO I=THIS%UVC1+3,THIS%UVC1+4
                DO J=1,NMNF
                    UQTP(NEXT+I,J)=WMT(ACC,I-THIS%UVC1-2+(J-1)*2)
                END DO
            END DO 
            
            DO J=1,NMNF
                 UQTP(NEXT+THIS%UVC1+5,J)=FMT(BCC,J)
              END DO
            DO J=1,NMNF
                 UQTP(NEXT+THIS%UVC1+6,J)=UMT(BCC,J)
            END DO   
   
        NEXT=NEXT+THIS%UVC2    
    END DO
    
    
    MARK=0
    MARK(THIS%CNTS(1:THIS%CNT))=1
    COUNTER=0
    DO I=1,THIS%NDOF
        IF(MARK(I).EQ.0) THEN
            COUNTER=COUNTER+1
            UQ(COUNTER,:)=UQTP(I,:)
        END IF
    END DO
    
END SUBROUTINE
!->+++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE OUTPUTFRE(THIS)
IMPLICIT NONE
CLASS(EIGEN) :: THIS

integer :: i,j

    OPEN(18,FILE='fre.dat')
    DO J=1, 10
        WRITE(18,995) THIS%FRE(J),THIS%CHN(J),THIS%CHARA(J)
    END DO 

      WRITE(18, *) 
      WRITE(18, *) 'Centrifugal Force:'
      DO 700 I=THIS%NELE, 1, -1
            WRITE(18, 997) THIS%ELEMEN(I)%NOD(2)%CTFGRL
700 CONTINUE 
        WRITE(18, 997) THIS%ELEMEN(1)%NOD(1)%CTFGRL
      CLOSE(18)  
      
995 format(400(f22.6,3X,I3,A1))
997 format(400f26.6)
END SUBROUTINE OUTPUTFRE	
!->+++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE OUTPUTORX(THIS)
IMPLICIT NONE
CLASS(EIGEN) :: THIS
INTEGER :: I,J,K,ACC,FLAG(THIS%NMNF),TT

        OPEN(77, FILE='ORIF.DAT')
        OPEN(78, FILE='ORIL.DAT')
        OPEN(79, FILE='ORIT.DAT')
        
        FLAG=0 
        ACC=0

        DO I=1,THIS%NMNF
            IF(THIS%CHARA(I).EQ.'F') THEN
                ACC=ACC+1
                FLAG(ACC)=I
            END IF
        END DO
        
        DO J=1,THIS%NELE
            WRITE(77,997) THIS%ELEMEN(J)%NOD(1)%RLUD,&
            (THIS%WM(THIS%NC+1-J,FLAG(I)*2-1),I=1,ACC)
        END DO
        WRITE(77,997) THIS%ELEMEN(THIS%NELE)%NOD(2)%RLUD,&
            (THIS%WM(1,FLAG(I)*2-1),I=1,ACC)
 
        FLAG=0
        ACC=0
        DO I=1,THIS%NMNF
            IF(THIS%CHARA(I).EQ.'L') THEN
                ACC=ACC+1
                FLAG(ACC)=I
            END IF
        END DO
        
        DO J=1,THIS%NELE
            WRITE(78,997) THIS%ELEMEN(J)%NOD(1)%RLUD,&
            (THIS%VM(THIS%NC+1-J,FLAG(I)*2-1),I=1,ACC)
        END DO
            WRITE(78,997) THIS%ELEMEN(THIS%NELE)%NOD(2)%RLUD,&
            (THIS%VM(1,FLAG(I)*2-1),I=1,ACC)

        FLAG=0
        ACC=0
        DO I=1,THIS%NMNF
            IF(THIS%CHARA(I).EQ.'T') THEN
                ACC=ACC+1
                FLAG(ACC)=I
            END IF
        END DO
    
        TT=0
        K=0
        DO J=1,THIS%NQ
            IF(TT.EQ.0) THEN
                K=K+1
                IF(K.GT.THIS%NELE) EXIT
                WRITE(79,997) THIS%ELEMEN(K)%NOD(1)%RLUD,&
                (THIS%FM(THIS%NQ+1-J,FLAG(I)),I=1,ACC)
                TT=1
            ELSE
                WRITE(79,997) (THIS%ELEMEN(K)%NOD(1)%RLUD+THIS%ELEMEN(K)%NOD(2)%RLUD)/2,&
                (THIS%FM(THIS%NQ+1-J,FLAG(I)),I=1,ACC)
                TT=0
            END IF
        END DO
        WRITE(79,997) THIS%ELEMEN(THIS%NELE)%NOD(2)%RLUD,&
        (THIS%FM(1,FLAG(I)),I=1,ACC)

      close(77)
      close(78) 
      close(79)



997 format(400f26.6) 
END SUBROUTINE OUTPUTORX
!->+++++++++++++++++++++++++++++++++++++++++++++
subroutine eigensolvn(THIS,IMK,NN,A,B)
implicit none
CLASS(EIGEN) :: THIS
INTEGER, INTENT(OUT):: IMK
INTEGER,INTENT(IN) :: NN
REAL(RDT),INTENT(IN) :: A(NN,NN),B(NN,NN)

INTEGER ::  NMD

INTEGER ITYPE, N, LDA, LDB, IL, IU, LDZ, LWORK 
CHARACTER JOBZ, RANGE, UPLO  
REAL(RDT) VL, VU, ABSTOL

INTEGER M, INFO
INTEGER, ALLOCATABLE:: IWORK(:), IFAIL(:)
REAL(RDT), ALLOCATABLE :: W(:), Z(:,:), WORK(:)

INTEGER ERR,i
    

    NMD=THIS%NMNF

!* DSYGV, DSYGVD, DSYGVX: LAPACK
!*  DSYGV, DSYGVD computes all the eigenvalues, and optionally, the eigenvectors
!*  of a real generalized symmetric-definite eigenproblem, of the form
!*  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and
!*  B are assumed to be symmetric and B is also positive definite.
      
      !CALL DSYGV( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK, &
      !     LWORK, INFO )         
      !CALL DSYGVD( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK, &
      !     LWORK, IWORK, LIWORK, INFO )           

      ITYPE=1
      JOBZ='V' 
      RANGE='I' 
      UPLO='U'
      N=NN
      LDA=NN
      LDB=NN
      VL=0.0D0
      VU=1000.0D0
      IL=1
      IU=NMD   
      ABSTOL=1.0D-6
      LDZ=NN 
      LWORK=35*NN 
      !>LWORK=8*NN  

      allocate(W(N), Z(LDZ, NMD), WORK(LWORK), &
                  IWORK(5*N), IFAIL(N), & 
                  stat=err)
	if(err.ne.0) then
		write(*, 999)
		stop
      end if
      
      work=0.0D0 
      iwork=0 
      CALL DSYGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB, &
            VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, &
            LWORK, IWORK, IFAIL, INFO )

      IF( INFO .EQ. 0 .AND. M .GE. 0 ) THEN
            IMK=0
            DO I=1,M
                THIS%VEGF(I)=SQRT(W(I))
            END DO
            THIS%VTEGF=Z 
      ELSE
            IMK=1
            RETURN          
      END IF 

      deallocate(W, Z, WORK, &
                  IWORK, IFAIL, & 
                  stat=err)
	if(err.ne.0) then
		write(*, 998)
		stop
      end if

998 FORMAT(3X,'DEALLOCATING ARRAY ERROR, PROGRAM STOPPED!')      
999 FORMAT(3X,'ALLOCATING ARRAY ERROR, PROGRAM STOPPED!')
end subroutine         
!->+++++++++++++++++++++++++++++++++++++++++++++
subroutine eigensolvg(THIS,IMK,NN,A,B)
implicit none
CLASS(EIGEN) :: THIS
INTEGER, INTENT(OUT):: IMK
INTEGER,INTENT(IN) :: NN
REAL(RDT),INTENT(IN) :: A(NN,NN),B(NN,NN)
!REAL(RDT), INTENT(OUT):: VEGR(NMD), VEGI(NMD), VTEG(NN, NMD)


INTEGER :: NMD
CHARACTER BALANC, JOBVL, JOBVR, SENSE
INTEGER N, LDA, LDB
INTEGER LDVL, LDVR, ILO, IHI
INTEGER LWORK
INTEGER INFO
DOUBLE PRECISION ABNRM, BBNRM
DOUBLE PRECISION ,  ALLOCATABLE:: ALPHAR(:), ALPHAI(:), BETA(:)
DOUBLE PRECISION , ALLOCATABLE:: VL(:, :), VR(:, :), LSCALE(:), RSCALE(:)
DOUBLE PRECISION , ALLOCATABLE:: RCONDE(:), RCONDV(:)
DOUBLE PRECISION , ALLOCATABLE:: WORK(:)
INTEGER, ALLOCATABLE:: IWORK(:)
LOGICAL, ALLOCATABLE:: BWORK(:)

INTEGER I,J
DOUBLE PRECISION :: TP1
DOUBLE PRECISION ,ALLOCATABLE :: VEGI(:),VEGR(:),VEGRT(:),VEGIT(:),TP2(:),VTEGT(:,:)

INTEGER ERR

!*  DGGEV, DGGEVX computes for a pair of N-by-N real nonsymmetric matrices (A,B)
!*  the generalized eigenvalues, and optionally, the left and/or right
!*  generalized eigenvectors.           
!>   CALL DGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI, &
!>           BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
                      
!*  DGEGV, This routine is deprecated and has been replaced by routine DGGEV.

    NMD=THIS%NMNF

      BALANC='N' 
      JOBVL='N'
      JOBVR='V'
      SENSE='E'
      
      N=NN
      LDA=NN
      LDB=NN
      LDVL=NN
      LDVR=NN  
      LWORK=2*NN*NN+8*NN+16
     

      allocate(ALPHAR(N), ALPHAI(N), BETA(N), &
                  VL(N, N), VR(N, N), LSCALE(N), RSCALE(N), &
                  RCONDE(N), RCONDV(N), &
                  WORK(LWORK), &
                  IWORK(N+6), &
                  BWORK(N), &
                  VEGI(N),VEGR(N),&
                  VEGRT(N),VEGIT(N), &
                  TP2(N),VTEGT(N,N),&
                  stat=err)
	if(err.ne.0) then
		write(*, 999)
		stop
      end if
      
  
      CALL DGGEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, B, LDB, &
                  ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, ILO, &
                  IHI, LSCALE, RSCALE, ABNRM, BBNRM, RCONDE, &
                  RCONDV, WORK, LWORK, IWORK, BWORK, INFO )

      IF( INFO .EQ. 0 ) THEN
          IMK=0 
          VEGIT(:)=ALPHAI(:)/BETA(:)
          VEGRT(:)=ALPHAR(:)/BETA(:)
          VTEGT=VR


          DO I=1,NN  
            DO J=I,NN
               IF(VEGRT(J).LT.0.0D0) THEN
                        VEGRT(J)=-VEGRT(J)
               END IF
               IF(VEGRT(J).LT.VEGRT(I)) THEN
                    TP1=VEGRT(I)
                    VEGRT(I)=VEGRT(J)
                    VEGRT(J)=TP1
                    
                    TP1=VEGIT(I)
                    VEGIT(I)=VEGIT(J)
                    VEGIT(J)=TP1    
                      
                    TP2=VTEGT(:,I)
                    VTEGT(:,I)=VTEGT(:,J)
                    VTEGT(:,J)=TP2  
                END IF
            END DO
          END DO 
          DO I=1,N
            IF(VEGRT(I).LT.0.0D0) THEN
                VEGRT(I)=-VEGRT(I)
                WRITE(*,*) 'CAUTION',I,VEGRT(I),'REVERSED'
            END IF
        END DO

!          VEGI=VEGIT(1:NMD)
!          VEGR=DSQRT(VEGRT(1:NMD)) 
          THIS%VEGF=DSQRT(VEGRT(1:NMD)) 
          THIS%VTEGF=VTEGT(1:NN,1:NMD)
      ELSE
            IMK=1
            RETURN          
      END IF 

      deallocate(ALPHAR, ALPHAI, BETA, &
                  VL, VR, LSCALE, RSCALE, &
                  RCONDE, RCONDV, &
                  WORK, &
                  IWORK, &
                  BWORK, &
                  VEGR,VEGI, &
                  VEGRT,VEGIT, &
                  TP2,VTEGT,&
                  stat=err)
	if(err.ne.0) then
		write(*, 998)
		stop
      end if

998 FORMAT(3X,'DEALLOCATING ARRAY ERROR, PROGRAM STOPPED!')      
999 FORMAT(3X,'ALLOCATING ARRAY ERROR, PROGRAM STOPPED!')
end SUBROUTINE         

SUBROUTINE GETVTEG(THIS,TH0,TH1C,TH1S,OMG)
IMPLICIT NONE
CLASS(EIGEN) :: THIS
REAL(RDT),INTENT(IN) :: TH0,TH1C,TH1S,OMG

INTEGER :: I,J,MODN

    CALL THIS%MLEGGAUSS_EIGEN(TH0,TH1C,TH1S,OMG)
    CALL THIS%FORMULATE_ST()
    CALL THIS%eigensolv()
    CALL THIS%RPOC()
    
    MODN=THIS%NMNF

      J=0
      do 580 i=1, MODN
            IF(THIS%CHARA(I).EQ.'F') THEN
                IF(THIS%CHN(I).LE.THIS%NMNW) THEN
                    J=J+1
                    THIS%VEG(J)=THIS%VEGF(I)
                    THIS%vteg(:, J)=THIS%vtegF(:, i)
                END IF
            END IF
580 continue

      do 590 i=1, MODN
            IF(THIS%CHARA(I).EQ.'L') THEN
                IF(THIS%CHN(I).LE.THIS%NMNV) THEN
                    J=J+1
                    THIS%VEG(J)=THIS%VEGF(I)
                    THIS%vteg(:, J)=THIS%vtegF(:, i)
                END IF
            END IF
590 continue


      do 600 i=1, MODN
            IF(THIS%CHARA(I).EQ.'T') THEN
                IF(THIS%CHN(I).LE.THIS%NMNT) THEN
                    J=J+1
                    THIS%VEG(J)=THIS%VEGF(I)
                    THIS%vteg(:, J)=THIS%vtegF(:, i)
                END IF
            END IF
600 continue

      do 610 i=1, MODN
            IF(THIS%CHARA(I).EQ.'A') THEN
                IF(THIS%CHN(I).LE.THIS%NMNA) THEN
                    J=J+1
                    THIS%VEG(J)=THIS%VEGF(I)
                    THIS%vteg(:, J)=THIS%vtegF(:, i)
                END IF
            END IF
610 continue

    THIS%VTEGT=TRANSPOSE(THIS%VTEG)

END SUBROUTINE


END MODULE
