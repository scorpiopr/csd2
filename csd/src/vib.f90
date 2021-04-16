MODULE VIB_CLASS
USE EIGEN_CLASS
IMPLICIT NONE
TYPE,PUBLIC :: VIB

    INTEGER :: ROTAN,MDOUT
    REAL(RDT) :: OMG0
    REAL(RDT),ALLOCATABLE :: OMGA(:)
    TYPE(EIGEN),ALLOCATABLE :: EIGENS(:)

CONTAINS

    PROCEDURE,PUBLIC :: SOLCASE

    PROCEDURE,PUBLIC :: CONSTRUCT_VIB
    PROCEDURE,PUBLIC :: DECONSTRUCT_VIB
    PROCEDURE,PUBLIC :: CAL_VIB
    PROCEDURE,PUBLIC :: OUTPUTFRE_VIB
    PROCEDURE,PUBLIC :: OUTPUTFRE_VIB2
    PROCEDURE,PUBLIC :: OUTPUTVIB2
END TYPE VIB
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++PUBLIC++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
CONTAINS
SUBROUTINE SOLCASE(THIS,INPU)
USE INPUT_CLASS
IMPLICIT NONE
CLASS(VIB) :: THIS
CLASS(INPUT) :: INPU

    WRITE(*,*) '+++++++++++++++++++++++++++++++++++++++++++++'
    WRITE(*,*) '构建计算模型'
    CALL THIS%CONSTRUCT_VIB(INPU)
!->++++++++++++++++++++++++++++++++++++++++++++
    WRITE(*,*) '计算共振特性'
    CALL THIS%CAL_VIB()
!->++++++++++++++++++++++++++++++++++++++++++++
    WRITE(*,*) '输出计算结果'
    CALL THIS%OUTPUTFRE_VIB2()
    CALL THIS%OUTPUTVIB2()
!->++++++++++++++++++++++++++++++++++++++++++++
    WRITE(*,*) '清理'
    CALL THIS%DECONSTRUCT_VIB()

END SUBROUTINE
!->++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE CONSTRUCT_VIB(THIS,INPU)
USE INPUT_CLASS
IMPLICIT NONE
CLASS(VIB) :: THIS
CLASS(INPUT) :: INPU

INTEGER :: I
REAL(RDT) :: TH0,TH1C,TH1S,OMG

    TH0=0.0D0
    TH1C=0.0D0
    TH1S=0.0D0
    OMG=INPU%MAT%OMG
    THIS%OMG0=OMG
    THIS%ROTAN=INPU%EGN%ROTAN
    ALLOCATE(THIS%OMGA(THIS%ROTAN),THIS%EIGENS(THIS%ROTAN))

    IF(THIS%ROTAN.EQ.1) THEN
         THIS%MDOUT=1 
        THIS%OMGA(1)=INPU%MAT%OMG
    ELSE
        DO I=1,THIS%ROTAN
            THIS%MDOUT=THIS%ROTAN/1.2
            THIS%OMGA(I)=INPU%MAT%OMG*1.2*(I-1)/(THIS%ROTAN-1)
        END DO  
    END IF

    DO I=1,THIS%ROTAN
        CALL THIS%EIGENS(I)%CONSTRUCT_EIGEN(INPU)
        CALL THIS%EIGENS(I)%MLEGGAUSS_EIGEN(TH0,TH1C,TH1S,THIS%OMGA(I))
    END DO  

END SUBROUTINE
!->+++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE CAL_VIB(THIS)
IMPLICIT NONE
CLASS(VIB) :: THIS

INTEGER :: I
INTEGER :: J,K

    DO I=1,THIS%ROTAN
		WRITE(*,*) THIS%OMGA(I)
        CALL THIS%EIGENS(I)%FORMULATE_ST()
        CALL THIS%EIGENS(I)%eigensolv()
        CALL THIS%EIGENS(I)%RPOC()
        CALL THIS%EIGENS(I)%OUTPUTORX()
    END DO  

END SUBROUTINE
!->+++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE OUTPUTFRE_VIB(THIS)
IMPLICIT NONE
CLASS(VIB) :: THIS

integer :: i,j

    OPEN(18,FILE='fre.dat')
    DO J=1, 10
        WRITE(18,995) (THIS%EIGENS(I)%FRE(J),THIS%EIGENS(I)%CHN(J),&
                                  THIS%EIGENS(I)%CHARA(J),I=1,THIS%ROTAN)
    END DO 

      WRITE(18, *) 
      WRITE(18, *) 'Centrifugal Force:'
      DO 700 J=THIS%EIGENS(1)%NELE, 1, -1
            WRITE(18, 997) (THIS%EIGENS(I)%ELEMEN(J)%NOD(2)%CTFGRL,I=1,THIS%ROTAN)
700 CONTINUE 
        WRITE(18, 997) (THIS%EIGENS(I)%ELEMEN(1)%NOD(1)%CTFGRL,I=1,THIS%ROTAN)
      CLOSE(18)  
      
995 format(400(f22.6,3X,I3,A1))
997 format(400f26.6)
END SUBROUTINE OUTPUTFRE_VIB	
!->+++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE OUTPUTFRE_VIB2(THIS)
IMPLICIT NONE
CLASS(VIB) :: THIS

integer :: i,j

    OPEN(18,FILE='fre.dat')

    WRITE(18, *) 'Frequency characteristics:'
    DO I=1,THIS%ROTAN
    WRITE(18, *) 
    WRITE(18, 993) 'OMG=',THIS%OMGA(I),'Rad/s',1.0*1.2*(I-1)/(THIS%ROTAN-1)
        DO J=1, 10
            WRITE(18,995) THIS%EIGENS(I)%FRE(J),THIS%EIGENS(I)%CHN(J),&
                                    THIS%EIGENS(I)%CHARA(J),THIS%EIGENS(I)%FREH(J),'Hz'
        END DO 

!        WRITE(18, *) 
!        WRITE(18, *) 'Centrifugal Force:'
!        DO J=THIS%EIGENS(1)%NELE, 1, -1
!            WRITE(18, 997) THIS%EIGENS(I)%ELEMEN(J)%NOD(2)%CTFGRL
!        END DO
!        WRITE(18, 997) THIS%EIGENS(I)%ELEMEN(1)%NOD(1)%CTFGRL
    END DO

    CLOSE(18)  
993 format(A4,f12.6,A10,f12.6)
995 format(400(f22.6,3X,I3,A1,f22.6,A2))
997 format(400f26.6)
END SUBROUTINE OUTPUTFRE_VIB2	
!->+++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE OUTPUTVIB2(THIS)
IMPLICIT NONE
CLASS(VIB) :: THIS
	
INTEGER :: I,J,K
    
    OPEN(77, FILE='VIB.DAT')
    DO I=1,THIS%ROTAN
        WRITE(77,997) THIS%OMGA(I)/THIS%OMG0,(THIS%EIGENS(I)%FRE(J),J=1,10),&
                                (THIS%OMGA(I)/THIS%OMGA(THIS%MDOUT)*J,J=1,8)
    END DO

    CLOSE(77)

997 format(400f26.6) 
END SUBROUTINE OUTPUTVIB2	
!->+++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE DECONSTRUCT_VIB(THIS)
IMPLICIT NONE
CLASS(VIB) :: THIS

    DEALLOCATE(THIS%OMGA,THIS%EIGENS)

END SUBROUTINE
END MODULE
