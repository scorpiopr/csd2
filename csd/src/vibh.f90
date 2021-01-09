MODULE VIBH_CLASS
USE VIB_CLASS
IMPLICIT NONE
TYPE,PUBLIC,EXTENDS(VIB) :: VIBH


CONTAINS

    PROCEDURE,PUBLIC :: SOLCASEH

    PROCEDURE,PRIVATE :: CONSTRUCT_VIBH
    PROCEDURE,PRIVATE :: CAL_VIBH

END TYPE VIBH
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++PUBLIC++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
CONTAINS
SUBROUTINE SOLCASEH(THIS,INPU)
USE INPUT_CLASS
IMPLICIT NONE
CLASS(VIBH) :: THIS
CLASS(INPUT) :: INPU

    WRITE(*,*) '+++++++++++++++++++++++++++++++++++++++++++++'
    WRITE(*,*) '构建计算模型'
    CALL THIS%CONSTRUCT_VIBH(INPU)
!->++++++++++++++++++++++++++++++++++++++++++++
    WRITE(*,*) '计算共振特性'
    CALL THIS%CAL_VIBH()
!->++++++++++++++++++++++++++++++++++++++++++++
    WRITE(*,*) '输出计算结果'
    CALL THIS%OUTPUTFRE_VIB2()
    CALL THIS%OUTPUTVIB2()
!->++++++++++++++++++++++++++++++++++++++++++++
    WRITE(*,*) '清理'
    CALL THIS%DECONSTRUCT_VIB()

END SUBROUTINE
!->++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE CONSTRUCT_VIBH(THIS,INPU)
USE INPUT_CLASS
IMPLICIT NONE
CLASS(VIBH) :: THIS
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
        CALL THIS%EIGENS(I)%CONSTRUCT_EIGENH(INPU)
        CALL THIS%EIGENS(I)%MLEGGAUSS_EIGEN(TH0,TH1C,TH1S,THIS%OMGA(I))
    END DO  

END SUBROUTINE
!->+++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE CAL_VIBH(THIS)
IMPLICIT NONE
CLASS(VIBH) :: THIS

INTEGER :: I

    DO I=1,THIS%ROTAN
        CALL THIS%EIGENS(I)%FORMULATE_ST()
        CALL THIS%EIGENS(I)%GIVECONSH()
        CALL THIS%EIGENS(I)%eigensolvh()
        CALL THIS%EIGENS(I)%PROCH()
        CALL THIS%EIGENS(I)%OUTPUTORF()
    END DO  

END SUBROUTINE
END MODULE
