    MODULE CFDCSD_CLASS

    USE INPUT_CLASS

    USE HOVER_CLASS
    USE HSQN_CLASS
    USE FORWARD_CLASS
    USE VIB_CLASS
    USE VIBH_CLASS
    USE FTRIM_CLASS
    TYPE,PUBLIC :: CFDCSD


        TYPE(INPUT) :: INPU

        TYPE(FTRIM) :: FIM
        TYPE(FTRIM),ALLOCATABLE :: FIMS(:)
        TYPE(HOVER) :: HOV
        TYPE(HSQN) :: HOVS
        TYPE(FORWARD) :: FWD
        TYPE(VIB) :: VIBS
        TYPE(VIBH) :: VIBSH

    CONTAINS

    PROCEDURE,PUBLIC :: CSDINI
    PROCEDURE,PUBLIC :: CSD
    PROCEDURE,PUBLIC :: RUNCSD
    
    END TYPE CFDCSD
    CONTAINS
    !->+++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE RUNCSD(THIS,INPUT_CFG)
    IMPLICIT NONE
    CLASS(CFDCSD) :: THIS
    CHARACTER(LEN=100),INTENT(IN) :: INPUT_CFG

    CALL THIS%CSDINI(INPUT_CFG)
    CALL THIS%CSD()

    END SUBROUTINE
    !->+++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE CSD(THIS)
    IMPLICIT NONE
    CLASS(CFDCSD) :: THIS
    INTEGER :: I,IMAX,N

    SELECT CASE(THIS%INPU%IPT%SOL)
    CASE(-1)
        THIS%INPU%IPT%SOL=0
        CALL THIS%VIBSH%SOLCASEH(THIS%INPU)
    CASE(0)
        CALL THIS%VIBS%SOLCASE(THIS%INPU)
    CASE(1)
        CALL THIS%HOV%SOLCASE(THIS%INPU)
        !CALL THIS%VIBS%SOLCASE(THIS%INPU)
        !CALL THIS%HOVS%SOLCASES(THIS%INPU)
    CASE(2)
        !DO N=1,THIS%INPU%IPT%NumofRotors
        IF(THIS%INPU%IPT%mktrim.EQ.-1) THEN
            CALL THIS%FWD%SOLCASE(THIS%INPU)
        ELSE
            !imax=21
            !allocate(this%fims(imax))
            !do i=5,imax!输出各种前进比下的结果
            !    this%inpu%ipt%mu=(i-1)*0.02
            !    call this%fims(i)%sol(this%inpu)
            !end do
            CALL THIS%FIM%SOL(THIS%INPU)
        END IF
        !END DO
    END SELECT

    END SUBROUTINE
    !->+++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE CSDINI(THIS,INPUT_CFG)
    IMPLICIT NONE
    CLASS(CFDCSD) :: THIS
    CHARACTER(LEN=100),INTENT(IN) :: INPUT_CFG
    INTEGER :: I,J,N,INF,OUF
    CHARACTER(LEN=32),ALLOCATABLE :: ARGO(:)

    
    OPEN(121,FILE=trim(adjustl(INPUT_CFG)))
    OPEN(131,FILE=trim('MAT.DAT') )
    CALL THIS%INPU%INITIAL2(131)

    CALL THIS%INPU%INITIAL1(121)
    CLOSE(121)
    CLOSE(131)
    IF(THIS%INPU%IPT%NOD.EQ.1) THEN
        CALL THIS%INPU%NODIM()
    END IF

    END SUBROUTINE
    !->+++++++++++++++++++++++++++++++++++++++++++++
    END MODULE


