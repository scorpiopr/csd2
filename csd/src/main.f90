    PROGRAM main
    USE CFDCSD_CLASS
    IMPLICIT NONE
    TYPE(CFDCSD) :: FSI
    INTEGER :: N
    CHARACTER(LEN=100) :: INPUT_CFG
    !CHARACTER(LEN=100),ALLOCATABLE :: ARGO(:)

    !CALL GETCOMAND(N,ARGO)
    
    !INPUT_CFG=ARGO(1)
    
    CALL GETCFGFILE(INPUT_CFG)
    
    CALL FSI%RUNCSD(INPUT_CFG)

    !->+++++++++++++++++++++++++++++++++++++++++++++
    CONTAINS
    SUBROUTINE GETCOMAND(N,ARGO)
    IMPLICIT NONE
    INTEGER,INTENT(OUT) :: N
    CHARACTER(len=100),DIMENSION(:),ALLOCATABLE,INTENT(OUT) :: ARGO

    INTEGER :: i
    CHARACTER(len=100) :: arg

    N=COMMAND_ARGUMENT_COUNT()

    IF(N.EQ.0) THEN
        WRITE(*,'(/,a)') ' Please provide an input file name, executing as ./cfdcsd input_file_name'
        STOP
    END IF

    ALLOCATE(ARGO(N))

    i = 1
    DO
        CALL get_command_argument(i, arg)
        IF (LEN_TRIM(arg) == 0) EXIT
        ARGO(I)=TRIM(arg)
        i = i+1
    END DO

    END SUBROUTINE
    !END PROGRAM
    
    !->+++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE GETCFGFILE(INPUT_CFG)
    IMPLICIT NONE
    CHARACTER(LEN=100) :: INPUT_CFG
    CHARACTER(LEN=100) :: INPUTCFG

    
    WRITE(*,*),'输入配置文件名：'
    READ(*,*),INPUTCFG
    
    INPUT_CFG=TRIM(ADJUSTL(INPUTCFG))
    
    END SUBROUTINE
    
    END PROGRAM
