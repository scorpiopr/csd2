MODULE HSQN_CLASS
USE HOVER_CLASS
USE STAB_CLASS
TYPE,PUBLIC :: HSQN

    INTEGER :: NHOV
    TYPE(HOVER),ALLOCATABLE :: HOVS(:)
    TYPE(STABSIS),ALLOCATABLE :: STAB(:)
    CONTAINS

    PROCEDURE,PUBLIC :: CONSTRUCT_HSQN
    PROCEDURE,PUBLIC :: SOLCASES
    PROCEDURE,PUBLIC :: CONSTRUCT_STAB

END TYPE
CONTAINS
!->+++++++++++++++++++++++++++++++++++++++++++++
subroutine CONSTRUCT_HSQN(THIS,INPU)    
IMPLICIT NONE
CLASS(HSQN) :: THIS
TYPE(INPUT) :: INPU
INTEGER :: I
REAL(RDT) :: TH75

    THIS%NHOV=50
    ALLOCATE(THIS%HOVS(THIS%NHOV))

    DO I=1,THIS%NHOV
        
!        INPU%IPT%LDTIP(3)=-(I-1)*10.0D0
        CALL THIS%HOVS(I)%CONSTRUCT_HOVER(INPU)    

        THIS%HOVS(I)%TH0=(i-20.0)*0.01
        TH75=THIS%HOVS(I)%TH0+THIS%HOVS(I)%BLADE%PHES75()
        CALL THIS%HOVS(I)%AEROF%LADIDHV(TH75) 

!        THIS%HOVS(I)%AEROF%CT=(i)*1.0*0.001/2
!        THIS%HOVS(I)%AEROF%CTC=(i-1)*1.0*0.001
!        THIS%HOVS(I)%AEROF%CT=0.01
!        THIS%HOVS(I)%AEROF%CTC=0.02
!        CALL THIS%HOVS(I)%AEROF%UNIFLOW()
    END DO

    open(33, file='howvf.dat',position='append')
    CLOSE (33, STATUS = 'DELETE')

END SUBROUTINE
!->+++++++++++++++++++++++++++++++++++++++++++++
subroutine SOLCASES(THIS,INPU)    
implicit none
CLASS(HSQN) :: THIS
TYPE(INPUT) :: INPU

INTEGER :: I

   WRITE(*,*) '+++++++++++++++++++++++++++++++++++++++++++++'
    WRITE(*,*) '构建计算模型'
    CALL THIS%CONSTRUCT_HSQN(INPU)
!->++++++++++++++++++++++++++++++++++++++++++++
    WRITE(*,*) '计算悬停特性'
    DO I=1,THIS%NHOV
        WRITE(*,*) 'TH0:',THIS%HOVS(I)%TH0
        CALL THIS%HOVS(I)%hoverdriver()
        WRITE(*,*) 
    END DO
!->++++++++++++++++++++++++++++++++++++++++++++
    WRITE(*,*) '稳定性分析'
    CALL THIS%CONSTRUCT_STAB(INPU)
    DO I=1,THIS%NHOV
        CALL THIS%STAB(I)%getmatrixam(THIS%HOVS(I))
        CALL THIS%STAB(I)%SOLVSm()
    END DO
!->++++++++++++++++++++++++++++++++++++++++++++
    WRITE(*,*) '输出悬停计算结果'
    DO I=1,THIS%NHOV
        CALL THIS%HOVS(I)%OUTPUT_HOVER()
        CALL THIS%STAB(I)%OUTPUT_STAB(THIS%HOVS(I))
    END DO
    
end subroutine SOLCASES 
!->+++++++++++++++++++++++++++++++++++++++++++++
subroutine CONSTRUCT_STAB(THIS,INPU)    
IMPLICIT NONE
CLASS(HSQN) :: THIS
CLASS(INPUT) :: INPU
INTEGER :: i

    ALLOCATE(THIS%STAB(THIS%NHOV))
    DO I=1,THIS%NHOV
        CALL THIS%STAB(I)%CONSTRUCT_HOVSTAB(THIS%HOVS(I),INPU) 
    END DO

END SUBROUTINE
!->+++++++++++++++++++++++++++++++++++++++++++++
END MODULE
