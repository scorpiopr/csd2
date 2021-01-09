    MODULE HOVER_CLASS
    USE EIGENBEAM_CLASS
    USE AERO_CLASS
    USE HTRIM_CLASS
    IMPLICIT NONE
    !->+++++++++++++++++++++++++++++++++++++++++++++
    TYPE,PUBLIC :: HOVER

        TYPE(EIGENBEAM) :: BLADE
        TYPE(AERO) :: AEROF
        TYPE(HTRIM) :: MTRIM

        REAL(RDT) :: TH0,OMG
        INTEGER :: SOL=1
        REAL(RDT) :: TH1S=0.0D0,TH1C=0.0D0,PSI=0.0D0
    CONTAINS

    PROCEDURE,PUBLIC :: SOLCASE
    PROCEDURE,PUBLIC :: CONSTRUCT_HOVER
    PROCEDURE,PUBLIC ::  hoverdriver
    PROCEDURE,PUBLIC ::  HOVER2
    PROCEDURE,PUBLIC :: HOVERTRIM
    PROCEDURE,PUBLIC :: OUTPUT_HOVER
    PROCEDURE,PUBLIC :: GET_DEFORM
    PROCEDURE,PUBLIC :: GET_RLUDAUD
    PROCEDURE,PRIVATE :: hoverini

    END TYPE HOVER
    !->+++++++++++++++++++++++++++++++++++++++++++++
    !->+++++++++++++++++++++++++++++++++++++++++++++
    !->+++++++++++++++PUBLIC++++++++++++++++++++++++++
    !->+++++++++++++++++++++++++++++++++++++++++++++
    !->+++++++++++++++++++++++++++++++++++++++++++++
    CONTAINS
    !->+++++++++++++++++++++++++++++++++++++++++++++
    subroutine SOLCASE(THIS,INPU)
    implicit none
    CLASS(HOVER) :: THIS
    TYPE(INPUT) :: INPU
    REAL(RDT) :: DROP
    INTEGER :: I,J,K

    WRITE(*,*) '+++++++++++++++++++++++++++++++++++++++++++++'
    WRITE(*,*) '构建计算模型'
    CALL THIS%CONSTRUCT_HOVER(INPU)
    !->++++++++++++++++++++++++++++++++++++++++++++
    WRITE(*,*) '计算悬停特性'
    IF(THIS%AEROF%MKFM .EQ. 6) THEN
    OPEN(444,FILE='cl1_rotor.dat')
    DO I=1,THIS%AEROF%NSTP
        DO J=1,NPTS
            READ(444,"(1000F)") DROP,DROP,(THIS%AEROF%CFDFPSI(K,J,I),K=1,3),&
								(THIS%AEROF%CFDMPSI(K,J,I),K=1,3)
                !IF(RLUDAUD(J).LT.0.5) THEN		
				!	THIS%AEROF%CFDFPSI(:,J,I)=0.0
				!	THIS%AEROF%CFDMPSI(:,J,I)=0.0			
				!END IF			
				!				
				!IF(DABS(THIS%AEROF%CFDFPSI(2,J,I)).GT.0.1) THEN
				!	THIS%AEROF%CFDFPSI(1:3,J,I)=0.0
				!END IF
				!IF(DABS(THIS%AEROF%CFDMPSI(1,J,I)).GT.0.1) THEN
				!	THIS%AEROF%CFDMPSI(1:3,J,I)=0.0
				!END IF					
								
								
        END DO
    END DO
    close(444)
    THIS%AEROF%MKFM=6
    END IF
    CALL THIS%hoverdriver()
    !->++++++++++++++++++++++++++++++++++++++++++++
    WRITE(*,*) '输出悬停计算结果'
    CALL THIS%OUTPUT_HOVER()
    CALL THIS%GET_DEFORM(0,INPU%IPT%NumofRotors)
    CALL THIS%GET_RLUDAUD()!输出展向气动点相对坐标
    end subroutine SOLCASE
    !->+++++++++++++++++++++++++++++++++++++++++++++
    subroutine CONSTRUCT_HOVER(THIS,INPU)
    implicit none
    CLASS(HOVER) :: THIS
    TYPE(INPUT) :: INPU

    THIS%TH0=INPU%IPT%TH0
    THIS%OMG=INPU%MAT%OMG

    CALL THIS%MTRIM%CONSTRUCT_HTRIM(INPU)
    CALL THIS%BLADE%CONSTRUCT_EIGENBEAM(THIS%TH0,THIS%TH1C,THIS%TH1S,THIS%OMG,INPU)
    CALL THIS%AEROF%CONSTRUCT_AERO(INPU,THIS%OMG,THIS%BLADE%R)
    !CALL THIS%AEROF%CONSTRUCT_AERO2(THIS%AEROF%NSTP,THIS%OMG,INPU)!计算PSIV

    CALL THIS%hoverini()

    end subroutine CONSTRUCT_HOVER
    !->+++++++++++++++++++++++++++++++++++++++++++++
    subroutine hoverini(THIS)
    implicit none
    CLASS(HOVER) :: THIS

    CALL THIS%AEROF%UNIFLOW()

    end subroutine hoverini
    !->+++++++++++++++++++++++++++++++++++++++++++++
    !-> the purpose:
    subroutine hoverdriver(THIS)
    implicit none
    CLASS(HOVER) :: THIS
    !->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !-> dummy arguments / interface:
    !->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !-> FUNCTIONS INVOKED:
    !->
    real(rdt), external:: LADIDHV, PHES75
    !->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !-> SUBROUTINES INVOKED:
    !->
    !->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !-> LOCAL VARIABLES:
    integer, parameter:: nsthv=1000
    integer i,j,k
    INTEGER :: FLAG
    !->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !-> main body:

    IF(THIS%MTRIM%MKTRIM.NE.-1) THEN
        CALL THIS%AEROF%GETTH0(THIS%TH0)
    END IF

    FLAG=0
    do i=1,nsthv
        CALL THIS%HOVER2()
        CALL THIS%HOVERTRIM(FLAG,THIS%MTRIM%MKTRIM)
        IF(FLAG.EQ.1) EXIT
    END DO

993 format(1000f26.12)
    end subroutine hoverdriver
    !->+++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE HOVER2(THIS)
    implicit none
    CLASS(HOVER) :: THIS
    REAL(RDT), EXTERNAL:: DNRM2
    character(80) filenam
    integer i, ii,IIRD
    REAL(RDT) NUQ
    REAL(RDT), ALLOCATABLE:: UQS(:),UQSX(:)
    REAL(RDT) :: epss=1.0D-11

    ALLOCATE(UQS(THIS%BLADE%NDOFC),UQSX(THIS%BLADE%NDOFC))

    IIRD=0
    DO II=1,100
        IIRD=IIRD+1
        CALL THIS%BLADE%UPDATE()
        CALL THIS%AEROF%MLEGGAUSS_AERO(THIS%OMG,THIS%BLADE,THIS%TH0,THIS%TH1C,THIS%TH1S,THIS%PSI,THIS%SOL)
        CALL THIS%BLADE%MLEGGAUSS_ST_BEAM(THIS%TH0,THIS%TH1C,THIS%TH1S,THIS%PSI,THIS%OMG,THIS%SOL)
        CALL THIS%BLADE%FORMULATE_ST()

        UQSX=THIS%BLADE%UQS
        call THIS%BLADE%GETSTATIC()
        UQS=THIS%BLADE%UQS
        NUQ=DNRM2(THIS%BLADE%NDOFC,UQS-UQSX,1)
        WRITE(*,*) 'ITER:',II,'  RES:',NUQ

        if( nuq .lt. epss) then
            exit
        end if
    END DO
    DEALLOCATE(UQS,UQSX)

997 format(1000f40.12)
    END SUBROUTINE HOVER2
    !->+++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE HOVERTRIM(THIS,FLAG,MKTRIM)
    IMPLICIT NONE
    CLASS(HOVER) :: THIS
    INTEGER,INTENT(OUT) :: FLAG
    INTEGER,INTENT(IN) :: MKTRIM

    REAL(RDT) :: DTH0
    real(rdt), parameter:: epshv=1.0d-6
    REAL(RDT) :: TH75,CTT

    TRIM:SELECT CASE(MKTRIM)
    CASE(-1)
        DTH0=0.0D0
        goto 99
    CASE(1)
        !TRIM 1
        DTH0=THIS%AEROF%GETTH02()-THIS%TH0

        IF( ABS(DTH0) .LT. epshv/100.0 ) THEN
            goto 99
        END IF
        THIS%TH0=THIS%TH0+DTH0
        CALL THIS%AEROF%UNIFLOW()
        !TRIM 2
    CASE(2)
        TH75=THIS%TH0+THIS%BLADE%PHES75()
        DTH0=THIS%AEROF%GETTH02()-TH75
        IF( ABS(DTH0) .LT. epshv/100.0 ) THEN
            goto 99
        END IF
        THIS%TH0=THIS%TH0+DTH0
        CALL THIS%AEROF%UNIFLOW()
        !TRIM 3
    CASE(3)
        TH75=THIS%TH0+THIS%BLADE%PHES75()
        DTH0=THIS%AEROF%GETTH02()-TH75
        IF( ABS(DTH0) .LT. epshv/100.0 ) THEN
            goto 99
        END IF
        THIS%TH0=THIS%TH0+DTH0
        TH75=THIS%TH0+THIS%BLADE%PHES75()
        CALL THIS%AEROF%LADIDHV(TH75)
        !TRIM 4
    CASE(4)
        TH75=THIS%TH0+THIS%BLADE%PHES75()
        CTT=THIS%AEROF%GETCT(TH75)
        DTH0=(THIS%AEROF%CT-CTT)*6.0/THIS%AEROF%alcsrf/THIS%AEROF%sgslid
        IF( ABS(DTH0) .LT. epshv/100.0 ) THEN
            goto 99
        END IF
        THIS%TH0=THIS%TH0+DTH0
        CALL THIS%AEROF%LADIDHV(TH75)
        !TRIM 5
    CASE(5)
        CALL THIS%BLADE%MLEGGAUSS_FR_BEAM(THIS%TH0,THIS%TH1C,THIS%TH1S,THIS%PSI,THIS%OMG,THIS%SOL)
        CALL THIS%BLADE%FORMULATE_FR()
        CTT=THIS%AEROF%nrbd*THIS%BLADE%BRFC(3)/THIS%AEROF%trf
        DTH0=(THIS%AEROF%CT-CTT)*6.0/THIS%AEROF%alcsrf/THIS%AEROF%sgslid
        IF( ABS(DTH0) .LT. epshv ) THEN
            goto 99
        END IF
        THIS%TH0=THIS%TH0+DTH0
        TH75=THIS%TH0+THIS%BLADE%PHES75()
        CALL THIS%AEROF%LADIDHV(TH75)
    END SELECT TRIM

    WRITE(*, 992) DTH0,THIS%TH0/PI*180
    RETURN

99  FLAG=1

992 format(100f20.10)
    END SUBROUTINE
    !->+++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE OUTPUT_HOVER(THIS)
    IMPLICIT NONE
    CLASS(HOVER) :: THIS

    open(33, file='howvf.dat',position='append')
    WRITE(33,995) '总距角(deg):',THIS%th0/PI*180,&
        '翼尖挥舞变形:',(THIS%BLADE%UQS(1)*SIN(THIS%th0)+THIS%BLADE%UQS(3)*COS(THIS%th0))/THIS%BLADE%R,&
        '翼尖摆振变形:',(THIS%BLADE%UQS(1)*COS(THIS%th0)-THIS%BLADE%UQS(3)*SIN(THIS%th0))/THIS%BLADE%R,&
        '翼尖扭转变形:',THIS%BLADE%UQS(5),&
        'CT:',THIS%AEROF%GETCT(THIS%th0),&
        'CT:',THIS%AEROF%CT
    CLOSE(33)


995 format((A13,f26.12))
    END SUBROUTINE
    !->+++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE GET_DEFORM(THIS,N,NumofRotors)
    IMPLICIT NONE
    CLASS(HOVER) :: THIS
    integer :: I,j,k,N,NumofRotors,NumofBladesEachRotor
    REAL(RDT):: RXUD
    real(rdt) :: REA(3), TED(3, 3), TRE(3, 3)
    
    RXUD=19.965
    NumofBladesEachRotor=THIS%AEROF%nrbd
    
    OPEN(32, FILE='predeform.dat')
    open(33, file='deform.dat')
    
    WRITE(32,*) NumofRotors
    WRITE(32,*) NumofBladesEachRotor
    !WRITE(32,*) THIS%nstp!一圈站位数
    WRITE(32,*) 1
    WRITE(32,*) 101
    WRITE(32,*) THIS%OMG
    WRITE(32,*) THIS%BLADE%R
    WRITE(32,*) THIS%BLADE%EL1
    
    do i=1,101
        RXUD=(RLUD(NNOD)-RLUD(1))*(i-1)/100+RLUD(1)!考虑根切
        call THIS%BLADE%rvndlct2(REA,TED,TRE,RXUD,THIS%OMG,0)
        WRITE(33,995) (RXUD-RLUD(1))/(RLUD(NNOD)-RLUD(1))!展向相对位置
        WRITE(33,995) REA(1:3)/THIS%BLADE%R
        WRITE(33,995) (TED(1:3,j),j=1,3)
        WRITE(33,995) (TRE(1:3,j),j=1,3)!肖宇论文Tpb
        WRITE(33,995) 
    end do
    
    CLOSE(32)
    CLOSE(33)


995 format(3(f26.12))

    END SUBROUTINE
    !->+++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE GET_RLUDAUD(THIS)
    IMPLICIT NONE
    CLASS(HOVER) :: THIS
    integer :: I
    
    OPEN(33, FILE='RLUDAUD.dat')
    DO I=1,NPTS
        WRITE(33, 995) RLUDAUD(I)
    END DO
    CLOSE(33)
    
995 format(400f26.12)
    
    END SUBROUTINE
    
    END MODULE
