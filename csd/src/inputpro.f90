    MODULE INPUTPRO_CLASS
    USE GlobalDataFun
    IMPLICIT NONE
    !->+++++++++++++++++++++++++++++++++++++++++++++
    TYPE,PUBLIC :: INPUTPRO

        INTEGER :: CHO
        REAL(RDT) :: CSM(4,4),TSM(6,6),TMM(6,6)
        REAL(RDT) :: BTA,TAO,RL

    CONTAINS
    PROCEDURE,PUBLIC :: INITIAL_MAT
    PROCEDURE,PUBLIC :: SET=>SET_V
    PROCEDURE,PUBLIC :: OUT=>OUT_V
    PROCEDURE,PUBLIC :: SET_M
    PROCEDURE,PUBLIC :: SET_S
    PROCEDURE,PUBLIC :: GETMAT

    END TYPE INPUTPRO
    !->+++++++++++++++++++++++++++++++++++++++++++++
    CONTAINS
    !->+++++++++++++++++++++++++++++++++++++++++++++
    !->+++++++++++++++++++++++++++++++++++++++++++++
    !->+++++++++++++++PUBLIC++++++++++++++++++++++++++
    !->+++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE INITIAL_MAT(THIS,IIN,IOUT,INPU)
    USE INPUT_CLASS
    IMPLICIT NONE
    CLASS(INPUTPRO) :: THIS
    INTEGER,INTENT(IN) :: IIN,IOUT
    CLASS(INPUT) :: INPU

    !INPU%MOD%CHO
    !-1 为欧拉梁，读取文件同控制参数，如果KM1和KM2大于0，则采用修改过的惯性矩
    !0和1，采用vabs格式，读取新文件，如果KM1和KM2大于0，则采用修改过的惯性矩
    !3和4，采用vabs格式，读取新文件，完全不修改惯性矩
    !5，为欧拉梁，读取外部文件
    IF(INPU%MOD%CHO.EQ.-1) THEN
        THIS%CHO=0
        CALL THIS%SET_M(INPU)
        CALL THIS%SET_S(INPU)
        IF(INPU%MOD%KM1.GE.0.0D0) THEN
            THIS%TMM(5,5)=(INPU%MOD%KM1*INPU%IPT%R)**2
        END IF
        IF(INPU%MOD%KM2.GE.0.0D0) THEN
            THIS%TMM(6,6)=(INPU%MOD%KM2*INPU%IPT%R)**2
        END IF
    ELSE IF(INPU%MOD%CHO.EQ.0.OR.INPU%MOD%CHO.EQ.1) THEN
        THIS%CHO=INPU%MOD%CHO
        CALL THIS%SET(IIN)
        CALL THIS%SET_M(INPU)
        IF(INPU%MOD%KM1.GE.0.0D0) THEN
            THIS%TMM(5,5)=(INPU%MOD%KM1*INPU%IPT%R)**2
        END IF
        IF(INPU%MOD%KM2.GE.0.0D0) THEN
            THIS%TMM(6,6)=(INPU%MOD%KM2*INPU%IPT%R)**2
        END IF
    ELSE IF(INPU%MOD%CHO.EQ.3) THEN
        THIS%CHO=0
        CALL THIS%SET(IIN)
    ELSE IF(INPU%MOD%CHO.EQ.4) THEN
        THIS%CHO=1
        CALL THIS%SET(IIN)
    ELSE IF(INPU%MOD%CHO.EQ.5) THEN
        THIS%CHO=0
        CALL THIS%GETMAT(IIN)
    END IF

    CALL THIS%OUT(IOUT)

    END SUBROUTINE
    !->+++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE SET_V(THIS,IIN)
    IMPLICIT NONE
    CLASS(INPUTPRO) :: THIS
    INTEGER,INTENT(IN) :: IIN

    INTEGER :: ISTAT,I,J
    character( len = 10000 ) :: str

    REWIND(IIN)
    DO
        read(IIN,'(a)',IOSTAT=ISTAT) str
        IF(ISTAT.GT.0) THEN
            WRITE(*,*) 'INPUT FILE ERROR'
            STOP
        ELSE IF(ISTAT.LT.0) THEN
            EXIT
        END IF

        IF(trim(adjustl(STR)).EQ.'Classical Stiffness Matrix (1-extension; 2-twist; 3,4-bending)') THEN
            read(IIN,'(a)') str
            DO I=1,4
                READ(IIN,*) (THIS%CSM(I,J),J=1,4)
            END DO
        END IF

        IF(trim(adjustl(STR)).EQ.'Timoshenko Stiffness Matrix (1-extension; 2,3-shear, 4-twist; 5,6-bending)') THEN
            read(IIN,'(a)') str
            DO I=1,6
                READ(IIN,*) (THIS%TSM(I,J),J=1,6)
            END DO
        END IF

        IF(trim(adjustl(STR)).EQ.'The 6X6 Mass Matrix') THEN
            read(IIN,'(a)') str
            DO I=1,6
                READ(IIN,*) (THIS%TMM(I,J),J=1,6)
            END DO
        END IF

    END DO

    END SUBROUTINE
    !->+++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE GETMAT(THIS,IIN)
    IMPLICIT NONE
    CLASS(INPUTPRO) :: THIS
    INTEGER,INTENT(IN) :: IIN

    REAL(RDT) ::RI,MASS,ZZT,EET,XC,XI,KP2,EIXXIN,EIZZIN,GJIN,TWISTI,TA0I,BSCT


    read(IIN,*) RI,MASS,ZZT,EET,XI,&
        KP2,EIZZIN,EIXXIN,GJIN,XC,&
        TWISTI,TA0I,BSCT

    THIS%TMM=0.0D0
    THIS%CSM=0.0D0

    THIS%TMM(1,1)=MASS
    THIS%TMM(2,2)=MASS
    THIS%TMM(3,3)=MASS

    THIS%TMM(5,5)=EET
    THIS%TMM(6,6)=ZZT
    THIS%TMM(4,4)=ZZT+EET

    THIS%TMM(1,6)=-XI*MASS
    THIS%TMM(3,4)=XI*MASS
    THIS%TMM(4,3)=XI*MASS
    THIS%TMM(6,1)=-XI*MASS

    THIS%CSM=0.0D0
    THIS%CSM(3,3)=EIXXIN
    THIS%CSM(4,4)=EIZZIN
    THIS%CSM(2,2)=GJIN
    THIS%CSM(1,1)=(THIS%CSM(3,3)+THIS%CSM(4,4))/KP2
    THIS%CSM(1,4)=-THIS%CSM(1,1)*XC

    THIS%BTA=TWISTI*PI/180
    THIS%RL=RI
    THIS%TAO=TA0I*PI/180

    END SUBROUTINE
    !->+++++++++++++++++++++++++++++++++++++++++++++
    !输出截面属性
    !->+++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE OUT_V(THIS,IOUT)
    IMPLICIT NONE
    CLASS(INPUTPRO) :: THIS
    INTEGER,INTENT(IN) :: IOUT

    INTEGER :: I,J

    WRITE(IOUT,*) THIS%CHO
    WRITE(IOUT,*) THIS%RL
    WRITE(IOUT,*) THIS%BTA
    WRITE(IOUT,*) THIS%TAO

    DO I=1,6
        WRITE(IOUT,997) (THIS%TMM(I,J),J=1,6)
    END DO

    IF(THIS%CHO.EQ.0) THEN
        DO I=1,4
            WRITE(IOUT,997) (THIS%CSM(I,J),J=1,4)
        END DO
    ELSE
        DO I=1,6
            WRITE(IOUT,997) (THIS%TSM(I,J),J=1,6)
        END DO
    END IF

997 FORMAT(100(E40.20,2x))
    END SUBROUTINE

    !->+++++++++++++++++++++++++++++++++++++++++++++
    !设置质量矩阵
    !->+++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE SET_M(THIS,INPU)
    USE INPUT_CLASS
    IMPLICIT NONE
    CLASS(INPUTPRO) :: THIS
    CLASS(INPUT) :: INPU
    
    INTEGER :: I
    REAL(RDT) :: MASS

    MASS=INPU%GETMASS()
    THIS%TMM=0.0
    THIS%TMM(1,1)=MASS
    THIS%TMM(2,2)=MASS
    THIS%TMM(3,3)=MASS

    END SUBROUTINE
    !->+++++++++++++++++++++++++++++++++++++++++++++
    !设置刚度矩阵
    !->+++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE SET_S(THIS,INPU)
    USE INPUT_CLASS
    IMPLICIT NONE
    CLASS(INPUTPRO) :: THIS
    CLASS(INPUT) :: INPU
    
    INTEGER :: I
    REAL(RDT) :: MASS

    MASS=INPU%GETMASS()
    THIS%CSM=0.0
    THIS%CSM(3,3)=MASS*INPU%MOD%EIEE*INPU%IPT%OMG**2*INPU%IPT%R**4
    THIS%CSM(4,4)=MASS*INPU%MOD%EIZZ*INPU%IPT%OMG**2*INPU%IPT%R**4
    THIS%CSM(2,2)=MASS*INPU%MOD%GJ*INPU%IPT%OMG**2*INPU%IPT%R**4
    THIS%CSM(1,1)=(THIS%CSM(3,3)+THIS%CSM(4,4))/&
        (INPU%MOD%KP2*&
        ((INPU%MOD%KM1*INPU%IPT%R)**2+&
        (INPU%MOD%KM2*INPU%IPT%R)**2))
    
    END SUBROUTINE
    END MODULE
