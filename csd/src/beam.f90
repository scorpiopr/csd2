    MODULE BEAM_CLASS
    USE ELEMENT_CLASS
    IMPLICIT NONE
    !->+++++++++++++++++++++++++++++++++++++++++++++
    TYPE,PUBLIC :: BEAM

        INTEGER :: SMD,NELE,NDOF,NDOFC
        INTEGER :: UVC0,UVC1,UVC2

        INTEGER :: BEAMTYPE

        INTEGER :: CNT
        INTEGER,ALLOCATABLE :: CNTS(:)

        REAL(RDT) :: R,EL1,BTAP,BTAJ
        REAL(RDT) :: THP,THP1D,THP2D

        REAL(RDT) :: LDTIP(3)
        INTEGER :: npt,nfp,nll

        REAL(RDT) :: KPTBAR,KFP
        REAL(RDT) :: KFPSPR,KLGDPR

        REAL(RDT) :: Cptbar,Clgdpr,Cfpspr

        REAL(RDT),ALLOCATABLE :: UMMC(:,:),UCMC(:,:),UKMC(:,:)
        REAL(RDT),ALLOCATABLE :: UFVC(:)

        REAL(RDT),ALLOCATABLE :: UMM(:,:),UCM(:,:),UKM(:,:)
        REAL(RDT),ALLOCATABLE :: UFV(:)

        REAL(RDT) :: BRFCIN(3),BRMTIN(3)
        REAL(RDT) :: BRFCAE(3),BRMTAE(3)
        REAL(RDT) :: BRFC(3),BRMT(3)

        REAL(RDT),ALLOCATABLE :: UQS(:),UQS1T(:),UQS2T(:)

        TYPE(ELEMENT),ALLOCATABLE :: ELEMEN(:)

    CONTAINS

    PROCEDURE,PUBLIC :: CONSTRUCT_BEAM
    PROCEDURE,PUBLIC :: MLEGGAUSS_ST_BEAM
    PROCEDURE,PUBLIC :: MLEGGAUSS_FR_BEAM
    PROCEDURE,PUBLIC :: FORMULATE_ST
    PROCEDURE,PUBLIC :: FORMULATE_FR
    PROCEDURE,PUBLIC :: GIVECONS
    PROCEDURE,PUBLIC :: DECONSTRUCT_BEAM
    PROCEDURE,PUBLIC :: UPDATE
    PROCEDURE,PUBLIC :: PHES75
    PROCEDURE,PUBLIC :: rvndlct
    PROCEDURE,PUBLIC :: rvndlct2
    PROCEDURE,PUBLIC :: GETTHP
    PROCEDURE,PUBLIC :: OUTPUTMATRIX
    PROCEDURE,PUBLIC :: BDSTLD
    PROCEDURE,PRIVATE :: CONTS4
    PROCEDURE,PRIVATE :: CONTS42

    END TYPE BEAM
    !->+++++++++++++++++++++++++++++++++++++++++++++
    CONTAINS
    !->+++++++++++++++++++++++++++++++++++++++++++++
    !->+++++++++++++++++++++++++++++++++++++++++++++
    !->+++++++++++++++PUBLIC++++++++++++++++++++++++++
    !->+++++++++++++++++++++++++++++++++++++++++++++
    !->+++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE CONSTRUCT_BEAM(THIS,INPU)
    USE INPUT_CLASS
    IMPLICIT NONE
    CLASS(BEAM) :: THIS
    CLASS(INPUT) :: INPU

    INTEGER :: I,J,TP,COUNTER

    THIS%BEAMTYPE=INPU%MAT%CHO

    THIS%NELE=INPU%MAT%NBPL-1
    IF(THIS%BEAMTYPE.EQ.0) THEN
        THIS%SMD=14
        THIS%UVC0=6
        THIS%UVC1=2
    ELSE IF(THIS%BEAMTYPE.EQ.1) THEN
        THIS%SMD=23
        THIS%UVC0=9
        THIS%UVC1=5
    END IF
    THIS%UVC2=THIS%UVC0+THIS%UVC1

    THIS%R=INPU%MAT%R
    THIS%EL1=INPU%MAT%EL1
    THIS%BTAP=INPU%IPT%BTAP

    THIS%BTAJ=INPU%MAT%BTA(THIS%NELE)

    THIS%LDTIP=INPU%IPT%LDTIP*PI/180D0

    ALLOCATE(THIS%ELEMEN(THIS%NELE))

    DO I=1,THIS%NELE
        IF(I.EQ.THIS%NELE) THEN
            THIS%ELEMEN(I)%IS_TIP=1
            THIS%ELEMEN(I)%LDT=THIS%LDTIP
        ELSE
            THIS%ELEMEN(I)%IS_TIP=0
            THIS%ELEMEN(I)%LDT=0.0D0
        END IF

        IF(I.GE.INPU%IPT%AFSEN) THEN
            THIS%ELEMEN(I)%IS_AERO=1
        ELSE
            THIS%ELEMEN(I)%IS_AERO=0
        END IF

        CALL THIS%ELEMEN(I)%CONSTRUCT_ELEMENT(I,INPU)
    END DO

    THIS%NDOF=THIS%UVC2*THIS%NELE+THIS%UVC0

    ALLOCATE(THIS%UMM(THIS%NDOF,THIS%NDOF),THIS%UCM(THIS%NDOF,THIS%NDOF))
    ALLOCATE(THIS%UKM(THIS%NDOF,THIS%NDOF),THIS%UFV(THIS%NDOF))

    IF(INPU%IPT%RHS.EQ.0) THEN
        !无铰式旋翼
        IF(THIS%BEAMTYPE.EQ.0) THEN
            THIS%CNT=6
            !去掉根部自由度
            ALLOCATE(THIS%CNTS(THIS%CNT))
            DO I=1,THIS%CNT
                THIS%CNTS(I)=THIS%NDOF-I+1
            END DO
        ELSE IF(THIS%BEAMTYPE.EQ.1) THEN
            !20自由度，每个梁单元去掉2个自由度
            THIS%CNT=8+THIS%NELE*2+1

            ALLOCATE(THIS%CNTS(THIS%CNT))
            THIS%CNTS(1)=7
            THIS%CNTS(2)=12
            COUNTER=1
            !每个梁单元去掉2个自由度
            DO I=1,THIS%NELE
                COUNTER=COUNTER+2
                THIS%CNTS(COUNTER:COUNTER+1)=THIS%CNTS(COUNTER-2:COUNTER-1)+THIS%UVC1+THIS%UVC0
            END DO
            !去掉根部自由度
            DO I=1,9
                THIS%CNTS(THIS%CNT-I+1)=THIS%NDOF-I+1
            END DO
        END IF
    ELSE IF(INPU%IPT%RHS.EQ.1) THEN
        !铰接式旋翼
        THIS%CNT=0
        !最多6个约束
        DO I=1,6
            IF(INPU%IPT%HINGES(I).EQ.0) EXIT
            THIS%CNT=THIS%CNT+1
        END DO
        ALLOCATE(THIS%CNTS(THIS%CNT))
        DO I=1,THIS%CNT
            THIS%CNTS(I)=THIS%NDOF-INPU%IPT%HINGES(I)+1
        END DO
    END IF

    !->+++++++++++++++++++++++++++++++++++++++++++++
    !将CNTS数组的值按从小到大排列
    DO I=1,THIS%CNT
        DO J=I,THIS%CNT
            IF(THIS%CNTS(J).LT.THIS%CNTS(I)) THEN
                TP=THIS%CNTS(J)
                THIS%CNTS(J)=THIS%CNTS(I)
                THIS%CNTS(I)=TP
            END IF
        END DO
    END DO

    THIS%NDOFC=THIS%NDOF-THIS%CNT

    ALLOCATE(THIS%UMMC(THIS%NDOFC,THIS%NDOFC),THIS%UCMC(THIS%NDOFC,THIS%NDOFC))
    ALLOCATE(THIS%UKMC(THIS%NDOFC,THIS%NDOFC),THIS%UFVC(THIS%NDOFC))

    ALLOCATE(THIS%UQS(THIS%NDOFC),THIS%UQS1T(THIS%NDOFC),THIS%UQS2T(THIS%NDOFC))

    THIS%R=INPU%MAT%R

    THIS%npt=THIS%NDOF+(INPU%IPT%nptct-1-INPU%IPT%CNDN)*THIS%uvc2+INPU%IPT%nhpit
    THIS%nfp=THIS%NDOF+(INPU%IPT%nfpct-1-INPU%IPT%CNDN)*THIS%uvc2+INPU%IPT%nhflp
    THIS%nll=THIS%NDOF+(INPU%IPT%nllct-1-INPU%IPT%CNDN)*THIS%uvc2+INPU%IPT%nhlag

    THIS%kptbar=INPU%IPT%kptbar
    THIS%kfp =INPU%IPT%kfp

    THIS%klgdpr =INPU%IPT%klgdpr
    THIS%kfpspr =INPU%IPT%kfpspr

    THIS%Cptbar =INPU%IPT%Cptbar
    THIS%Clgdpr =INPU%IPT%Clgdpr
    THIS%Cfpspr =INPU%IPT%Cfpspr

    THIS%UMMC=0.0D0
    THIS%UCMC=0.0D0
    THIS%UKMC=0.0D0
    THIS%UFVC=0.0D0
    THIS%UQS=0.0D0
    THIS%UQS1T=0.0D0
    THIS%UQS2T=0.0D0

    END SUBROUTINE
    !->++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE FORMULATE_ST(THIS)
    IMPLICIT NONE
    CLASS(BEAM) :: THIS
    INTEGER :: I
    INTEGER :: UVC,SMD
    INTEGER :: MAH2(THIS%SMD)
    REAL(RDT), ALLOCATABLE:: TRLG(:,:),RLGM(:,:),RLGC(:,:),RLGK(:,:),TFC(:),FC(:)

    IF(THIS%BEAMTYPE.EQ.0) THEN
        MAH2= (/9,10,1,2,11,12,3,4,13,7,5,14,8,6/)
    ELSE IF(THIS%BEAMTYPE.EQ.1) THEN
        !            MAH2= (/13,14,1,2,15,16,3,4,17,9,5,18,10,6,19,11,7,20,12,8/)
        MAH2= (/15,16,1,2,17,18,3,4,19,10,5,20,11,6,21,12,7,22,13,8,23,14,9/)
    END IF

    SMD=THIS%SMD

    ALLOCATE(TRLG(SMD,SMD),TFC(SMD))

    THIS%UMM=0.0D0
    THIS%UCM=0.0D0
    THIS%UKM=0.0D0
    THIS%UFV=0.0D0

    TRLG=0.0D0
    TFC=0.0D0

    UVC=0
    DO I=THIS%NELE,1,-1

        TRLG(MAH2(1:SMD),MAH2(1:SMD))=THIS%ELEMEN(I)%RLGM(1:SMD,1:SMD)
        THIS%UMM(UVC+1:UVC+SMD,UVC+1:UVC+SMD)=THIS%UMM(UVC+1:UVC+SMD,UVC+1:UVC+SMD)+TRLG(:,:)

        TRLG(MAH2(1:SMD),MAH2(1:SMD))=THIS%ELEMEN(I)%RLGC(1:SMD,1:SMD)
        THIS%UCM(UVC+1:UVC+SMD,UVC+1:UVC+SMD)=THIS%UCM(UVC+1:UVC+SMD,UVC+1:UVC+SMD)+TRLG(:,:)

        TRLG(MAH2(1:SMD),MAH2(1:SMD))=THIS%ELEMEN(I)%RLGK(1:SMD,1:SMD)
        THIS%UKM(UVC+1:UVC+SMD,UVC+1:UVC+SMD)=THIS%UKM(UVC+1:UVC+SMD,UVC+1:UVC+SMD)+TRLG(:,:)

        TFC(MAH2(1:SMD))=THIS%ELEMEN(I)%FC(1:SMD)
        THIS%UFV(UVC+1:UVC+SMD)=THIS%UFV(UVC+1:UVC+SMD)+TFC(:)

        UVC=UVC+THIS%UVC2
    END DO

    DEALLOCATE(TRLG,TFC)

    !ADD THE EFFECT OF THE CONTROL LINK计入操纵线系的影响

    THIS%UKM(THIS%npt, THIS%npt)=THIS%UKM(THIS%npt, THIS%npt)+THIS%kptbar
    THIS%UKM(THIS%npt, THIS%nfp)=THIS%UKM(THIS%npt, THIS%nfp)+THIS%kfp*THIS%kptbar
    THIS%UKM(THIS%nfp, THIS%npt)=THIS%UKM(THIS%nfp, THIS%npt)+THIS%kfp*THIS%kptbar
    THIS%UKM(THIS%nfp, THIS%nfp)=THIS%UKM(THIS%nfp, THIS%nfp)+THIS%kfp*THIS%kfp*THIS%kptbar

    !ADD THE EFFECT OF THE DAMPER 计入阻尼的影响

    THIS%UCM(THIS%npt, THIS%npt)=THIS%UCM(THIS%npt, THIS%npt)+THIS%Cptbar
    THIS%UCM(THIS%nll, THIS%nll)=THIS%UCM(THIS%nll, THIS%nll)+THIS%Clgdpr
    THIS%UCM(THIS%nfp, THIS%nfp)=THIS%UCM(THIS%nfp, THIS%nfp)+THIS%Cfpspr

    !ADD THE EFFECT OF THE SPRING 计入弹簧的影响

    THIS%UKM(THIS%nll, THIS%nll)=THIS%UKM(THIS%nll, THIS%nll)+THIS%klgdpr
    THIS%UKM(THIS%nfp, THIS%nfp)=THIS%UKM(THIS%nfp, THIS%nfp)+THIS%kfpspr

    CALL THIS%GIVECONS()
    !CALL THIS%OUTPUTMATRIX()
    !pause
    END SUBROUTINE
    !->++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE FORMULATE_FR(THIS)
    IMPLICIT NONE
    CLASS(BEAM) :: THIS
    INTEGER :: I


    THIS%BRFCIN=0.0D0
    THIS%BRMTIN=0.0D0
    THIS%BRFCAE=0.0D0
    THIS%BRMTAE=0.0D0

    DO I=THIS%NELE,1,-1

        THIS%BRFCIN=THIS%BRFCIN+THIS%ELEMEN(I)%FC1
        THIS%BRMTIN=THIS%BRMTIN+THIS%ELEMEN(I)%MT1

        THIS%BRFCAE=THIS%BRFCAE+THIS%ELEMEN(I)%FC2
        THIS%BRMTAE=THIS%BRMTAE+THIS%ELEMEN(I)%MT2

    END DO

    THIS%BRFC=THIS%BRFCIN+THIS%BRFCAE
    THIS%BRMT=THIS%BRMTIN+THIS%BRMTAE


    END SUBROUTINE
    !->++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE MLEGGAUSS_ST_BEAM(THIS,TH0,TH1C,TH1S,PSI,OMG,SOL)
    IMPLICIT NONE
    CLASS(BEAM) :: THIS
    INTEGER,INTENT(IN)::SOL
    REAL(RDT),INTENT(IN)::TH0,TH1C,TH1S,PSI,OMG

    INTEGER :: I

    REAL(RDT),EXTERNAL:: FABSC,triseries

    CALL THIS%GETTHP(TH0, TH1C, TH1S, PSI, OMG)

    THIS%ELEMEN(THIS%NELE)%NOD(2)%CTFGRL=0.0D0
    DO I=THIS%NELE,1,-1
        CALL THIS%ELEMEN(I)%MLEGGAUSS_ST(THIS%R,THIS%EL1,THIS%BTAP,THIS%THP,THIS%THP1D,THIS%THP2D,OMG,SOL)
        IF(I.GE.2) THEN
            THIS%ELEMEN(I-1)%NOD(2)%CTFGRL=THIS%ELEMEN(I)%NOD(1)%CTFGRL
        END IF
    END DO

    END SUBROUTINE
    !->++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE MLEGGAUSS_FR_BEAM(THIS,TH0,TH1C,TH1S,PSI,OMG,SOL)
    IMPLICIT NONE
    CLASS(BEAM) :: THIS
    INTEGER,INTENT(IN)::SOL
    REAL(RDT),INTENT(IN)::TH0,TH1C,TH1S,PSI,OMG

    INTEGER :: I

    REAL(RDT),EXTERNAL:: FABSC,triseries

    CALL THIS%GETTHP(TH0, TH1C, TH1S, PSI, OMG)

    DO I=THIS%NELE,1,-1
        CALL THIS%ELEMEN(I)%MLEGGAUSS_FR(THIS%R,THIS%EL1,THIS%BTAP,THIS%THP,THIS%THP1D,THIS%THP2D,OMG,SOL)
    END DO

    END SUBROUTINE
    !->++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE BDSTLD(THIS,AIM,TH0,TH1C,TH1S,PSI,OMG,STM)
    USE MBEAMAERO
    IMPLICIT NONE
    CLASS(BEAM) :: THIS
    REAL(RDT),INTENT(IN)::AIM
    REAL(RDT),INTENT(IN)::TH0,TH1C,TH1S,PSI,OMG
    REAL(RDT),INTENT(OUT)::STM(3)

    INTEGER :: I

    REAL(RDT) SS,LL,RXUD,he

    INTEGER NBPLI,NBPLIP1

    REAL(RDT),EXTERNAL:: FABSC,triseries

    CALL THIS%GETTHP(TH0, TH1C, TH1S, PSI, OMG)

    rxud=AIM
    call rxlct(nbpli, ss, ll, NNOD, rxud, rlud)

    ll=ll*THIS%R
    HE=THIS%ELEMEN(NBPLI)%NOD(1)%RLUD*THIS%R-THIS%EL1

    NBPLIP1=NBPLI+1
    CALL THIS%ELEMEN(NBPLI)%GETLV(HE,SS,LL,THIS%BTAP,THIS%EL1,THIS%THP,THIS%THP1D,THIS%THP2D,OMG)
    CALL THIS%ELEMEN(NBPLI)%UMINT(STM)

    END SUBROUTINE
    !->+++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE GIVECONS(THIS)
    IMPLICIT NONE
    CLASS(BEAM) :: THIS
    INTEGER :: NDOF,NDOFC,NELE

    INTEGER :: CCDFSN(THIS%CNT)

    NDOF=THIS%NDOF
    NDOFC=THIS%NDOFC
    NELE=THIS%NELE

    CCDFSN=THIS%CNTS(1:THIS%CNT)
    CALL THIS%CONTS4(THIS%UKMC,NDOF,NDOFC,0,CCDFSN,THIS%UKM)
    CALL THIS%CONTS4(THIS%UMMC,NDOF,NDOFC,0,CCDFSN,THIS%UMM)
    CALL THIS%CONTS4(THIS%UCMC,NDOF,NDOFC,0,CCDFSN,THIS%UCM)
    CALL THIS%CONTS42(THIS%UFVC,NDOF,NDOFC,0,CCDFSN,THIS%UFV)

    END SUBROUTINE
    !->+++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE DECONSTRUCT_BEAM(THIS)
    IMPLICIT NONE
    CLASS(BEAM) :: THIS

    DEALLOCATE(THIS%UMMC,THIS%UCMC,THIS%UKMC)
    DEALLOCATE(THIS%UFVC)
    DEALLOCATE(THIS%UMM,THIS%UCM,THIS%UKM)
    DEALLOCATE(THIS%UFV)
    DEALLOCATE(THIS%UQS,THIS%UQS1T,THIS%UQS2T)
    DEALLOCATE(THIS%ELEMEN)

    END SUBROUTINE
    !->++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE UPDATE(THIS)
    IMPLICIT NONE
    CLASS(BEAM) :: THIS

    INTEGER :: I,COUNTER,ST,ED
    REAL(RDT),ALLOCATABLE :: USP(:),US1TP(:),US2TP(:)
    INTEGER,ALLOCATABLE :: MARK(:)

    ALLOCATE(USP(THIS%NDOF),US1TP(THIS%NDOF),US2TP(THIS%NDOF),MARK(THIS%NDOF))
    MARK=0
    MARK(THIS%CNTS(1:THIS%CNT))=1

    USP=0.0D0
    US1TP=0.0D0
    US2TP=0.0D0

    COUNTER=0
    DO I=1,THIS%NDOF
        IF(MARK(I).EQ.0) THEN
            COUNTER=COUNTER+1
            USP(I)=THIS%UQS(COUNTER)
            US1TP(I)=THIS%UQS1T(COUNTER)
            US2TP(I)=THIS%UQS2T(COUNTER)
        ELSE
            USP(I)=0.0D0
            US1TP(I)=0.0D0
            US2TP(I)=0.0D0
        END IF
    END DO

    DO I=1,THIS%NELE
        ST=(I-1)*THIS%UVC2+1
        ED=I*THIS%UVC2+THIS%UVC0
        CALL THIS%ELEMEN(THIS%NELE-I+1)%UPDATE_ELE(THIS%UVC2+THIS%UVC0,USP(ST:ED),US1TP(ST:ED),US2TP(ST:ED))
    END DO

    END SUBROUTINE
    !->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine rvndlct(THIS,rnd, vnd,npt,RXUD,cdlct,OMG,PSI,CHO,CHO2)
    !->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !-> modules used:
    USE MBEAMAERO
    implicit none
    CLASS(BEAM) :: THIS
    !->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !-> dummy arguments / interface:
    INTEGER,INTENT(IN) :: NPT,CHO,CHO2
    REAL(RDT),INTENT(IN) :: OMG
    REAL(RDT),INTENT(IN)::PSI,RXUD, cdlct(3, npt)
    real(rdt), intent(out):: rnd(3, npt), vnd(3, npt)

    !->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !-> FUNCTIONS INVOKED:
    !->
    REAL(RDT), EXTERNAL:: VIP, FWM, FABSC, triseries
    !->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !-> SUBROUTINES INVOKED:
    !->
    !->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !-> LOCAL VARIABLES:
    integer nbpli,nbplip1,I
    real(rdt) ss, ll
    real(rdt) he
    REAL(RDT) R0(3),r0b(3), r0e(3)

    real(rdt) rc(3, npt), vc(3, npt),rd(3, npt),vd(3, npt),TP1(3,npt),TP2(3,npt)

    integer mk,dm
    real(rdt) trn(3, 3),ter(3, 3), tnr(3, 3), tne(3, 3), tre(3,3)

    !->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !-> main body:
    mk=1
    dm=3

    call rxlct(nbpli, ss, ll, NNOD, rxud, rlud*THIS%R)

    ll=ll
    HE=THIS%ELEMEN(NBPLI)%NOD(1)%RLUD*THIS%R-THIS%EL1

    NBPLIP1=NBPLI+1
    CALL THIS%ELEMEN(NBPLI)%GETLV(HE,SS,LL,THIS%BTAP,THIS%EL1,THIS%THP,THIS%THP1D,THIS%THP2D,OMG)


    R0=THIS%EL1*(/1.0D0, 0.0D0, 0.0D0/)
    CALL MVMU2(R0B,3,3,THIS%ELEMEN(NBPLI)%TBR,R0)
    CALL MVMU2(R0E,3,3,THIS%ELEMEN(NBPLI)%TEB,R0B)

    IF(THIS%ELEMEN(NBPLI)%IS_TIP.EQ.1) THEN
        CALL MAMU(TP1, 3, 3, npt, THIS%ELEMEN(NBPLI)%TEB, cdlct)
    ELSE
        TP1=cdlct
    END IF
    IF (CHO2 == 1) THEN
        call THIS%ELEMEN(NBPLI)%rvclct(rc, vc, dm, npt, TP1,1)
    ELSE
        call THIS%ELEMEN(NBPLI)%rvclct(rc, vc, dm, npt, TP1)
    END IF

    call avadd(rc, r0e, dm, npt, dm, mk)
    call avadd(vc, THIS%ELEMEN(NBPLI)%vbe, dm, npt, dm, mk)

    CALL MAMU(TER,dm,dm,dm,THIS%ELEMEN(NBPLI)%TEB,THIS%ELEMEN(NBPLI)%TBR)

    CALL MAMU(rd, 3, 3, npt, TRANSPOSE(TER), rc)
    CALL MAMU(vd, 3, 3, npt, TRANSPOSE(TER), vc)

    RD(1,2:npt)=RD(1,1)

    IF(CHO.EQ.0) THEN
        CALL TXYZ(trn, dm, psi)
        TNR=TRANSPOSE(TRN)

        CALL MAMU(rnd, 3, 3, npt, TNR, rd)
        CALL MAMU(vnd, 3, 3, npt, TNR, vd)
    ELSE IF(CHO.EQ.1) THEN
        RETURN
    END IF

    end subroutine rvndlct
    !->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine rvndlct2(THIS,REA,TED,TRE,RXUD,OMG,CHO2)
    !->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !-> modules used:
    USE MBEAMAERO
    implicit none
    CLASS(BEAM) :: THIS
    !->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !-> dummy arguments / interface:
    INTEGER,INTENT(IN) :: CHO2
    REAL(RDT),INTENT(IN) :: OMG
    REAL(RDT),INTENT(IN):: RXUD
    real(rdt), intent(out):: REA(3), TED(3, 3), TRE(3, 3)

    !->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !-> FUNCTIONS INVOKED:
    !->
    REAL(RDT), EXTERNAL:: VIP, FWM, FABSC, triseries
    !->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !-> SUBROUTINES INVOKED:
    !->
    !->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !-> LOCAL VARIABLES:
    integer nbpli,nbplip1,I
    real(rdt) ss, ll
    real(rdt) he
    REAL(RDT) R0(3),r0b(3), r0e(3)


    integer mk,dm
    real(rdt) trn(3, 3),ter(3, 3), tnr(3, 3), tne(3, 3)

    !->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !-> main body:
    mk=1
    dm=3

    !call rxlct(nbpli, ss, ll, NNOD, rxud, rlud*THIS%R)
    call rxlct(nbpli, ss, ll, NNOD, rxud, rlud)!LSJ
    
    !ll=ll
    ll=ll*THIS%R!LSJ
    HE=THIS%ELEMEN(NBPLI)%NOD(1)%RLUD*THIS%R-THIS%EL1

    NBPLIP1=NBPLI+1
    !CALL THIS%ELEMEN(NBPLI)%GETLV(HE,SS,LL,THIS%BTAP,THIS%EL1,THIS%THP,THIS%THP1D,THIS%THP2D,OMG)
    CALL THIS%ELEMEN(NBPLI)%GETLV(HE,SS,LL,THIS%BTAP,THIS%EL1,0.0D0,0.0D0,0.0D0,OMG)!不输出操纵量


    R0=THIS%EL1*(/1.0D0, 0.0D0, 0.0D0/)
    CALL MVMU2(R0B,3,3,THIS%ELEMEN(NBPLI)%TBR,R0)
    CALL MVMU2(R0E,3,3,THIS%ELEMEN(NBPLI)%TEB,R0B)

    IF (CHO2 == 1) THEN
        call THIS%ELEMEN(NBPLI)%rvclct2(REA, TED, dm, 1)
    ELSE
        call THIS%ELEMEN(NBPLI)%rvclct2(REA, TED, dm)
    END IF

    CALL MAMU(TER,dm,dm,dm,THIS%ELEMEN(NBPLI)%TEB,THIS%ELEMEN(NBPLI)%TBR)
    TRE=TRANSPOSE(TER)

    end subroutine rvndlct2
    !->+++++++++++++++++++++++++++++++++++++++++++++
    !->+++++++++++++++++++++++++++++++++++++++++++++
    !->+++++++++++++++PRIVATE+++++++++++++++++++++++++
    !->+++++++++++++++++++++++++++++++++++++++++++++
    !->+++++++++++++++++++++++++++++++++++++++++++++
    !->++++++++++++++++++++++++++++++++++++++++++++
    REAL(RDT) FUNCTION PHES75(THIS)
    USE MBEAMAERO
    IMPLICIT NONE
    CLASS(BEAM) :: THIS
    INTEGER NBPLI,NBPLIP1
    REAL(RDT) :: SS,LL

    CALL THIS%UPDATE()
    DO NBPLI=NNOD-1,1,-1
        NBPLIP1=NBPLI+1
        IF(RLUD(NBPLI) .LE. 0.75) THEN
            LL=RLUD(NBPLIP1)-RLUD(NBPLI)
            SS=(0.75-RLUD(NBPLI))/LL
            LL=LL*THIS%R
            EXIT
        END IF
    END DO

    CALL THIS%ELEMEN(NBPLI)%GETLVSS(SS,LL)

    PHES75=THIS%ELEMEN(NBPLI)%BTAT+THIS%ELEMEN(NBPLI)%SPH

    END FUNCTION
    !->++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE GETTHP(THIS,TH0,TH1C,TH1S,PSI,OMG)
    USE MBEAMAERO
    IMPLICIT NONE
    CLASS(BEAM) :: THIS
    REAL(RDT),INTENT(IN) :: TH0,TH1C,TH1S,PSI,OMG

    REAL(RDT),EXTERNAL:: FABSC,triseries

    THIS%THP=triseries(TH0, TH1C, TH1S, PSI)
    THIS%THP1D=OMG*FABSC(TH1S,-TH1C,PSI)
    THIS%THP2D=-OMG**2*(THIS%THP-TH0)

    END SUBROUTINE
    !->++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE OUTPUTMATRIX(THIS)
    IMPLICIT NONE
    CLASS(BEAM) :: THIS
    INTEGER :: I,J

    open(11,file='tk.dat')
    open(12,file='tm.dat')
    open(13,file='tc.dat')
    open(14,file='tf.dat')
    do i=1,THIS%NDOFC
        write(11,995) (THIS%ukmc(i,j),j=1,THIS%NDOFC)
        write(12,995) (THIS%ummc(i,j),j=1,THIS%NDOFC)
        write(13,995) (THIS%ucmc(i,j),j=1,THIS%NDOFC)
        write(14,995) THIS%ufvc(i)
    end do

    close(11)
    close(12)
    close(13)
    close(14)

995 format(400f26.12)

    END SUBROUTINE

    SUBROUTINE CONTS4(THIS,UMMC,UQD,UQDC,UNC,CCDFSN,UMM)
    implicit none
    CLASS(BEAM) :: THIS
    INTEGER,INTENT(IN)::UQD,UQDC,UNC,CCDFSN(UQD-UQDC)
    REAL(RDT),INTENT(OUT)::UMMC(UQDC,UQDC)
    REAL(RDT),INTENT(IN)::UMM(UQD,UQD)

    INTEGER CCTN
    INTEGER	 I, J, PT
    INTEGER ERR

    REAL(RDT), ALLOCATABLE:: UMMT(:, :)

    allocate(UMMT(UQD,UQD),&
        stat=err)
    if(err.ne.0) then
        write(*, 999)
        stop
    end if

    UMMT=UMM

    CCTN=UQD-UQDC!根部约束个数

    PT=UNC+CCDFSN(1)

    DO 10 I=1, CCTN-1
        DO 20 J=UNC+CCDFSN(I)+1,UNC+CCDFSN(I+1)-1
            UMMT(PT,:)=UMMT(J,:)
            UMMT(:,PT)=UMMT(:,J)

            PT=PT+1
20      CONTINUE
10  CONTINUE

    DO 30 J=UNC+CCDFSN(CCTN)+1,UQD
        UMMT(PT,:)=UMMT(J,:)
        UMMT(:,PT)=UMMT(:,J)

        PT=PT+1
30  CONTINUE


    !->write(*,*) PT
    !>pause 888

    UMMC=UMMT(1:UQDC,1:UQDC)

    deallocate(UMMT,  &
        stat=err)
    if(err.ne.0) then
        write(*, 999)
        stop
    end if

998 FORMAT(3X,'DEALLOCATING ARRAY ERROR, PROGRAM STOPPED!')
999 FORMAT(3X,'ALLOCATING ARRAY ERROR, PROGRAM STOPPED!')
    END



    SUBROUTINE CONTS42(THIS,UFVC,UQD,UQDC,UNC,CCDFSN,UFV)
    implicit none
    CLASS(BEAM) :: THIS
    INTEGER,INTENT(IN)::UQD,UQDC,UNC,CCDFSN(UQD-UQDC)
    REAL(RDT),INTENT(OUT)::UFVC(UQDC)
    REAL(RDT),INTENT(IN)::UFV(UQD)

    INTEGER CCTN
    INTEGER	 I, J, PT
    INTEGER ERR

    REAL(RDT), ALLOCATABLE:: UFVT(:)

    allocate(UFVT(UQD), &
        stat=err)
    if(err.ne.0) then
        write(*, 999)
        stop
    end if


    UFVT=UFV

    CCTN=UQD-UQDC

    PT=UNC+CCDFSN(1)

    DO 10 I=1, CCTN-1
        DO 20 J=UNC+CCDFSN(I)+1,UNC+CCDFSN(I+1)-1
            UFVT(PT)=UFV(J)

            PT=PT+1
20      CONTINUE
10  CONTINUE

    DO 30 J=UNC+CCDFSN(CCTN)+1,UQD
        UFVT(PT)=UFV(J)

        PT=PT+1
30  CONTINUE


    !->write(*,*) PT
    !>pause 888

    UFVC=UFVT(1:UQDC)

    deallocate(UFVT, &
        stat=err)
    if(err.ne.0) then
        write(*, 999)
        stop
    end if

998 FORMAT(3X,'DEALLOCATING ARRAY ERROR, PROGRAM STOPPED!')
999 FORMAT(3X,'ALLOCATING ARRAY ERROR, PROGRAM STOPPED!')
    END
    END MODULE
