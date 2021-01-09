    MODULE EIGENBEAM_CLASS
    USE EIGEN_CLASS
    IMPLICIT NONE
    !->+++++++++++++++++++++++++++++++++++++++++++++
    TYPE,PUBLIC,EXTENDS(EIGEN) :: EIGENBEAM

        INTEGER :: NATMOD

        INTEGER :: NHF,NHF2

        REAL(RDT), ALLOCATABLE:: q2t(:), q1t(:), q(:)
        REAL(RDT), ALLOCATABLE:: q2tp(:), q1tp(:), qp(:)

        REAL(RDT), ALLOCATABLE:: mmtinv(:,:), mmt(:, :), cmt(:, :), kmt(:, :), fv(:), fvx(:)

        REAL(RDT), ALLOCATABLE:: yst(:), y(:),yp(:)

        REAL(RDT) :: rdcm, rdck

    CONTAINS

    PROCEDURE,PUBLIC :: CONSTRUCT_EIGENBEAM
    PROCEDURE,PUBLIC :: INITIAL_UQSX
    PROCEDURE,PUBLIC :: GETMKCFT
    PROCEDURE,PUBLIC :: GETUQS
    PROCEDURE,PUBLIC :: GETSTATIC

    END TYPE EIGENBEAM
    !->+++++++++++++++++++++++++++++++++++++++++++++
    !->+++++++++++++++++++++++++++++++++++++++++++++
    !->+++++++++++++++PUBLIC++++++++++++++++++++++++++
    !->+++++++++++++++++++++++++++++++++++++++++++++
    !->+++++++++++++++++++++++++++++++++++++++++++++
    CONTAINS
    !->+++++++++++++++++++++++++++++++++++++++++++++
    subroutine CONSTRUCT_EIGENBEAM(THIS,TH0,TH1C,TH1S,OMG,INPU)
    USE INPUT_CLASS
    implicit none
    CLASS(EIGENBEAM) :: THIS
    TYPE(INPUT) :: INPU
    REAL(RDT),INTENT(IN) :: TH0,TH1C,TH1S,OMG

    CALL THIS%CONSTRUCT_EIGEN(INPU)

    THIS%natmod=INPU%FWD%natmod

    THIS%rdcm=INPU%IPT%rdcm
    THIS%rdck=INPU%IPT%rdck

    if( this%natmod .eq. 0 ) then
        THIS%nhf=THIS%NDOFC
        THIS%nhf2=2*THIS%NDOFC
    else if( this%natmod .eq. 1 ) then
        THIS%nhf=THIS%NMN
        THIS%nhf2=2*THIS%NMN
    end if

    ALLOCATE(THIS%q2t(THIS%nhf), THIS%q1t(THIS%nhf), THIS%q(THIS%nhf))
    ALLOCATE(THIS%q2tp(THIS%nhf), THIS%q1tp(THIS%nhf), THIS%qp(THIS%nhf))

    ALLOCATE(THIS%MMTINV(THIS%nhf,THIS%nhf), THIS%mmt(THIS%nhf, THIS%nhf),&
        THIS%cmt(THIS%nhf, THIS%nhf), THIS%kmt(THIS%nhf, THIS%nhf), &
        THIS%fv(THIS%nhf), THIS%fvx(THIS%nhf))

    ALLOCATE(THIS%yst(THIS%nhf2), THIS%y(THIS%nhf2), THIS%yp(THIS%nhf2))

    CALL THIS%GETVTEG(TH0,TH1C,TH1S,OMG)

    END SUBROUTINE
    !->+++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE INITIAL_UQSX(THIS)
    IMPLICIT NONE
    CLASS(EIGENBEAM) :: THIS

    IF(THIS%NATMOD.EQ.0) THEN

        call THIS%GETSTATIC()

        call mtxinv2(THIS%MMTINV, 1, THIS%NHF, THIS%UMMC)

        THIS%yst=0.0D0
        THIS%yst(1:THIS%nhf)=THIS%UQS
        THIS%Y=0.0D0
        THIS%YP=0.0D0

    ELSE IF(THIS%NATMOD.EQ.1) THEN

        CALL THIS%GETMKCFT()
        THIS%yst=0.0D0
        call linsolver2(THIS%yst, 1, THIS%NMN, THIS%kmt, THIS%fv)
        THIS%Y=0.0D0
        THIS%YP=0.0D0
        call mtxinv2(THIS%MMTINV, 1, THIS%NMN, THIS%mmt)

        CALL MVMU2(THIS%UQS, THIS%NDOFC, THIS%NMN, THIS%VTEG, THIS%YST)
        THIS%UQS1T=0.0D0
        THIS%UQS2T=0.0D0

    END IF

    END SUBROUTINE
    !->+++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE GETMKCFT(THIS)
    IMPLICIT NONE
    CLASS(EIGENBEAM) :: THIS

    CALL MATMLYTP(THIS%mmt, THIS%NMN, THIS%NDOFC, THIS%NDOFC, THIS%NMN, THIS%VTEGT, THIS%UMMC, THIS%VTEG)
    CALL MATMLYTP(THIS%kmt, THIS%NMN, THIS%NDOFC, THIS%NDOFC, THIS%NMN, THIS%VTEGT, THIS%UKMC, THIS%VTEG)
    CALL MATMLYTP(THIS%cmt, THIS%NMN, THIS%NDOFC, THIS%NDOFC, THIS%NMN, THIS%VTEGT, THIS%UCMC, THIS%VTEG)
    CALL MVMU2(THIS%fv, THIS%NMN, THIS%NDOFC, THIS%VTEGT, -THIS%UFVC)

    END SUBROUTINE
    !->+++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE GETUQS(THIS)
    IMPLICIT NONE
    CLASS(EIGENBEAM) :: THIS

    IF(THIS%NATMOD.EQ.1) THEN
        CALL MVMU2(THIS%UQS, THIS%NDOFC, THIS%NMN, THIS%VTEG, THIS%Y)
        CALL MVMU2(THIS%UQS1T, THIS%NDOFC, THIS%NMN, THIS%VTEG, THIS%Y(THIS%NMN+1))
        CALL MVMU2(THIS%UQS2T, THIS%NDOFC, THIS%NMN, THIS%VTEG, THIS%YP(THIS%NMN+1))
    ELSE IF(THIS%NATMOD.EQ.0) THEN
        THIS%UQS=THIS%Y(1:THIS%NHF)
        THIS%UQS1T=THIS%YP(1:THIS%NHF)
        THIS%UQS2T=THIS%YP(THIS%NHF+1:THIS%NHF2)
    END IF

    END SUBROUTINE
    !->+++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE GETSTATIC(THIS)
    IMPLICIT NONE
    CLASS(EIGENBEAM) :: THIS

    IF(THIS%NATMOD.EQ.0) THEN
        call linsolver2(THIS%UQS, 1, THIS%NDOFC, THIS%UKMC, -THIS%UFVC )
    ELSE
        CALL THIS%GETMKCFT()
        call linsolver2(THIS%Q, 1, THIS%NMN, THIS%kmt, THIS%fv )
        CALL MVMU2(THIS%UQS, THIS%NDOFC, THIS%NMN, THIS%VTEG, THIS%Q)
    END IF

    THIS%Q1T=0.0D0
    THIS%Q2T=0.0D0

    THIS%UQS1T=0.0D0
    THIS%UQS2T=0.0D0

    END SUBROUTINE
    END MODULE
