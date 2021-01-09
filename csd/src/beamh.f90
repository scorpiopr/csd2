    MODULE BEAMH_CLASS
    USE BEAM_CLASS
    IMPLICIT NONE
    !->+++++++++++++++++++++++++++++++++++++++++++++
    TYPE,PUBLIC,EXTENDS(BEAM) :: BEAMH

        INTEGER :: NDOFCC

        REAL(RDT),ALLOCATABLE :: UMMCC(:,:),UCMCC(:,:),UKMCC(:,:)
        REAL(RDT),ALLOCATABLE :: UFVCC(:)

    CONTAINS

    PROCEDURE,PUBLIC :: CONSTRUCT_BEAMH
    PROCEDURE,PUBLIC :: GIVECONSH
    !    PROCEDURE,PUBLIC :: DECONSTRUCT_BEAM
    PROCEDURE,PRIVATE :: CONTS4_H
    PROCEDURE,PRIVATE :: CONTS42_H

    END TYPE BEAMH
    !->+++++++++++++++++++++++++++++++++++++++++++++
    CONTAINS
    SUBROUTINE CONSTRUCT_BEAMH(THIS,INPU)
    USE INPUT_CLASS
    IMPLICIT NONE
    CLASS(BEAMH) :: THIS
    CLASS(INPUT) :: INPU

    CALL THIS%CONSTRUCT_BEAM(INPU)

    THIS%NDOFCC=INPU%MAT%NBPL*2-1
    ALLOCATE(THIS%UMMCC(THIS%NDOFCC,THIS%NDOFCC),&
        THIS%UCMCC(THIS%NDOFCC,THIS%NDOFCC),&
        THIS%UKMCC(THIS%NDOFCC,THIS%NDOFCC),&
        THIS%UFVCC(THIS%NDOFCC))

    END SUBROUTINE
    !->+++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE GIVECONSH(THIS)
    IMPLICIT NONE
    CLASS(BEAMH) :: THIS
    INTEGER :: NDOFC,NDOFCC,NNOD

    NDOFCC=THIS%NDOFCC
    NDOFC=THIS%NDOFC
    NNOD=THIS%NELE+1

    CALL THIS%CONTS4_H(THIS%UMMCC,NDOFCC,NDOFC,NNOD,THIS%UMMC)
    CALL THIS%CONTS4_H(THIS%UKMCC,NDOFCC,NDOFC,NNOD,THIS%UKMC)
    CALL THIS%CONTS4_H(THIS%UCMCC,NDOFCC,NDOFC,NNOD,THIS%UCMC)
    CALL THIS%CONTS42_H(THIS%UFVCC,NDOFCC,NDOFC,NNOD,THIS%UFVC)

    END SUBROUTINE
    !->+++++++++++++++++++++++++++++++++++++++++++++

    SUBROUTINE CONTS4_H(THIS,UMMCC,UQDCC,UQDC,NBPL,UMMC)
    implicit none
    CLASS(BEAMH) :: THIS
    INTEGER,INTENT(IN)::UQDC,NBPL
    INTEGER,INTENT(IN)::UQDCC
    REAL(RDT),INTENT(IN)::UMMC(UQDC,UQDC)
    REAL(RDT),INTENT(OUT)::UMMCC(UQDCC,UQDCC)



    INTEGER CCTN
    INTEGER	 I, J, PT
    INTEGER ERR

    REAL(RDT), ALLOCATABLE:: UMMT(:, :)

    allocate(UMMT(UQDC,UQDC), &
        stat=err)
    if(err.ne.0) then
        write(*, 999)
        stop
    end if

    UMMT=UMMC

    PT=1
    J=3
    UMMT(PT:PT+1,:)=UMMT(J:J+1,:)
    UMMT(:,PT:PT+1)=UMMT(:,J:J+1)

    DO 10 I=1,NBPL-1
        IF(I.EQ.NBPL-1) THEN
            J=I*8+2
            PT=PT+2
            UMMT(PT,:)=UMMT(J,:)
            UMMT(:,PT)=UMMT(:,J)
        ELSE
            J=I*8+3
            PT=PT+2
            UMMT(PT:PT+1,:)=UMMT(J:J+1,:)
            UMMT(:,PT:PT+1)=UMMT(:,J:J+1)
        END IF
10  CONTINUE


    IF(UQDCC.NE.PT) THEN
        WRITE(*,*) UQDCC,PT
        PAUSE
        STOP
    END IF
    UMMCC(1:UQDCC,1:UQDCC)=UMMT(1:UQDCC,1:UQDCC)

    deallocate(UMMT, &
        stat=err)
    if(err.ne.0) then
        write(*, 999)
        stop
    end if

998 FORMAT(3X,'DEALLOCATING ARRAY ERROR, PROGRAM STOPPED!')
999 FORMAT(3X,'ALLOCATING ARRAY ERROR, PROGRAM STOPPED!')
    END

    SUBROUTINE CONTS42_H(THIS,UFVCC,UQDCC,UQDC,NBPL,UFVC)
    implicit none
    CLASS(BEAMH) :: THIS
    INTEGER,INTENT(IN)::UQDC,NBPL
    INTEGER,INTENT(IN)::UQDCC
    REAL(RDT),INTENT(IN)::UFVC(UQDC)
    REAL(RDT),INTENT(OUT)::UFVCC(UQDCC)

    INTEGER CCTN
    INTEGER	 I, J, PT
    INTEGER ERR

    REAL(RDT), ALLOCATABLE:: UFVT(:)

    allocate(UFVT(UQDC), &
        stat=err)
    if(err.ne.0) then
        write(*, 999)
        stop
    end if


    UFVT=UFVC

    PT=1
    J=3

    UFVT(PT:PT+1)=UFVT(J:J+1)

    DO 10 I=1,NBPL-1
        IF(I.EQ.NBPL-1) THEN
            J=I*8+2
            PT=PT+2
            UFVT(PT)=UFVT(J)
        ELSE
            J=I*8+3
            PT=PT+2
            UFVT(PT:PT+1)=UFVT(J:J+1)
        END IF
10  CONTINUE


    IF(UQDCC.NE.PT) THEN
        WRITE(*,*) UQDCC,PT
        PAUSE
        STOP
    END IF


    UFVCC(1:UQDCC)=UFVT(1:UQDCC)

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

