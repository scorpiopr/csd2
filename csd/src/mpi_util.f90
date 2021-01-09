!MODULE mpi_class
!implicit none
!include "mpif.h"
!TYPE,PUBLIC :: MPI
!
!    INTEGER :: IERR,MYCOM
!    INTEGER :: myid, NPROC
!    INTEGER :: status,tag
!    INTEGER,ALLOCATABLE :: IREQR(:),IREQS(:)
!
!    CONTAINS
!
!    PROCEDURE,PUBLIC :: MPI_INITIAL
!    PROCEDURE,PUBLIC :: BCASTC
!    PROCEDURE,PUBLIC :: BCASTS
!    PROCEDURE,PUBLIC :: BCAST
!    PROCEDURE,PUBLIC :: REDUCE
!    PROCEDURE,PUBLIC :: REDUCEI
!    PROCEDURE,PUBLIC :: REDUCED
!    PROCEDURE,PUBLIC :: SEND_RECV
!    PROCEDURE,PUBLIC :: SEND
!    PROCEDURE,PUBLIC :: SENDI
!    PROCEDURE,PUBLIC :: ISEND
!    PROCEDURE,PUBLIC :: RECV
!    PROCEDURE,PUBLIC :: RECVI
!    PROCEDURE,PUBLIC :: IRECV
!    PROCEDURE,PUBLIC :: FINALIZE
!    PROCEDURE,PUBLIC :: BARRIER
!    PROCEDURE,PUBLIC :: GATHER
!    PROCEDURE,PUBLIC :: WAITS
!    PROCEDURE,PUBLIC :: WAITR
!    PROCEDURE,PUBLIC :: ISEND_RECV
!END TYPE MPI
!
!CONTAINS
!SUBROUTINE MPI_INITIAL(THIS)
!IMPLICIT NONE
!CLASS(MPI) :: THIS
!
!    CALL MPI_INIT(THIS%IERR)
!    CALL MPI_COMM_DUP(MPI_COMM_WORLD,THIS%MYCOM,THIS%IERR)
!    call MPI_COMM_RANK(THIS%MYCOM, THIS%myid, THIS%ierr)  
!    call MPI_COMM_SIZE(THIS%MYCOM, THIS%NPROC, THIS%ierr)
!
!    ALLOCATE(THIS%IREQR(THIS%NPROC),THIS%IREQS(THIS%NPROC))
!
!    WRITE(*,*) 'PROCESS  ',THIS%MYID, ' OF ', THIS%NPROC, ' IS RUNNING  '
!
!END SUBROUTINE
!!->+++++++++++++++++++++++++++++++++++++++++++++
!SUBROUTINE BCAST(THIS,N,I)
!IMPLICIT NONE
!CLASS(MPI) :: THIS
!INTEGER,INTENT(IN) :: N,I
!
!    call MPI_BCAST(N,1,MPI_INTEGER,I,THIS%MYCOM,THIS%ierr)
!
!END SUBROUTINE
!!->+++++++++++++++++++++++++++++++++++++++++++++
!SUBROUTINE BCASTS(THIS,N,M,I)
!IMPLICIT NONE
!CLASS(MPI) :: THIS
!INTEGER,INTENT(IN) :: M,I
!INTEGER,INTENT(IN) :: N(M)
!
!    call MPI_BCAST(N,M,MPI_INTEGER,I,THIS%MYCOM,THIS%ierr)
!
!END SUBROUTINE
!!->+++++++++++++++++++++++++++++++++++++++++++++
!SUBROUTINE BCASTC(THIS,N)
!IMPLICIT NONE
!CLASS(MPI) :: THIS
!CHARACTER,INTENT(IN) :: N
!
!    call MPI_BCAST(N,1,MPI_CHARACTER,0,THIS%MYCOM,THIS%ierr)
!
!END SUBROUTINE
!!->+++++++++++++++++++++++++++++++++++++++++++++
!SUBROUTINE REDUCE(THIS,in,out)
!IMPLICIT NONE
!CLASS(MPI) :: THIS
!REAL(8),INTENT(IN) :: IN
!REAL(8),INTENT(OUT) :: OUT
!
!    call MPI_ALLREDUCE(IN,OUT,1,MPI_DOUBLE_PRECISION,MPI_SUM,THIS%MYCOM,THIS%ierr) 
!
!END SUBROUTINE
!!->+++++++++++++++++++++++++++++++++++++++++++++
!SUBROUTINE REDUCEI(THIS,in,out)
!IMPLICIT NONE
!CLASS(MPI) :: THIS
!INTEGER,INTENT(IN) :: IN
!INTEGER,INTENT(OUT) :: OUT
!
!    call MPI_ALLREDUCE(IN,OUT,1,MPI_INTEGER,MPI_SUM,THIS%MYCOM,THIS%ierr) 
!
!END SUBROUTINE
!!->+++++++++++++++++++++++++++++++++++++++++++++
!SUBROUTINE REDUCED(THIS,N,in,out)
!IMPLICIT NONE
!CLASS(MPI) :: THIS
!INTEGER,INTENT(IN) :: N
!REAL(8),INTENT(IN) :: IN(N)
!REAL(8),INTENT(OUT) :: OUT(N)
!    
!    CALL THIS%BARRIER()
!    call MPI_ALLREDUCE(IN,OUT,N,MPI_DOUBLE_PRECISION,MPI_SUM,THIS%MYCOM,THIS%ierr) 
!
!END SUBROUTINE
!!->+++++++++++++++++++++++++++++++++++++++++++++
!SUBROUTINE FINALIZE(THIS)
!IMPLICIT NONE
!CLASS(MPI) :: THIS
!
!    CALL MPI_COMM_FREE(THIS%MYCOM,THIS%IERR)
!    call MPI_FINALIZE(THIS%ierr)
!
!END SUBROUTINE
!!->+++++++++++++++++++++++++++++++++++++++++++++
!SUBROUTINE SEND_RECV(THIS,WS,N,WD,SOURCE,DENSTINY)
!IMPLICIT NONE
!CLASS(MPI) :: THIS
!INTEGER,INTENT(IN) :: N,SOURCE,DENSTINY
!REAL(8),INTENT(IN) :: WS(N)
!REAL(8),INTENT(OUT) :: WD(N)
!
!    THIS%tag=99
!    IF(THIS%MYID.EQ.SOURCE)    THEN
!        call MPI_SEND(WS, N, MPI_DOUBLE_PRECISION, DENSTINY, THIS%tag,THIS%MYCOM, THIS%IERR )
!    END IF
!    IF(THIS%MYID.EQ.DENSTINY)     THEN
!        call MPI_RECV( WD, N, MPI_DOUBLE_PRECISION, SOURCE, THIS%tag,THIS%MYCOM, THIS%status, THIS%IERR )
!    END IF
!
!END SUBROUTINE
!!->+++++++++++++++++++++++++++++++++++++++++++++
!SUBROUTINE ISEND_RECV(THIS,WS,N,WD,SOURCE,DENSTINY,RANK)
!IMPLICIT NONE
!CLASS(MPI) :: THIS
!INTEGER,INTENT(IN) :: N,SOURCE,DENSTINY,RANK
!REAL(8),INTENT(IN) :: WS(N)
!REAL(8),INTENT(OUT) :: WD(N)
!
!    THIS%tag=99
!    IF(THIS%MYID.EQ.SOURCE)    THEN
!        call MPI_ISEND(WS, N, MPI_DOUBLE_PRECISION, DENSTINY, THIS%tag,THIS%MYCOM, THIS%IREQS(RANK),THIS%IERR )
!    END IF
!    IF(THIS%MYID.EQ.DENSTINY)     THEN
!        call MPI_IRECV( WD, N, MPI_DOUBLE_PRECISION, SOURCE, THIS%tag,THIS%MYCOM, THIS%IREQR(RANK), THIS%IERR )
!    END IF
!
!END SUBROUTINE
!!->+++++++++++++++++++++++++++++++++++++++++++++
!SUBROUTINE SEND(THIS,WS,N,SOURCE,DENSTINY)
!IMPLICIT NONE
!CLASS(MPI) :: THIS
!INTEGER,INTENT(IN) :: N,SOURCE,DENSTINY
!REAL(8),INTENT(IN) :: WS(N)
!
!    THIS%tag=99
!    IF(THIS%MYID.EQ.SOURCE)    THEN
!        call MPI_SEND(WS, N, MPI_DOUBLE_PRECISION, DENSTINY, THIS%tag,THIS%MYCOM, THIS%IERR )
!    END IF
!
!END SUBROUTINE
!!->+++++++++++++++++++++++++++++++++++++++++++++
!SUBROUTINE SENDI(THIS,WS,SOURCE,DENSTINY)
!IMPLICIT NONE
!CLASS(MPI) :: THIS
!INTEGER,INTENT(IN) :: SOURCE,DENSTINY
!INTEGER,INTENT(IN) :: WS
!
!    THIS%tag=99
!    IF(THIS%MYID.EQ.SOURCE)    THEN
!        call MPI_SEND(WS, 1, MPI_INT, DENSTINY, THIS%tag,THIS%MYCOM, THIS%IERR )
!    END IF
!
!END SUBROUTINE
!!->+++++++++++++++++++++++++++++++++++++++++++++
!SUBROUTINE RECV(THIS,WD,N,SOURCE,DENSTINY)
!IMPLICIT NONE
!CLASS(MPI) :: THIS
!INTEGER,INTENT(IN) :: N,SOURCE,DENSTINY
!REAL(8),INTENT(OUT) :: WD(N)
!
!    THIS%tag=99
!    IF(THIS%MYID.EQ.DENSTINY)     THEN
!        call MPI_RECV( WD, N, MPI_DOUBLE_PRECISION, SOURCE, THIS%tag,THIS%MYCOM, THIS%status, THIS%IERR )
!    END IF
!
!END SUBROUTINE
!!->+++++++++++++++++++++++++++++++++++++++++++++
!SUBROUTINE RECVI(THIS,WD,SOURCE,DENSTINY)
!IMPLICIT NONE
!CLASS(MPI) :: THIS
!INTEGER,INTENT(IN) :: SOURCE,DENSTINY
!INTEGER,INTENT(OUT) :: WD
!
!    THIS%tag=99
!    IF(THIS%MYID.EQ.DENSTINY)     THEN
!        call MPI_RECV( WD, 1, MPI_INT, SOURCE, THIS%tag,THIS%MYCOM, THIS%status, THIS%IERR )
!    END IF
!
!END SUBROUTINE
!!->+++++++++++++++++++++++++++++++++++++++++++++
!SUBROUTINE ISEND(THIS,WS,N,DENSTINY,RANK)
!IMPLICIT NONE
!CLASS(MPI) :: THIS
!INTEGER,INTENT(IN) :: N,DENSTINY,RANK
!REAL(8),INTENT(IN) :: WS(N)
!
!    THIS%tag=99
!    call MPI_ISEND(WS, N, MPI_DOUBLE_PRECISION, DENSTINY, THIS%tag,THIS%MYCOM, THIS%IREQS(RANK),THIS%IERR )
!
!END SUBROUTINE
!!->+++++++++++++++++++++++++++++++++++++++++++++
!SUBROUTINE IRECV(THIS,WD,N,SOURCE,RANK)
!IMPLICIT NONE
!CLASS(MPI) :: THIS
!INTEGER,INTENT(IN) :: N,SOURCE,RANK
!REAL(8),INTENT(OUT) :: WD(N)
!
!    THIS%tag=99
!    WRITE(*,*)
!    call MPI_IRECV( WD, N, MPI_DOUBLE_PRECISION, SOURCE, THIS%tag,THIS%MYCOM, THIS%IREQR(RANK),THIS%IERR )
!END SUBROUTINE
!!->+++++++++++++++++++++++++++++++++++++++++++++
!SUBROUTINE BARRIER(THIS)
!IMPLICIT NONE
!CLASS(MPI) :: THIS
!
!        CALL MPI_Barrier(THIS%MYCOM, THIS%IERR)
!
!END SUBROUTINE
!!->+++++++++++++++++++++++++++++++++++++++++++++
!SUBROUTINE WAITS(THIS,RANK)
!IMPLICIT NONE
!CLASS(MPI) :: THIS
!INTEGER :: RANK
!
!        CALL MPI_WAIT(THIS%IREQS(RANK), THIS%status,THIS%IERR)
!
!END SUBROUTINE
!!->+++++++++++++++++++++++++++++++++++++++++++++
!SUBROUTINE WAITR(THIS,RANK)
!IMPLICIT NONE
!CLASS(MPI) :: THIS
!INTEGER :: RANK
!
!        CALL MPI_WAIT(THIS%IREQR(RANK), THIS%status,THIS%IERR)
!
!END SUBROUTINE
!!->+++++++++++++++++++++++++++++++++++++++++++++
!SUBROUTINE GATHER(THIS,N,SBUF,RBUF)
!IMPLICIT NONE
!CLASS(MPI) :: THIS
!INTEGER,INTENT(IN) :: N
!REAL(8) :: SBUF(N),RBUF(N)
!
!        CALL MPI_Allgather(SBUF, N, MPI_DOUBLE_PRECISION, rbuf, N, MPI_DOUBLE_PRECISION, THIS%MYCOM,THIS%IERR)
!
!END SUBROUTINE
!
!end MODULE
