MODULE LB234_CLASS
USE GlobalDataFun
IMPLICIT NONE
TYPE,PUBLIC :: LB2

    real(rdt),allocatable :: afastdt(:),qstdt(:),cncfstdt(:),fpnstdt(:),fpmstdt(:)
    real(rdt),allocatable :: xct(:),yct(:),zct(:)
    real(rdt),allocatable :: d1t(:),d2t(:),d3t(:),d4t(:),d5t(:),d6t(:),d7nt(:),d7mt(:)
    real(rdt),allocatable :: delta_xx1stdt(:),delta_xx2stdt(:),delta_xx3stdt(:),delta_xx4stdt(:)
    real(rdt),allocatable :: cnvt(:),taovt(:)
    real(rdt),allocatable :: cvstdt(:)

    CONTAINS

    PROCEDURE,PUBLIC :: LBINI2
    PROCEDURE,PUBLIC :: LBDE1

end TYPE LB2


TYPE,PUBLIC :: LB3

    real(rdt),allocatable :: afastdtx(:),qstdtx(:)
    real(rdt),allocatable :: xctx(:),yctx(:),zctx(:)
    real(rdt),allocatable :: d1tx(:),d2tx(:),d3tx(:),d4tx(:),d5tx(:)
    real(rdt),allocatable :: delta_xx1stdtx(:),delta_xx2stdtx(:),delta_xx3stdtx(:),delta_xx4stdtx(:)
    real(rdt),allocatable :: cncfstdtx(:),fpnstdtx(:),fpmstdtx(:)
    real(rdt),allocatable :: d6tx(:),d7ntx(:),d7mtx(:)
    real(rdt),allocatable :: cnvtx(:),taovtx(:)
    real(rdt),allocatable :: cvstdtx(:)
    
    real(rdt),allocatable :: afastdtxx(:),qstdtxx(:)
    real(rdt),allocatable :: xctxx(:),yctxx(:),zctxx(:)
    real(rdt),allocatable :: d1txx(:),d2txx(:),d3txx(:),d4txx(:),d5txx(:)
    real(rdt),allocatable :: delta_xx1stdtxx(:),delta_xx2stdtxx(:),delta_xx3stdtxx(:),delta_xx4stdtxx(:)
    real(rdt),allocatable :: cncfstdtxx(:),fpnstdtxx(:),fpmstdtxx(:)
    real(rdt),allocatable :: d6txx(:),d7ntxx(:),d7mtxx(:)
    real(rdt),allocatable :: cnvtxx(:),taovtxx(:)
    real(rdt),allocatable :: cvstdtxx(:)

    CONTAINS

    PROCEDURE,PUBLIC :: LBINI3
    PROCEDURE,PUBLIC :: LBDE2

end TYPE LB3

TYPE,PUBLIC :: LB4

    real(rdt) :: afastdo,qstdo,cncfstdo,fpnstdo,fpmstdo
    real(rdt) :: xco,yco,zco
    real(rdt) :: d1o,d2o,d3o,d4o,d5o,d6o,d7no,d7mo
    real(rdt) :: delta_xx1stdo,delta_xx2stdo,delta_xx3stdo,delta_xx4stdo
    real(rdt) :: cnvo,taovo
    real(rdt) :: cvstdo

end TYPE LB4

CONTAINS 

SUBROUTINE LBINI2(THIS,NSTPT)
IMPLICIT NONE
CLASS(LB2) :: THIS
INTEGER,INTENT(IN) :: NSTPT

INTEGER :: ERR

    ALLOCATE(THIS%AFASTDT(NSTPT),THIS%QSTDT(NSTPT),THIS%XCT(NSTPT),&
        THIS%YCT(NSTPT),THIS%ZCT(NSTPT),THIS%D1T(NSTPT),&
        THIS%D2T(NSTPT),THIS%D3T(NSTPT),THIS%D4T(NSTPT),&
        THIS%D5T(NSTPT),THIS%DELTA_XX1STDT(NSTPT),THIS%DELTA_XX2STDT(NSTPT),&
        THIS%DELTA_XX3STDT(NSTPT),THIS%DELTA_XX4STDT(NSTPT),&
        THIS%CNCFSTDT(NSTPT),THIS%FPMSTDT(NSTPT),THIS%FPNSTDT(NSTPT),&
        THIS%D6T(NSTPT),THIS%D7MT(NSTPT),THIS%D7NT(NSTPT),&
        THIS%CNVT(NSTPT),THIS%TAOVT(NSTPT),THIS%CVSTDT(NSTPT),&
        STAT=ERR)
    
	IF(ERR.NE.0) THEN
		WRITE(*,998)
		STOP
	END IF

    THIS%AFASTDT=0.0D0
    THIS%QSTDT=0.0D0
    THIS%CNCFSTDT=0.0D0
    THIS%FPMSTDT=0.0D0
    THIS%FPNSTDT=0.0D0
    THIS%CVSTDT=0.0D0
    
    THIS%XCT=0.0D0
    THIS%YCT=0.0D0
    THIS%ZCT=0.0D0
    
    THIS%D1T=0.0D0
    THIS%D2T=0.0D0
    THIS%D3T=0.0D0
    THIS%D4T=0.0D0
    THIS%D5T=0.0D0
    THIS%D6T=0.0D0
    THIS%D7MT=0.0D0
    THIS%D7NT=0.0D0
    
    THIS%DELTA_XX1STDT=0.0D0
    THIS%DELTA_XX2STDT=0.0D0
    THIS%DELTA_XX3STDT=0.0D0
    THIS%DELTA_XX4STDT=0.0D0
    
    THIS%CNVT=0.0D0
    THIS%TAOVT=0.0D0

998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
END SUBROUTINE



subroutine LBINI3(THIS,nstpt)
IMPLICIT NONE
integer,intent(in) :: nstpt
CLASS(LB3) :: THIS
integer :: err

    allocate(THIS%afastdtx(nstpt),THIS%qstdtx(nstpt),THIS%xctx(nstpt),THIS%yctx(nstpt),THIS%zctx(nstpt),&
        THIS%d1tx(nstpt),THIS%d2tx(nstpt),THIS%d3tx(nstpt),THIS%d4tx(nstpt),THIS%d5tx(nstpt),&
        THIS%delta_xx1stdtx(nstpt),THIS%delta_xx2stdtx(nstpt),THIS%delta_xx3stdtx(nstpt),THIS%delta_xx4stdtx(nstpt),&
        THIS%cncfstdtx(nstpt),THIS%fpmstdtx(nstpt),THIS%fpnstdtx(nstpt),&
        THIS%d6tx(nstpt),THIS%d7mtx(nstpt),THIS%d7ntx(nstpt),&
        THIS%cnvtx(nstpt),THIS%taovtx(nstpt),THIS%cvstdtx(nstpt),&
        THIS%afastdtxx(nstpt),THIS%qstdtxx(nstpt),THIS%xctxx(nstpt),THIS%yctxx(nstpt),THIS%zctxx(nstpt),&
        THIS%d1txx(nstpt),THIS%d2txx(nstpt),THIS%d3txx(nstpt),THIS%d4txx(nstpt),THIS%d5txx(nstpt),&
        THIS%delta_xx1stdtxx(nstpt),THIS%delta_xx2stdtxx(nstpt),THIS%delta_xx3stdtxx(nstpt),THIS%delta_xx4stdtxx(nstpt),&
        THIS%cncfstdtxx(nstpt),THIS%fpnstdtxx(nstpt),THIS%fpmstdtxx(nstpt),&
        THIS%d6txx(nstpt),THIS%d7ntxx(nstpt),THIS%d7mtxx(nstpt),&
        THIS%cnvtxx(nstpt),THIS%taovtxx(nstpt),THIS%cvstdtxx(nstpt),&
        stat=err)
    
	IF(ERR.NE.0) THEN
		WRITE(*,998)
		STOP
	END IF

!*************************1****************************
    
   THIS%afastdtx=0.0D0
   THIS%qstdtx=0.0D0
   THIS%cncfstdtx=0.0D0
   THIS%fpmstdtx=0.0D0
   THIS%fpnstdtx=0.0D0
   THIS%cvstdtx=0.0D0

    
   THIS%xctx=0.0D0
   THIS%yctx=0.0D0
   THIS%zctx=0.0D0

   THIS%d1tx=0.0D0
   THIS%d2tx=0.0D0
   THIS%d3tx=0.0D0
   THIS%d4tx=0.0D0
   THIS%d5tx=0.0D0
   THIS%d6tx=0.0D0
   THIS%d7mtx=0.0D0
   THIS%d7ntx=0.0D0

   THIS%delta_xx1stdtx=0.0D0
   THIS%delta_xx2stdtx=0.0D0
   THIS%delta_xx3stdtx=0.0D0
   THIS%delta_xx4stdtx=0.0D0

   THIS%cnvtx=0.0D0
   THIS%taovtx=0.0D0
!*************************2****************************
   THIS%afastdtxx=0.0D0
   THIS%qstdtxx=0.0D0
   THIS%cncfstdtxx=0.0D0
   THIS%fpmstdtxx=0.0D0
   THIS%fpnstdtxx=0.0D0
   THIS%cvstdtxx=0.0D0
  

   THIS%xctxx=0.0D0
   THIS%yctxx=0.0D0
   THIS%zctxx=0.0D0

   THIS%d1txx=0.0D0
   THIS%d2txx=0.0D0
   THIS%d3txx=0.0D0
   THIS%d4txx=0.0D0
   THIS%d5txx=0.0D0
   THIS%d6txx=0.0D0
   THIS%d7mtxx=0.0D0
   THIS%d7ntxx=0.0D0

   THIS%delta_xx1stdtxx=0.0D0
   THIS%delta_xx2stdtxx=0.0D0
   THIS%delta_xx3stdtxx=0.0D0
   THIS%delta_xx4stdtxx=0.0D0

   THIS%cnvtxx=0.0D0
   THIS%taovtxx=0.0D0
    
   
998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
END SUBROUTINE


SUBROUTINE LBDE1(THIS)
IMPLICIT NONE
CLASS(LB2) :: THIS
INTEGER :: ERR

    DEALLOCATE(THIS%AFASTDT,THIS%QSTDT,THIS%XCT,THIS%YCT,THIS%ZCT,&
        THIS%D1T,THIS%D2T,THIS%D3T,THIS%D4T,THIS%D5T,&
        THIS%DELTA_XX1STDT,THIS%DELTA_XX2STDT,THIS%DELTA_XX3STDT,THIS%DELTA_XX4STDT,&
        THIS%CNCFSTDT,THIS%FPMSTDT,THIS%FPNSTDT,&
        THIS%D6T,THIS%D7MT,THIS%D7NT,&
        THIS%CNVT,THIS%TAOVT,THIS%CVSTDT,&
        STAT=ERR)
    
	IF(ERR.NE.0) THEN
		WRITE(*,999)
		STOP
	END IF

    
998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
END SUBROUTINE

subroutine LBDE2(THIS)
IMPLICIT NONE
CLASS(LB3) :: THIS
integer :: err

    deallocate(THIS%afastdtx,THIS%qstdtx,THIS%xctx,THIS%yctx,THIS%zctx,&
        THIS%d1tx,THIS%d2tx,THIS%d3tx,THIS%d4tx,THIS%d5tx,&
        THIS%delta_xx1stdtx,THIS%delta_xx2stdtx,THIS%delta_xx3stdtx,THIS%delta_xx4stdtx,&
        THIS%cncfstdtx,THIS%fpmstdtx,THIS%fpnstdtx,&
        THIS%d6tx,THIS%d7mtx,THIS%d7ntx,&
        THIS%cnvtx,THIS%taovtx,THIS%cvstdtx,&
        THIS%afastdtxx,THIS%qstdtxx,THIS%xctxx,THIS%yctxx,THIS%zctxx,&
        THIS%d1txx,THIS%d2txx,THIS%d3txx,THIS%d4txx,THIS%d5txx,&
        THIS%delta_xx1stdtxx,THIS%delta_xx2stdtxx,THIS%delta_xx3stdtxx,THIS%delta_xx4stdtxx,&
        THIS%cncfstdtxx,THIS%fpnstdtxx,THIS%fpmstdtxx,&
        THIS%d6txx,THIS%d7ntxx,THIS%d7mtxx,&
        THIS%cnvtxx,THIS%taovtxx,THIS%cvstdtxx,&
        stat=err)
    
	IF(ERR.NE.0) THEN
		WRITE(*,999)
		STOP
	END IF

   
   
998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
END SUBROUTINE

END MODULE