MODULE SOLVER_CLASS
USE EIGENBEAM_CLASS
USE AERO_CLASS
IMPLICIT NONE
!->+++++++++++++++++++++++++++++++++++++++++++++
TYPE,PUBLIC :: SOLVER

    REAL(RDT) :: BETA,GAM
    INTEGER :: sel, nit
    INTEGER :: mrksl

    CONTAINS

    PROCEDURE,PUBLIC :: CONSTRUCT_SOLVER
    PROCEDURE,PUBLIC :: SOLVERS

    PROCEDURE,PRIVATE :: NewmarkOneStep
    PROCEDURE,PRIVATE :: NewmarkOneStep2
    PROCEDURE,PRIVATE :: odedrivers
    PROCEDURE,PRIVATE ::  statespace
    PROCEDURE,PRIVATE ::  rk41
    PROCEDURE,PRIVATE :: RK42
    PROCEDURE,PRIVATE ::  rk43
    PROCEDURE,PRIVATE :: RK44
    PROCEDURE,PRIVATE ::  rkg
    PROCEDURE,PRIVATE :: rkFehlberg
    PROCEDURE,PRIVATE :: ptint
    PROCEDURE,PRIVATE :: mckf

END TYPE SOLVER
CONTAINS
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++PUBLIC++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE CONSTRUCT_SOLVER(THIS,INPU)    
implicit none
CLASS(SOLVER) :: THIS
TYPE(INPUT) :: INPU

    THIS%BETA=INPU%FWD%BETA
    THIS%GAM=INPU%FWD%GAM
    THIS%sel=INPU%FWD%sel
    THIS%nit=INPU%FWD%nit
    THIS%mrksl=INPU%FWD%mrksl

END SUBROUTINE

SUBROUTINE SOLVERS(THIS,BLADE,AEROF,TH0,TH1C,TH1S,OMG,SOL,III,TEND,TSDT,HH)
IMPLICIT NONE 
CLASS(SOLVER) :: THIS
CLASS(EIGENBEAM) :: BLADE
CLASS(AERO) :: AEROF
INTEGER,INTENT(IN) :: SOL
REAL(RDT),INTENT(IN) :: TH0,TH1C,TH1S,OMG

INTEGER,INTENT(IN) :: III
REAL(RDT),INTENT(IN) :: TEND,TSDT,HH
REAL(RDT),EXTERNAL :: DNRM2

real(rdt) res

         GOTO (911, 912, 913, 914) THIS%mrksl!三种求解方法选择	     
                
            
911    call THIS%odedrivers(BLADE,AEROF,&
                            TH0,TH1C,TH1S,OMG,SOL,BLADE%y, &
                            BLADE%nhf2,TSDT,HH,BLADE%YST,&
                            res, 1) 

          call THIS%statespace(BLADE,AEROF,&
                            TH0,TH1C,TH1S,OMG,SOL,tend, BLADE%y, BLADE%yp)
          write(*,*) iii,BLADE%y(1)
          goto 960     
                

912      BLADE%q2tp=BLADE%yp(BLADE%nhf+1:BLADE%nhf2)
            BLADE%q1tp=BLADE%y(BLADE%nhf+1:BLADE%nhf2)
            BLADE%qp=BLADE%y(1:BLADE%nhf)
            call THIS%NewmarkOneStep(BLADE,AEROF,tSDT,HH,&
            TH0,TH1C,TH1S,OMG,SOL)

            BLADE%y(1:BLADE%nhf)=BLADE%q
            BLADE%y(BLADE%nhf+1:BLADE%nhf2)=BLADE%q1t  
            BLADE%yp(1:BLADE%nhf)=BLADE%q1t
            BLADE%yp(BLADE%nhf+1:BLADE%nhf2)=BLADE%q2t   
            goto 960
          
            !->/// nit=1
913       continue      
            BLADE%q2tp=BLADE%yp(BLADE%nhf+1:BLADE%nhf2)
            BLADE%q1tp=BLADE%y(BLADE%nhf+1:BLADE%nhf2)
            BLADE%qp=BLADE%y(1:BLADE%nhf)             
            
            call THIS%NewmarkOneStep2(BLADE,AEROF,tSDT,HH,&
            TH0,TH1C,TH1S,OMG,SOL)

            BLADE%y(1:BLADE%nhf)=BLADE%q
            BLADE%y(BLADE%nhf+1:BLADE%nhf2)=BLADE%q1t  
            BLADE%yp(1:BLADE%nhf)=BLADE%q1t
            BLADE%yp(BLADE%nhf+1:BLADE%nhf2)=BLADE%q2t 
!          write(*,*) iii,THIS%BLADE%y(1) 
            goto 960
            
914       CONTINUE
            call THIS%ptint(BLADE,AEROF,&
                        TH0,TH1C,TH1S,OMG,SOL,BLADE%y, &
                        BLADE%nhf2,tsdt, hh,BLADE%YST)
!          write(*,*) iii,BLADE%y(1)
            call THIS%statespace(BLADE,AEROF,&
                            TH0,TH1C,TH1S,OMG,SOL,tend, BLADE%y, BLADE%yp)
            goto 960
            
960       continue
        !IF(DNRM2(BLADE%NDOFC,BLADE%y/BLADE%R,1).GE.1E2) THEN
        !    WRITE(*,*) BLADE%y/BLADE%R
        !    PAUSE
        !END IF

END SUBROUTINE
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++PRIVATE++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
subroutine ptint(THIS,BLADE,AEROF,TH0,TH1C,TH1S,OMG,SOL,yout,N,T,h,yinn)
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:	
!use Maccpti
implicit none
CLASS(SOLVER) :: THIS
CLASS(EIGENBEAM) :: BLADE
CLASS(AERO) :: AEROF
INTEGER,INTENT(IN) :: SOL
REAL(RDT),INTENT(IN) :: TH0,TH1C,TH1S,OMG

INTEGER,INTENT(IN ) :: N
real(rdt),INTENT(IN)  :: yinn(BLADE%nhf2)
real(rdt),INTENT(OUT) :: yout(BLADE%nhf2)
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
real(rdt),intent(in) :: h,t
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->        
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
real(rdt),allocatable :: ha(:,:),fa(:),ptim(:,:),ptic(:,:),ptik(:,:),ptif(:),ta(:,:)
real(rdt),allocatable :: yt(:)
integer :: col
INTEGER err
REAL(RDT) :: PPSI


!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

    col=0

    ALLOCATE(ha(BLADE%nhf2,BLADE%nhf2),fa(BLADE%nhf2),&
              ptim(BLADE%nhf,BLADE%nhf),ptic(BLADE%nhf,BLADE%nhf),&
              ptik(BLADE%nhf,BLADE%nhf),ptif(BLADE%nhf),&
              ta(BLADE%nhf2,BLADE%nhf2),yt(BLADE%nhf2),&
            STAT=ERR)

	IF(ERR.NE.0) THEN
		WRITE(*,998)
		STOP
	END IF  


      if( BLADE%natmod .eq. 1 ) then 
            goto 131
      end if 
      
     BLADE%UQS=yinn(1:BLADE%nhf)
      BLADE%UQS1T=yinn(BLADE%nhf+1:BLADE%nhf2)

      !>   ppsi=tt*ACOS(-ONED)/180  
      ppsi=omg*t
    CALL BLADE%UPDATE()
    CALL AEROF%MLEGGAUSS_AERO(OMG,BLADE,TH0,TH1C,TH1S,PPSI,SOL)
    CALL BLADE%MLEGGAUSS_ST_BEAM(TH0,TH1C,TH1S,PPSI,OMG,SOL)
    CALL BLADE%FORMULATE_ST()

      BLADE%ucmc=BLADE%ucmc+&
                                        BLADE%rdcm*BLADE%ummc+&
                                        BLADE%rdck*BLADE%ukmc

      ptim=BLADE%UMMC
      ptic=BLADE%UCMC
      ptik=BLADE%UKMC
      ptif=-BLADE%UFVC
      
      call calpticoe(ha,fa,ptim,ptic,ptik,ptif,BLADE%nhf,col)
      call calptitran(ta,ha,BLADE%nhf2,h,20)
      
!      call ptix2p(yinn,ptim,ptic,nhf,col)
      call ptisol(yout,yt,BLADE%nhf2,ta,ha,fa,yinn,h)
!        call accpti(ta2pai,y2pai,ta,yt,neqd)
!      call ptip2x(yout,ptim,ptic,nhf,col)
      
      goto 200
      
131 continue

      CALL MVMU2(BLADE%UQS, BLADE%NDOFC, BLADE%NMN, BLADE%VTEG, yinn)
      CALL MVMU2(BLADE%UQS1T, BLADE%NDOFC, BLADE%NMN, BLADE%VTEG, yinn(BLADE%NMN+1))

      !>   ppsi=tt*ACOS(-ONED)/180  
      ppsi=omg*t
    CALL BLADE%UPDATE()
    CALL AEROF%MLEGGAUSS_AERO(OMG,BLADE,TH0,TH1C,TH1S,PPSI,SOL)
    CALL BLADE%MLEGGAUSS_ST_BEAM(TH0,TH1C,TH1S,PPSI,OMG,SOL)
    CALL BLADE%FORMULATE_ST()
      
      BLADE%ucmc=BLADE%ucmc+&
                                        BLADE%rdcm*BLADE%ummc+&
                                        BLADE%rdck*BLADE%ukmc

      CALL MATMLYTP(ptim, BLADE%NMN, BLADE%NDOFC, BLADE%NDOFC, BLADE%NMN, BLADE%VTEGT, BLADE%UMMC, BLADE%VTEG)
      CALL MATMLYTP(ptic, BLADE%NMN, BLADE%NDOFC, BLADE%NDOFC, BLADE%NMN, BLADE%VTEGT, BLADE%UCMC, BLADE%VTEG)
      CALL MATMLYTP(ptik, BLADE%NMN, BLADE%NDOFC, BLADE%NDOFC, BLADE%NMN, BLADE%VTEGT, BLADE%UKMC, BLADE%VTEG)    
      CALL MVMU2(ptif, BLADE%NMN, BLADE%NDOFC, BLADE%VTEGT, -BLADE%UFVC) 
      
      

      call calpticoe(ha,fa,ptim,ptic,ptik,ptif,BLADE%NMN,col)
      call calptitran(ta,ha,BLADE%nhf2,h,20)
      
!      call ptix2p(yinn,ptim,ptic,NMND,col)
      call ptisol(yout,yt,BLADE%nhf2,ta,ha,fa,yinn,h)
!        call accpti(ta2pai,y2pai,ta,yt,neqd)
!      call ptip2x(yout,ptim,ptic,NMND,col)
      
      
!      tp1(1:NMND)=yout(NMND+1:nhf2)
!
!      CALL MVMU2(UQS2T, UQDC, NMND, VTEG, tp1(NMND+1)) 
      
   
 
      
200 continue   

            BLADE%yst=BLADE%y

	DEALLOCATE(ha,fa,&
              ptim,ptic,ptik,ptif,&
              ta,&
            STAT=ERR)

	IF(ERR.NE.0) THEN
		WRITE(*,999)
		STOP
	END IF     
	 
998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
end subroutine
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine odedrivers(THIS,BLADE,AEROF,TH0,TH1C,TH1S,OMG,SOL,y, N,x, hh, yst,res,mdsl)
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
implicit none
CLASS(SOLVER) :: THIS
CLASS(EIGENBEAM) :: BLADE
CLASS(AERO) :: AEROF
INTEGER,INTENT(IN) :: SOL
REAL(RDT),INTENT(IN) :: TH0,TH1C,TH1S,OMG
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
integer, intent(in):: N,mdsl
real(RDT), intent(in):: x, hh
real(RDT), intent(inout):: y(n),yst(n),res
REAL(RDT) :: PPSI
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->       
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
integer iflag

integer INFO(15)
integer nvx(5),nroot,nmm,err
integer LIW1, LIW2, LIW3, LIW4
integer LRW1, LRW2, LRW3, LRW4, LRW5

integer IDID, MSTATE, MINT, IERFLG

integer IPAR(1) 
real(RDT) RPAR(1)

real(RDT) RTOL(1), ATOL(1)
real(RDT) relerr, abserr

real(rdt) xend

external JACX

integer, allocatable:: iwork(:)
REAL(RDT), ALLOCATABLE:: workx(:)

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> main body:


        GOTO (911, 912, 913, 914, 915, 916, 917,918,919) mdsl	  


!        	rk41(y, n, x, hh, yst)	

911       call THIS%rk41(BLADE,AEROF,&
                            TH0,TH1C,TH1S,OMG,SOL,y,N, &
                            X, hh, yst)	
            yst=y   
          goto 960
            
912       call THIS%rk42(BLADE,AEROF,TH0,TH1C,TH1S,OMG,SOL,Y, n, x, HH, yST)
          yst=y
          goto 960 	  

913     call THIS%rk43(BLADE,AEROF,TH0,TH1C,TH1S,OMG,SOL,yst, n, x, HH)          
          y=yst
          goto 960

914       call THIS%rk44(BLADE,AEROF,TH0,TH1C,TH1S,OMG,SOL,yst, n, x, HH, res)      
          y=yst
          goto 960
          
915     call THIS%rkg(BLADE,AEROF,TH0,TH1C,TH1S,OMG,SOL,y, n, x, hh, yst, 0)
          yst=y
          goto 960 
          
916    call THIS%rkFehlberg(BLADE,AEROF,TH0,TH1C,TH1S,OMG,SOL,y,  n, x, hh, yst, res) 
          yst=y   
          goto 960
          
917       continue
918       continue
919       continue    
      xend=x+hh

      relerr=1.0D-4
      abserr=1.0D-4
      iflag=1     
        nroot=0      
        nvx(1)=7*N+33
        nvx(2)=22*N+100
        nvx(3)=21*N+130
        nvx(4)=N*N+10*N+250
        nvx(5)=N*N+17*N+250+2*NROOT
        nmm=nvx(5)
      
    ALLOCATE( workx(nmm), iwork(56+n), &      
                  STAT=ERR)
	IF(ERR.NE.0) THEN
		WRITE(*, 998)
		STOP
	END IF
	

      goto (801, 802, 803) mdsl-6


801 continue  
      RTOL=1.0D-6
      ATOL=1.0D-4
      INFO(1) = 0 
      INFO(2) = 0 
      INFO(3) = 0      
      LIW3=34 
      LRW3=7*N+33

!      call DDERKF(THIS%statespace, n, x, yst, xend, INFO, RTOL, ATOL, IDID, &
!             WORKX, LRW3, IWORK, LIW3, RPAR, IPAR)      
    STOP
      y=yst  
      goto 890
     

802 continue 
      RTOL=1.0D-3
      ATOL=1.0D-3
      INFO(1) = 0 
      INFO(2) = 0 
      INFO(3) = 0
      INFO(4) = 0  
      LIW1=51  
      LRW1=21*N+130    
        
!      call DDEABM(THIS%statespace, n, x, yst, xend, INFO, RTOL, ATOL, IDID, &
!             WORKX, LRW1, IWORK, LIW1, RPAR, IPAR)
STOP
      y=yst
      goto 890


803 continue  
      RTOL=1.0D-5
      ATOL=1.0D-5
      INFO(1) = 0 
      INFO(2) = 0 
      INFO(3) = 0
      INFO(4) = 0  
      INFO(5) = 0
      INFO(6) = 0 
      LIW2=N+56
      LRW2=N*N+10*N+250     
      
!      call DDEBDF(THIS%statespace, n, x, yst, xend, INFO, RTOL, ATOL, IDID, &
!             WORKX, LRW2, IWORK, LIW2, RPAR, IPAR, JACX)      
             STOP
      y=yst
      goto 890
      

890 continue      
 

      DEALLOCATE(workx, iwork, &      
                  STAT=ERR)
                  
	IF(ERR.NE.0) THEN
		WRITE(*, 999)
		STOP
	END IF
	
960      continue

990 FORMAT (11X, 'T', 14X, 'Y(1)', 11X, 'Y(2)', 11X, 'Y(3)')
997 FORMAT (I7, 5F15.5) 
998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')    
end SUBROUTINE

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> time integration by Newmark-beta method, one/single step.
!-> Mmt*q2t+Cmt*q1t+Kmt*q=fv
subroutine NewmarkOneStep2(THIS,BLADE,AEROF,tSDT,HH,TH0,TH1C,TH1S,OMG,SOL)
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
implicit none
CLASS(SOLVER) :: THIS
CLASS(EIGENBEAM) :: BLADE
CLASS(AERO) :: AEROF
real(RDT), intent(in):: tSDT,HH
INTEGER,INTENT(IN) :: SOL
REAL(RDT),INTENT(IN) :: TH0,TH1C,TH1S,OMG
REAL(RDT) :: PPSI
!-> LOCAL VARIABLES:
integer, parameter:: ione=1
integer i, j
real(RDT), parameter:: oned=1.0D0
real(RDT) c1, c2, c3, c4, c5, c6
REAL(RDT) :: T,DT

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> main body:      
            
        T=TSDT
        DT=HH
      c2=oned/THIS%beta/dt   
      c1=c2/dt
      c2=-c2
      c3=oned-0.5d0/THIS%beta
      c4=THIS%gam/THIS%beta
      c5=oned-c4
      c6=(oned-0.5d0*c4)*dt
      c4=c4/dt     
      
      BLADE%q2t=BLADE%q2tp 
      BLADE%q1t=BLADE%q1tp 
      BLADE%q=BLADE%qp   
      goto(111, 112) THIS%sel 

111 continue
      call THIS%mckf(BLADE,AEROF,t+dt,TH0,TH1C,TH1S,OMG,SOL) 
      BLADE%kmt=BLADE%kmt+c4*BLADE%cmt+c1*BLADE%mmt
      call mtxinv2(BLADE%kmt, 1, BLADE%nhf, BLADE%kmt)         
      do 50 i=1, THIS%nit   
            BLADE%fvx=BLADE%fv

            call DGEMV('N', BLADE%nhf, BLADE%nhf, -oned, BLADE%cmt, BLADE%nhf, -c4*BLADE%q+c5*BLADE%q1t+c6*BLADE%q2t, ione, oned, BLADE%fvx, ione)
            call DGEMV('N', BLADE%nhf, BLADE%nhf, -oned, BLADE%mmt, BLADE%nhf, -c1*BLADE%q+c2*BLADE%q1t+c3*BLADE%q2t, ione, oned, BLADE%fvx, ione)          
 
            call MVMU2(BLADE%q, BLADE%nhf, BLADE%nhf, BLADE%kmt, BLADE%fvx) 
               
            BLADE%q1t=c4*(BLADE%q-BLADE%qp)+c5*BLADE%q1tp+c6*BLADE%q2tp            
            BLADE%q2t=c1*(BLADE%q-BLADE%qp)+c2*BLADE%q1tp+c3*BLADE%q2tp   
                    
50  continue
      return

112 continue
      do 100 i=1, THIS%nit   
            call THIS%mckf(BLADE,AEROF,t+dt,TH0,TH1C,TH1S,OMG,SOL) 
            BLADE%kmt=BLADE%kmt+c4*BLADE%cmt+c1*BLADE%mmt

            call DGEMV('N', BLADE%nhf, BLADE%nhf, -oned, BLADE%cmt, BLADE%nhf, -c4*BLADE%q+c5*BLADE%q1t+c6*BLADE%q2t, ione, oned, BLADE%fv, ione)
            call DGEMV('N', BLADE%nhf, BLADE%nhf, -oned, BLADE%mmt, BLADE%nhf, -c1*BLADE%q+c2*BLADE%q1t+c3*BLADE%q2t, ione, oned, BLADE%fv, ione)   

            call linsolver2(BLADE%q, 1, BLADE%nhf, BLADE%kmt, BLADE%fv) 
            
            BLADE%q1t=c4*(BLADE%q-BLADE%qp)+c5*BLADE%q1tp+c6*BLADE%q2tp            
            BLADE%q2t=c1*(BLADE%q-BLADE%qp)+c2*BLADE%q1tp+c3*BLADE%q2tp  
            
100 continue
      return
end subroutine NewmarkOneStep2
!->+++++++++++++++++++++++++++++++++++++++++++++

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> time integration by Newmark-beta method, one/single step.
!-> Mmt*q2t+Cmt*q1t+Kmt*q=fv
subroutine NewmarkOneStep(THIS,BLADE,AEROF,tSDT,HH,TH0,TH1C,TH1S,OMG,SOL)
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
implicit none
CLASS(SOLVER) :: THIS
CLASS(EIGENBEAM) :: BLADE
CLASS(AERO) :: AEROF
real(RDT), intent(in):: tSDT,HH
INTEGER,INTENT(IN) :: SOL
REAL(RDT),INTENT(IN) :: TH0,TH1C,TH1S,OMG

integer, parameter:: ione=1
integer i, j
real(RDT), parameter:: oned=1.0D0
real(RDT) c1, c2, c3, c4
REAL(RDT) :: T,DT
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> main body:      
            
        T=TSDT
        DT=HH 
            
      c1=(0.5d0-THIS%beta)*dt*dt
      c2=THIS%beta*dt*dt
      c3=(oned-THIS%gam)*dt
      c4=THIS%gam*dt
     
      BLADE%q2t=BLADE%q2tp 
      BLADE%q1t=BLADE%q1tp 
      BLADE%q=BLADE%qp              
     
      call THIS%mckf(BLADE,AEROF,t+dt,TH0,TH1C,TH1S,OMG,SOL) 
      goto(111, 112) THIS%sel 

111 continue
      do 50 i=1, THIS%nit   
            BLADE%q=BLADE%qp+dt*BLADE%q1tp+c1*BLADE%q2tp
            BLADE%q1t=BLADE%q1tp+c3*BLADE%q2tp
            
            call DGEMV('N', BLADE%nhf, BLADE%nhf, -oned, BLADE%cmt, BLADE%nhf, BLADE%q1t, ione, oned, BLADE%fv, ione)
            call DGEMV('N', BLADE%nhf, BLADE%nhf, -oned, BLADE%kmt, BLADE%nhf, BLADE%q, ione, oned, BLADE%fv, ione)
            call mtxinv2(BLADE%mmt, 1, BLADE%nhf, BLADE%mmt+c2*BLADE%kmt+c4*BLADE%cmt)  
            call MVMU2(BLADE%q2t, BLADE%nhf, BLADE%nhf, BLADE%mmt, BLADE%fv) 
            BLADE%q=BLADE%qp+dt*BLADE%q1tp+c1*BLADE%q2tp+c2*BLADE%q2t  
            BLADE%q1t=BLADE%q1tp+c3*BLADE%q2tp+c4*BLADE%q2t
 
50  continue
      return

112 continue
      do 100 i=1, THIS%nit   
            BLADE%q=BLADE%qp+dt*BLADE%q1tp+c1*BLADE%q2tp
            BLADE%q1t=BLADE%q1tp+c3*BLADE%q2tp
            
            call DGEMV('N', BLADE%nhf, BLADE%nhf, -oned, BLADE%cmt, BLADE%nhf, BLADE%q1t, ione, oned, BLADE%fv, ione)
            call DGEMV('N', BLADE%nhf, BLADE%nhf, -oned, BLADE%kmt, BLADE%nhf, BLADE%q, ione, oned, BLADE%fv, ione)
            call linsolver2(BLADE%q2t, 1, BLADE%nhf, BLADE%mmt+c2*BLADE%kmt+c4*BLADE%cmt, BLADE%fv)  
            BLADE%q=BLADE%qp+dt*BLADE%q1tp+c1*BLADE%q2tp+c2*BLADE%q2t  
            BLADE%q1t=BLADE%q1tp+c3*BLADE%q2tp+c4*BLADE%q2t
100 continue
      return
end subroutine NewmarkOneStep
!->+++++++++++++++++++++++++++++++++++++++++++++







!->+++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> get M, C, K, and f for the newmark method.
subroutine mckf(THIS,BLADE,AEROF,t,TH0,TH1C,TH1S,OMG,SOL)
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
implicit none
CLASS(SOLVER) :: THIS
CLASS(EIGENBEAM) :: BLADE
CLASS(AERO) :: AEROF
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
real(RDT), intent(in):: t
INTEGER,INTENT(IN) :: SOL
REAL(RDT),INTENT(IN) :: TH0,TH1C,TH1S,OMG
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->       MVMU2, MATMLYTP, FORFLTMATX          
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
real(RDT), parameter:: oned=1.0D0
REAL(RDT) :: PPSI
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> main body:      

      if( BLADE%natmod .eq. 1 ) then 
            goto 131
      end if 
     
      BLADE%UQS=BLADE%q
      BLADE%UQS1T=BLADE%q1t
      BLADE%UQS2T=BLADE%q2t
      
      ppsi=omg*t
    CALL BLADE%UPDATE()
    CALL AEROF%MLEGGAUSS_AERO(OMG,BLADE,TH0,TH1C,TH1S,PPSI,SOL)
    CALL BLADE%MLEGGAUSS_ST_BEAM(TH0,TH1C,TH1S,PPSI,OMG,SOL)
    CALL BLADE%FORMULATE_ST()

      BLADE%ucmc=BLADE%ucmc+&
                                        BLADE%rdcm*BLADE%ummc+&
                                        BLADE%rdck*BLADE%ukmc
      
      BLADE%mmt=BLADE%UMMC
      BLADE%cmt=BLADE%UCMC
      BLADE%kmt=BLADE%UKMC
      BLADE%fv=-BLADE%UFVC
      return
      
131 continue
      CALL MVMU2(BLADE%UQS, BLADE%NDOFC, BLADE%NMN, BLADE%VTEG, BLADE%q)
      CALL MVMU2(BLADE%UQS1T, BLADE%NDOFC, BLADE%NMN, BLADE%VTEG, BLADE%q1t)
      CALL MVMU2(BLADE%UQS2T, BLADE%NDOFC, BLADE%NMN, BLADE%VTEG, BLADE%q2t)

      ppsi=omg*t
    CALL BLADE%UPDATE()
    CALL AEROF%MLEGGAUSS_AERO(OMG,BLADE,TH0,TH1C,TH1S,PPSI,SOL)
    CALL BLADE%MLEGGAUSS_ST_BEAM(TH0,TH1C,TH1S,PPSI,OMG,SOL)
    CALL BLADE%FORMULATE_ST()

      BLADE%ucmc=BLADE%ucmc+&
                                        BLADE%rdcm*BLADE%ummc+&
                                        BLADE%rdck*BLADE%ukmc
       
       CALL BLADE%GETMKCFT()
      return
           
end subroutine mckf
!->+++++++++++++++++++++++++++++++++++++++++++++

subroutine statespace(THIS,BLADE,AEROF,TH0,TH1C,TH1S,OMG,SOL,TT, YVS, YVSPM)
implicit none
CLASS(SOLVER) :: THIS
CLASS(EIGENBEAM) :: BLADE
CLASS(AERO) :: AEROF

INTEGER,INTENT(IN) :: SOL
REAL(RDT),INTENT(IN) :: TH0,TH1C,TH1S,OMG

REAL(RDT),INTENT(IN):: TT
REAL(RDT),INTENT(IN):: YVS(BLADE%nhf2)
REAL(RDT),INTENT(OUT):: YVSPM(BLADE%nhf2)
REAL(RDT) :: PPSI

real(RDT), parameter:: ONED=1.0D0

REAL(RDT) ALPHAT, BETAT
INTEGER INCX, INCY


      if( BLADE%natmod .eq. 1 ) then 
            goto 131
      end if 
      
     BLADE%UQS=YVS(1:BLADE%NHF)
      BLADE%UQS1T=YVS(BLADE%NHF+1:BLADE%NHF2)

      !>   ppsi=tt*ACOS(-ONED)/180  
      ppsi=omg*tt
    CALL BLADE%UPDATE()
    CALL AEROF%MLEGGAUSS_AERO(OMG,BLADE,TH0,TH1C,TH1S,PPSI,SOL)
    CALL BLADE%MLEGGAUSS_ST_BEAM(TH0,TH1C,TH1S,PPSI,OMG,SOL)
    CALL BLADE%FORMULATE_ST()

      BLADE%ucmc=BLADE%ucmc+&
                                        BLADE%rdcm*BLADE%ummc+&
                                        BLADE%rdck*BLADE%ukmc

      ALPHAT=1.0D0
      BETAT=1.0D0
      INCX=1
      INCY=1  
      call DGEMV('N',BLADE%NDOFC,BLADE%NDOFC,-ALPHAT,BLADE%UCMC,BLADE%NDOFC,BLADE%UQS1T,INCX,-BETAT,BLADE%UFVC,INCY)
      call DGEMV('N',BLADE%NDOFC,BLADE%NDOFC,-ALPHAT,BLADE%UKMC,BLADE%NDOFC,BLADE%UQS,INCX,BETAT,BLADE%UFVC,INCY)  
      
      call MVMU2(BLADE%UQS2T,BLADE%NDOFC,BLADE%NDOFC,BLADE%MMTINV,BLADE%UFVC)       
      
      YVSPM(1:BLADE%NHF)=BLADE%UQS1T    
      YVSPM(BLADE%NHF+1:BLADE%nhf2)=BLADE%UQS2T
      
      goto 200
      
131 continue

      CALL MVMU2(BLADE%UQS, BLADE%NDOFC, BLADE%NMN, BLADE%VTEG, YVS)
      CALL MVMU2(BLADE%UQS1T, BLADE%NDOFC, BLADE%NMN, BLADE%VTEG, YVS(BLADE%NMN+1))

      !>   ppsi=tt*ACOS(-ONED)/180  
      PPSI=omg*tt
    CALL BLADE%UPDATE()
    CALL AEROF%MLEGGAUSS_AERO(OMG,BLADE,TH0,TH1C,TH1S,PPSI,SOL)
    CALL BLADE%MLEGGAUSS_ST_BEAM(TH0,TH1C,TH1S,PPSI,OMG,SOL)
    CALL BLADE%FORMULATE_ST()

      BLADE%ucmc=BLADE%ucmc+&
                                        BLADE%rdcm*BLADE%ummc+&
                                        BLADE%rdck*BLADE%ukmc
            
      ALPHAT=1.0D0
      BETAT=1.0D0
      INCX=1
      INCY=1  
      call DGEMV('N',BLADE%NDOFC,BLADE%NDOFC,-ALPHAT,BLADE%UCMC,BLADE%NDOFC,BLADE%UQS1T,INCX,-BETAT,BLADE%UFVC,INCY)
      call DGEMV('N',BLADE%NDOFC,BLADE%NDOFC,-ALPHAT,BLADE%UKMC,BLADE%NDOFC,BLADE%UQS,INCX,BETAT,BLADE%UFVC,INCY)  
  
      CALL MVMU2(YVSPM(BLADE%NMN+1), BLADE%NMN, BLADE%NDOFC, BLADE%VTEGT, BLADE%UFVC) 

      CALL MVMU2(YVSPM(BLADE%NMN+1),BLADE%NMN,BLADE%NMN,BLADE%MMTINV,YVSPM(BLADE%NMN+1)) 
       
      YVSPM(1:BLADE%NMN)=YVS(BLADE%NMN+1:BLADE%nhf2)

      CALL MVMU2(BLADE%UQS2T, BLADE%NDOFC, BLADE%NMN, BLADE%VTEG, YVSPM(BLADE%NMN+1)) 
 
      
200 continue   
      
      
end subroutine statespace


!* Numerical Recipes subroutine for fourth order Runge-Kutta integration:
!*
!* Given values for the variables y(1:n) and their derivatives dydx(i:n) 
!* known at x, use the fourth-order Runge-Kutta method to advance the 
!* solution over an interval h and return the incremented variables as 
!* yout(1:n), which need not be a distinct array from y. The user supplies 
!* the subroutine derivs(x,y,dydx), which returns derivatives dydx at x.
!*
subroutine rk41(THIS,BLADE,AEROF,TH0,TH1C,TH1S,OMG,SOL,yout, n, x, h, y)
implicit none
CLASS(SOLVER) :: THIS
CLASS(EIGENBEAM) :: BLADE
CLASS(AERO) :: AEROF

INTEGER,INTENT(IN) :: SOL
REAL(RDT),INTENT(IN) :: TH0,TH1C,TH1S,OMG

integer, intent(in):: n
real(RDT), intent(in):: x, h
real(RDT), intent(in):: y(n)
real(RDT), intent(out):: yout(n)

integer i
integer err
REAL(RDT) h6, hh, xh
REAL(RDT), ALLOCATABLE:: dym(:), dyt(:), yt(:), dydx(:) 

      ALLOCATE(dym(n), dyt(n), yt(n), dydx(n), &
            STAT=ERR)

      IF(ERR.NE.0) THEN
	      WRITE(*, 998)
	      STOP
      END IF

    call THIS%statespace(BLADE,AEROF,&
                    TH0,TH1C,TH1S,OMG,SOL,x, y, dydx)

      hh=h*0.5D0
      h6=h/6.0D0
      xh=x+hh

      do i=1,n
            yt(i)=y(i)+hh*dydx(i)
      enddo
    call THIS%statespace(BLADE,AEROF,&
                    TH0,TH1C,TH1S,OMG,SOL,xh, yt, dyt)
      do i=1,n
            yt(i)=y(i)+hh*dyt(i)
      enddo
    call THIS%statespace(BLADE,AEROF,&
                    TH0,TH1C,TH1S,OMG,SOL,xh, yt, dym)
      do i=1,n
            yt(i)=y(i)+h*dym(i)
            dym(i)=dyt(i)+dym(i)
      enddo
    call THIS%statespace(BLADE,AEROF,&
                    TH0,TH1C,TH1S,OMG,SOL,x+h, yt, dyt)
      do i=1,n
            yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2d0*dym(i))
      enddo

      DEALLOCATE(dym, dyt, yt, dydx, &
            STAT=ERR)

      IF(ERR.NE.0) THEN
	      WRITE(*, 999)
	      STOP
      END IF

998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')      
end subroutine rk41


subroutine RK42(THIS,BLADE,AEROF,TH0,TH1C,TH1S,OMG,SOL,yout, n, x, hx, y)
implicit none
CLASS(SOLVER) :: THIS
CLASS(EIGENBEAM) :: BLADE
CLASS(AERO) :: AEROF

INTEGER,INTENT(IN) :: SOL
REAL(RDT),INTENT(IN) :: TH0,TH1C,TH1S,OMG

integer, intent(in):: n
real(RDT), intent(in):: x, hx
real(RDT), intent(in):: y(n)
real(RDT), intent(out):: yout(n)

integer, parameter:: rkstp=4
real(RDT) one 
real(RDT) cc(rkstp), bb(rkstp), aa(rkstp, rkstp)

integer err
integer i, j
REAL(RDT) xx
REAL(RDT), ALLOCATABLE:: ykk(:), ypkt(:), ypkk(:, :)

      ALLOCATE(ykk(n), ypkt(n), ypkk(n, rkstp), &
            STAT=ERR)

      IF(ERR.NE.0) THEN
	      WRITE(*, 998)
	      STOP
      END IF
      
      one=1.0D0     
      bb=one*(/1.0, 2.0, 2.0, 1.0/)/6.0  
      cc=one*(/0.0, 0.5, 0.5, 1.0/)    
      aa=0.0D0
      aa(2, 1)=one*0.5
      aa(3, 2)=one*0.5
      aa(4, 3)=one                                    
      
      ypkk=0
      yout=y 
      do 10 i=1, rkstp
            xx=x+cc(i)*hx
            ykk=y
            do 20 j=1, rkstp-1
                  ykk=ykk+ypkk(:, j)*aa(i, j)  
20        continue      
        call THIS%statespace(BLADE,AEROF,&
                    TH0,TH1C,TH1S,OMG,SOL,xx, ykk, ypkt)
            ypkk(:, i)=hx*ypkt(:)     
            yout=yout+bb(i)*ypkk(:, i)                                          
10  continue         

      DEALLOCATE(ykk, ypkt, ypkk, &
            STAT=ERR)

      IF(ERR.NE.0) THEN
	      WRITE(*, 998)
	      STOP
      END IF

998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')          
end subroutine RK42



SUBROUTINE RK43(THIS,BLADE,AEROF,TH0,TH1C,TH1S,OMG,SOL,X, n, T, h)
implicit none
CLASS(SOLVER) :: THIS
CLASS(EIGENBEAM) :: BLADE
CLASS(AERO) :: AEROF

INTEGER,INTENT(IN) :: SOL
REAL(RDT),INTENT(IN) :: TH0,TH1C,TH1S,OMG

integer, intent(in):: n
real(RDT), intent(in):: T, H
real(RDT), intent(inout):: X(n)

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->      
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
integer err
REAL(RDT) h2, tt
REAL(RDT), ALLOCATABLE:: xx(:), dydx(:), f1(:), f2(:), f3(:), f4(:)

      ALLOCATE(xx(n), dydx(n), f1(n), f2(n), f3(n), f4(n), &
            STAT=ERR)

      IF(ERR.NE.0) THEN
	      WRITE(*, 998)
	      STOP
      END IF
      
      H2 = 0.5*H
  
      tt=t
      xx=x   
    call THIS%statespace(BLADE,AEROF,&
                    TH0,TH1C,TH1S,OMG,SOL,tt, xx, dydx)
      F1 = H*dydx    

      tt=T + H2
      xx=X + 0.5*F1 
    call THIS%statespace(BLADE,AEROF,&
                    TH0,TH1C,TH1S,OMG,SOL,tt, xx, dydx)
      F2 = H*dydx       
 
      tt=T + H2
      xx=X + 0.5*F2   
    call THIS%statespace(BLADE,AEROF,&
                    TH0,TH1C,TH1S,OMG,SOL,tt, xx, dydx)

      F3 = H*dydx  

      tt=T + H
      xx=X + F3 
    call THIS%statespace(BLADE,AEROF,&
                    TH0,TH1C,TH1S,OMG,SOL,tt, xx, dydx)
      F4 = H*dydx       
       
      X = X + (F1 + F2 + F2  + F3 + F3 + F4)/6.0
       
      DEALLOCATE(xx, dydx, f1, f2, f3, f4, &
            STAT=ERR)

      IF(ERR.NE.0) THEN
	      WRITE(*, 998)
	      STOP
      END IF

998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')       
END SUBROUTINE







!C
!C RUNGE-KUTTA-FEHLBERG METHOD FOR SOLVING AN INITIAL VALUE PROBLEM (RK45,F)  
!C      
SUBROUTINE RK44(THIS,BLADE,AEROF,TH0,TH1C,TH1S,OMG,SOL,X, n, T, h,EST)
IMPLICIT NONE
CLASS(SOLVER) :: THIS
CLASS(EIGENBEAM) :: BLADE
CLASS(AERO) :: AEROF

INTEGER,INTENT(IN) :: SOL
REAL(RDT),INTENT(IN) :: TH0,TH1C,TH1S,OMG

integer, intent(in):: n
real(RDT), intent(in):: T, H
real(RDT), intent(inout):: X(n)
real(RDT), intent(out):: EST
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
REAL(RDT), EXTERNAL:: DNRM2
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->      
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES.
integer err
REAL(RDT) C21,C31,C32, C41,C42,C43, C51,C52,C53,C54
REAL(RDT) C61,C62,C63,C64,C65, A1,A3,A4,A5, B1,B3,B4,B5, B6, C40
DATA C21,C31,C32, C41,C42,C43, C51,C52,C53,C54, &
            C61,C62,C63,C64,C65, A1,A3,A4,A5, B1,B3,B4,B5, B6, C40 &
            /0.25,0.09375,0.28125, &
            0.87938097405553,-3.2771961766045,3.3208921256258, &
            2.0324074074074,-8.0,7.1734892787524,-0.20589668615984, &  
            -0.2962962962963,2.0,-1.3816764132554,0.45297270955166,-0.275, &
            0.11574074074074,0.54892787524366,0.5353313840156,-0.2, &
            0.11851851851852,0.51898635477583,0.50613149034201,-0.18, &
            0.036363636363636, 0.92307692307692/   
REAL(RDT) tt
REAL(RDT), ALLOCATABLE:: f1(:), f2(:), f3(:), f4(:), f5(:), f6(:)
REAL(RDT), ALLOCATABLE:: xx(:), x5(:), dydx(:)

      ALLOCATE(f1(n), f2(n), f3(n), f4(n), f5(n), f6(n), &
            xx(n), x5(n), dydx(n), & 
            STAT=ERR)

      IF(ERR.NE.0) THEN
	      WRITE(*, 998)
	      STOP
      END IF
      
      tt=t
      xx=x   
    call THIS%statespace(BLADE,AEROF,&
                    TH0,TH1C,TH1S,OMG,SOL,tt, xx, dydx)
      F1 = H*dydx 
      
      tt=T+ 0.25*H
      xx=X + C21*F1  
    call THIS%statespace(BLADE,AEROF,&
                    TH0,TH1C,TH1S,OMG,SOL,tt, xx, dydx)
      F2 = H*dydx 
      
      tt=T+0.375*H
      xx=X + C31*F1 + C32*F2
    call THIS%statespace(BLADE,AEROF,&
                    TH0,TH1C,TH1S,OMG,SOL,tt, xx, dydx)
      F3 = H*dydx 
      
      tt=T+C40*H
      xx=X + C41*F1 + C42*F2 + C43*F3
    call THIS%statespace(BLADE,AEROF,&
                    TH0,TH1C,TH1S,OMG,SOL,tt, xx, dydx)
      F4 = H*dydx 
                                        
      tt=T+H 
      xx=X + C51*F1 + C52*F2 + C53*F3 + C54*F4
    call THIS%statespace(BLADE,AEROF,&
                    TH0,TH1C,TH1S,OMG,SOL,tt, xx, dydx)
      F5 = H*dydx 
  
      tt=T+0.5*H
      xx=X + C61*F1 + C62*F2 + C63*F3 + C64*F4 + C65*F5
    call THIS%statespace(BLADE,AEROF,&
                    TH0,TH1C,TH1S,OMG,SOL,tt, xx, dydx)
      F6 = H*dydx 
        
      X5 = X + B1*F1 + B3*F3 + B4*F4 + B5*F5 + B6*F6      
      X  = X + A1*F1 + A3*F3 + A4*F4 + A5*F5    
      EST = DNRM2(n, X - X5, 1)    


      DEALLOCATE(f1, f2, f3, f4, f5, f6, &
            xx, x5, dydx, & 
            STAT=ERR)

      IF(ERR.NE.0) THEN
	      WRITE(*, 998)
	      STOP
      END IF
      
998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')    
END  SUBROUTINE


!->+++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> First-Order ode solver, dy/dt=f(y,t), explicit Runge-Kutta methods, one step.
!-> ch=1: runge-kutta-gill method.
!-> ch=0: classical runge-kutta method.
!->       coded by XiaoYu, 05162009.
!->       edited by rotor, 05192009.
subroutine rkg(THIS,BLADE,AEROF,TH0,TH1C,TH1S,OMG,SOL,yout, n, x, h, y, ch)
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
IMPLICIT NONE
CLASS(SOLVER) :: THIS
CLASS(EIGENBEAM) :: BLADE
CLASS(AERO) :: AEROF

INTEGER,INTENT(IN) :: SOL
REAL(RDT),INTENT(IN) :: TH0,TH1C,TH1S,OMG

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
integer, intent(in):: n,ch
real(RDT), intent(in):: x,h,y(n)
real(RDT), intent(out):: yout(n)

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->      
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
integer i
integer err
REAL(RDT) sqrttwo
REAL(RDT) w(4),c(4),a(4,4)
REAL(RDT), ALLOCATABLE:: yp(:,:), yt(:)

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

      ALLOCATE(yp(n,4), yt(n), &
            STAT=ERR)

      IF(ERR.NE.0) THEN
            WRITE(*, 998)
            STOP
      END IF

      c=(/0.0D0, 1.0D0, 1.0D0, 2.0D0/)
      c(:)=c(:)/2.0D0
                
      if(ch.eq.1) then
            sqrttwo=2.0D0
            sqrttwo=SQRT(sqrttwo)
            w(1)=1.0D0
            w(2)=2.0D0-sqrttwo
            w(3)=2.0D0+sqrttwo
            w(4)=1.0D0
            w(:)=w(:)/6.0D0

            a=0.0D0
            a(2,1)=1.0D0
            a(3,1)=sqrttwo-1.0D0
            a(3,2)=2.0D0-sqrttwo
            a(4,2)=-sqrttwo
            a(4,3)=2.0D0+sqrttwo
            a(:, :)=a(:, :)/2.0D0
      else if(ch.eq.0) then
            w=(/1.0D0, 2.0D0, 2.0D0, 1.0D0/)
            w(:)=w(:)/6.0D0

            a=0.0D0
            a(2,1)=1.0D0
            a(3,2)=1.0D0
            a(4,3)=2.0D0
            a(:, :)=a(:, :)/2.0D0
      else
            write(*,*) "error"
      end if
    
      yp=0.0D0
      yt=y

      do 10 i=1,4
            call THIS%statespace(BLADE,AEROF,&
                            TH0,TH1C,TH1S,OMG,SOL,x+h*c(i), yt, yp(:,i))
            yp(:,i)=yp(:,i)*h
            if(i.ne.4) then
                  yt(:)=y(:)+a(i+1,1)*yp(:,1)+a(i+1,2)*yp(:,2)+a(i+1,3)*yp(:,3)
            else if(i.eq.4) then
                  yout(:)=y(:)+w(1)*yp(:,1)+w(2)*yp(:,2)+w(3)*yp(:,3)+w(4)*yp(:,4)
            end if
10  continue
 

      DEALLOCATE(yp, yt, &
            STAT=ERR)

      IF(ERR.NE.0) THEN
            WRITE(*, 999)
            STOP
      END IF

998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')      
end subroutine rkg
!->+++++++++++++++++++++++++++++++++++++++++++++

subroutine rkFehlberg(THIS,BLADE,AEROF,TH0,TH1C,TH1S,OMG,SOL,yout,  n, x, hx, y, res)
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
IMPLICIT NONE
CLASS(SOLVER) :: THIS
CLASS(EIGENBEAM) :: BLADE
CLASS(AERO) :: AEROF

INTEGER,INTENT(IN) :: SOL
REAL(RDT),INTENT(IN) :: TH0,TH1C,TH1S,OMG

integer, intent(in):: n
real(RDT), intent(in):: x, hx
real(RDT), intent(in):: y(n)
real(RDT), intent(out):: yout(n)
real(RDT), intent(out):: res

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
REAL(RDT), EXTERNAL:: DNRM2
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->      
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES.
real(rdt) one
integer, parameter:: rkstp=6
real(rdt) cc(rkstp), bb(rkstp), bh(rkstp), aa(rkstp, rkstp)

integer err
integer i, j
REAL(RDT) xx
REAL(RDT), ALLOCATABLE:: ykk(:), ypkt(:), yhgh(:), ypkk(:, :)

      ALLOCATE(ykk(n), ypkt(n), yhgh(n), ypkk(n, rkstp), &
            STAT=ERR)

      IF(ERR.NE.0) THEN
	      WRITE(*, 998)
	      STOP
      END IF
      
      one=1.0D0 
      cc=one*(/0.0, 1.0/4.0,  3.0/8.0, 12.0/13.0, 1.0, 1.0/2.0/)                                     
      bb=one*(/25.0/216.0, 0.0, 1408.0/2565.0, 2197.0/4104.0, -1.0/5.0, 0.0/) 
      bh=one*(/16.0/135.0, 0.0, 6656.0/12825.0, 28561.0/56430.0, -9.0/50.0, 2.0/55.0/)     
      aa(:, 1)=one*(/0.0, 1.0/4.0, 3.0/32.0, 1932.0/2197.0, 439.0/216.0, -8.0/27.0/)
      aa(:, 2)=one*(/0.0, 0.0, 9.0/32.0, -7200.0/2197.0, -8.0, 2.0/)
      aa(:, 3)=one*(/0.0, 0.0, 0.0, 7296.0/2197.0, 3680.0/513.0, -3544.0/2565.0/)
      aa(:, 4)=one*(/0.0, 0.0, 0.0, 0.0, -845.0/4104.0, 1859.0/4104.0/)
      aa(:, 5)=one*(/0.0, 0.0, 0.0, 0.0, 0.0, -11.0/40.0/)
      aa(:, 6)=one*(/0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)    
     
                
      ypkk=0
      yout=y 
      yhgh=y 
      do 10 i=1, rkstp
            xx=x+cc(i)*hx
            ykk=y
            do 20 j=1, rkstp-1
                  ykk=ykk+ypkk(:, j)*aa(i, j)  
20        continue      
            call THIS%statespace(BLADE,AEROF,&
                            TH0,TH1C,TH1S,OMG,SOL,xx, ykk, ypkt)
            ypkk(:, i)=hx*ypkt(:)     
            yout=yout+bb(i)*ypkk(:, i)  
            yhgh=yhgh+bh(i)*ypkk(:, i)                                         
10  continue         
     
      res=DNRM2(n, yout-yhgh, 1)

      DEALLOCATE(ykk, ypkt, yhgh, ypkk, &
            STAT=ERR)

      IF(ERR.NE.0) THEN
	      WRITE(*, 998)
	      STOP
      END IF      

998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')         
end subroutine rkFehlberg

END MODULE
