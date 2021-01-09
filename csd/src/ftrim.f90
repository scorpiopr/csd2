MODULE FTRIM_CLASS
USE INPUT_CLASS
USE HTRIM_CLASS
USE FORWARD_CLASS
IMPLICIT NONE
!->+++++++++++++++++++++++++++++++++++++++++++++
TYPE,PUBLIC,EXTENDS(HTRIM) :: FTRIM

    INTEGER :: trimn
    TYPE(FORWARD) :: FOR
    TYPE(FORWARD),ALLOCATABLE :: FORTRIM(:)

    CONTAINS

    PROCEDURE,PUBLIC :: SOL

    PROCEDURE,PRIVATE :: CONSTRUCT_FTRIM
    PROCEDURE,PRIVATE :: INI_TRIM
    PROCEDURE,PRIVATE :: TRIMX
    PROCEDURE,PRIVATE :: proptrim
    PROCEDURE,PRIVATE :: fvnlin
    PROCEDURE,PRIVATE :: trimvct2
    PROCEDURE,PRIVATE :: OUTPUT

END TYPE FTRIM
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++PUBLIC++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
CONTAINS
!->+++++++++++++++++++++++++++++++++++++++++++++
subroutine SOL(THIS,INPU)    
implicit none
CLASS(FTRIM) :: THIS
TYPE(INPUT) :: INPU

    WRITE(*,*) '+++++++++++++++++++++++++++++++++++++++++++++'
    WRITE(*,*) '构建计算模型'
    CALL THIS%CONSTRUCT_FTRIM(INPU)
!->++++++++++++++++++++++++++++++++++++++++++++
    WRITE(*,*) '前飞初始化'
    CALL THIS%INI_TRIM()
!->++++++++++++++++++++++++++++++++++++++++++++
    WRITE(*,*) '计算前飞特性'
    CALL THIS%TRIMX()
!->++++++++++++++++++++++++++++++++++++++++++++
    WRITE(*,*) '分析力'
    CALL THIS%FOR%HUBFORCE()
    WRITE(*,*) THIS%FOR%HUBFC0,THIS%FOR%HUBMT0
!->++++++++++++++++++++++++++++++++++++++++++++
    WRITE(*,*) '输出前飞响应'
    CALL THIS%FOR%OUTWVF()
    CALL THIS%FOR%OUTVMU()

end subroutine SOL
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++PRIVATE++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE CONSTRUCT_FTRIM(THIS,INPU)
IMPLICIT NONE
CLASS(FTRIM) :: THIS
TYPE(INPUT) :: INPU

INTEGER :: I

    CALL THIS%CONSTRUCT_HTRIM(INPU)
    
    if( THIS%mktrim .eq. 0 ) then
        THIS%trimn=4
    else if( THIS%mktrim .eq. 1 ) then
        THIS%trimn=3
    else if( THIS%mktrim .eq. 2 ) then
        THIS%trimn=1
    end if

    CALL THIS%FOR%CONSTRUCT_FORWARD(INPU)
    ALLOCATE(THIS%FORTRIM(THIS%TRIMN))
    DO I=1,THIS%trimn
        CALL THIS%FORTRIM(I)%CONSTRUCT_FORWARD(INPU)
    END DO

END SUBROUTINE
!->+++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE INI_TRIM(THIS)
IMPLICIT NONE
CLASS(FTRIM) THIS

real(rdt) ladtpp, bt1c, bt1s
real(rdt) th75, th1s, th1c

    if(this%isinput.eq.0) then
        ladtpp=THIS%FOR%AEROF%LAMI+THIS%FOR%AEROF%LAMD
        bt1s=0.0D0
        bt1c=0.0D0 
        THIS%FOR%AEROF%DELTA=0.01
        THIS%FOR%AEROF%ct=THIS%cw/COS(THIS%FOR%AEROF%delta) 
        THIS%FOR%AEROF%ctc=2*THIS%FOR%AEROF%ct
        THIS%FOR%AEROF%alp=THIS%FOR%AEROF%thetapa+THIS%FOR%AEROF%delta
        CALL THIS%FOR%UPDATEINDUCEVF(THIS%FOR%OMG,THIS%FOR%BLADE%R)
!        write(*,*) th75, th1s, th1c, mu, ct, thtw, ladtpp, bt1c, bt1s, sgslid, alcsrf, lockn
        CALL THIS%trimest(th75, th1s, th1c, THIS%FOR%AEROF%mu, THIS%FOR%AEROF%ct, &
                                        ladtpp, bt1c, bt1s, THIS%FOR%AEROF%sgslid, THIS%FOR%AEROF%alcsrf)
        THIS%FOR%TH0=TH75 
        THIS%FOR%TH1S=TH1S
        THIS%FOR%TH1C=TH1C  
    end if

    CALL THIS%FOR%fwflginit()

END SUBROUTINE
!->+++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE TRIMX(THIS)
IMPLICIT NONE
CLASS(FTRIM) :: THIS

REAL(RDT), ALLOCATABLE:: UQS(:),UQSX(:)
REAL(RDT) NUQ,NUQ2
REAL(RDT) :: trimv(THIS%TRIMN),xvtm(THIS%TRIMN),xtrimv(THIS%TRIMN)
INTEGER I, II,IIRD
INTEGER :: is_converge
REAL(RDT), EXTERNAL:: DNRM2

    ALLOCATE(UQS(THIS%FOR%BLADE%NDOFC),UQSX(THIS%FOR%BLADE%NDOFC))
	
	open(102,file='trim_res.dat')
	CLOSE(102,STATUS="DELETE")
	!open(103,file='cont_res.dat')
	!CLOSE(103,STATUS="DELETE")
	OPEN(102,FILE='trim_res.dat')
	OPEN(103,FILE='cont_res.dat',POSITION='APPEND')
	
    do iird=1, THIS%FOR%nrd 


        DO I=1,THIS%TRIMN
                CALL THIS%FORTRIM(I)%COPY2(THIS%FOR)
        END DO

        UQSX=THIS%FOR%UQSA(:, 1)
        call THIS%FOR%sovlerswrap()
        UQS=THIS%FOR%UQSA(:, 1)
        NUQ=DNRM2(THIS%FOR%BLADE%NDOFC,UQS-UQSX,1)
                
        CALL  THIS%FOR%HUBFORCE()

        call THIS%trimvct2(trimv, THIS%FOR%AEROF%mu, THIS%FOR%AEROF%alp, &
                                        THIS%FOR%AEROF%delta, THIS%FOR%hubfc0, THIS%FOR%hubmt0)   

            NUQ2=ABS(THIS%FOR%AEROF%CT-THIS%FOR%HUBFC0(3))
                
            THIS%FOR%AEROF%ct=THIS%FOR%HUBFC0(3) 
            THIS%FOR%AEROF%CTC=2.0D0*THIS%FOR%AEROF%ct

            TRIM1:SELECT CASE(THIS%mktrim) 
            case(0)
                xvtm(1)=THIS%FOR%TH0
                xvtm(2)=THIS%FOR%TH1C
                xvtm(3)=THIS%FOR%TH1S
                xvtm(4)=THIS%FOR%AEROF%delta
            case(1)
                xvtm(1)=THIS%FOR%TH0
                xvtm(2)=THIS%FOR%TH1C
                xvtm(3)=THIS%FOR%TH1S
            case(2)
                xvtm(1)=THIS%FOR%TH0
            end SELECT TRIM1

            xtrimv(1:THIS%trimn)=trimv(1:THIS%trimn)

            call THIS%proptrim(xvtm, xtrimv)

            TRIM2:SELECT CASE(THIS%mktrim) 
            case(0)
                THIS%FOR%TH0=xvtm(1)
                THIS%FOR%TH1C=xvtm(2)
                THIS%FOR%TH1S=xvtm(3)
                THIS%FOR%AEROF%delta=xvtm(4)
            case(1)
                THIS%FOR%TH0=xvtm(1)
                THIS%FOR%TH1C=xvtm(2)
                THIS%FOR%TH1S=xvtm(3)
            case(2)
                THIS%FOR%TH0=xvtm(1)
            end SELECT TRIM2

            THIS%FOR%AEROF%alp=THIS%FOR%AEROF%thetapa+THIS%FOR%AEROF%delta
                
            !if( NUQ2 .lt. 1.0D-7) then
            !    is_converge=1
            !    EXIT 
            !end if
            CALL THIS%OUTPUT(IIRD,NUQ,NUQ2)
            
            CALL THIS%FOR%UPDATEINDUCEVF(THIS%FOR%OMG,THIS%FOR%BLADE%R)
    END DO 
	close(102)
	close(103)
END SUBROUTINE
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
subroutine proptrim(THIS,xv, fxv) 
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
implicit none
CLASS(FTRIM) :: THIS
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
real(rdt), intent(inout):: xv(THIS%TRIMN)
real(rdt), intent(in):: fxv(THIS%TRIMN)

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->        
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
integer err
integer i
real(rdt)   dlxi
real(rdt), allocatable:: jacb(:, :), xvtp(:), fvtp(:)

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> main body:



      ALLOCATE( jacb(THIS%TRIMN, THIS%TRIMN), xvtp(THIS%TRIMN), fvtp(THIS%TRIMN), &
            STAT=ERR)
                  
	IF(ERR.NE.0) THEN
		WRITE(*, 998)
		STOP
	END IF
	
	
      do 5 i=1, THIS%TRIMN
            xvtp=xv 
            dlxi=THIS%rlxv*xv(i) 
            if( ABS(dlxi) .lt. 1.0d-5 ) then
                  jacb(:,i)=0.0D0
            else
                xvtp(i)=xvtp(i)+dlxi
                CALL THIS%fvnlin(THIS%FORTRIM(I),fvtp, xvtp)  
                 
                jacb(:, i)=(fvtp-fxv)/dlxi 
!                if(i.eq.1) write(*,*) jacb(:,1)
!                pause
            end if
            !> write(*, *) dlxi, jacb(:, i) 
            !>   dlxi = 0 ??? 
5    continue

      call newtonrap_sp(xv, THIS%TRIMN, fxv, jacb, THIS%rlxr)

      DEALLOCATE( jacb, xvtp, fvtp, &
            STAT=ERR)
                  
	IF(ERR.NE.0) THEN
		WRITE(*, 999)
		STOP
	END IF

998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
end subroutine proptrim
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
subroutine fvnlin(THIS,FORX,fvtp, xvtp)  
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
implicit none
CLASS(FTRIM) :: THIS
CLASS(FORWARD) :: FORX
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
real(rdt), intent(in):: xvtp(THIS%trimn)
real(rdt), intent(out):: fvtp(THIS%trimn)

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->        
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
integer ii
real(rdt) :: trimv(THIS%trimn)

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> main body:
            
      if( THIS%mktrim .eq. 0 ) then
            FORX%TH0=xvtp(1)
            FORX%TH1C=xvtp(2)
            FORX%TH1S=xvtp(3)
            FORX%AEROF%delta=xvtp(4)
      else if( THIS%mktrim .eq. 1 ) then
            FORX%TH0=xvtp(1)
            FORX%TH1C=xvtp(2)
            FORX%TH1S=xvtp(3)
      else if( THIS%mktrim .eq. 2 ) then
            FORX%TH0=xvtp(1)
      end if

    FORX%AEROF%alp=FORX%AEROF%thetapa+FORX%AEROF%delta
!    CALL FORX%UPDATEINDUCEVF(FORX%OMG,FORX%BLADE%R)

      call FORX%sovlerswrap()
     CALL FORX%HUBFORCE()

      call THIS%trimvct2(trimv, FORX%AEROF%mu, FORX%AEROF%alp, FORX%AEROF%delta, FORX%hubfc0, FORX%hubmt0)   

      fvtp(1:THIS%TRIMN)=trimv(1:THIS%TRIMN)
      
end subroutine fvnlin
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
subroutine trimvct2(THIS,trimv, mu, alp, delta, fchbv, mthbv)
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
implicit none
CLASS(FTRIM) :: THIS
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
real(RDT), intent(in):: mu, alp, delta
real(RDT), intent(in):: fchbv(3), mthbv(3)
real(RDT), intent(out):: trimv(4)

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->       MVMU2, TXYZ, VCP
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
real(rdt) rcgnd(3), racnd(3), fcvnd(3), mtvnd(3), fchbtgnd(3), mthbtgnd(3)
real(rdt) fvcgnd(3), fvacnd(3)
real(rdt) fvcgnr(3), fvacnr(3)
real(rdt) fcnr(3), mtnr(3)
real(rdt) mtacnr(3), mtcgnr(3)
real(rdt) tnfcg(3, 3), tnfac(3, 3)
real(rdt) cdf,ZERO

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> main body:
    ZERO=0.0
      !>rcgnd=(/-rcg(1), rcg(2), -rcg(3)/)/rr
      !>racnd=(/rac(1), rac(2), -rac(3)/)/rr
      
      !>rcgnd=rcg/rr
      !>racnd=rac/rr 
      !>fcvnd=fchbv/trf
      !>mtvnd=mthbv/trf/rr 
     
      fchbtgnd=THIS%fchbtg/THIS%FOR%AEROF%trf
      mthbtgnd=THIS%mthbtg/THIS%FOR%AEROF%mrf
      
      rcgnd=THIS%rcgb
      racnd=THIS%racb
      fcvnd=fchbv
      mtvnd=mthbv

!      fchbtgnd=fchbtg
!      mthbtgnd=mthbtg
      
      fvcgnd=(/ZERO, ZERO, -this%cw/)
      cdf=0.5d0*THIS%fcdfpdka*(mu/COS(alp))**2.0D0!公式推导论文式6.49
      fvacnd=(/cdf, ZERO, ZERO/)   

      CALL TXYZ(tnfcg, 2, -delta)
      CALL TXYZ(tnfac, 2, -alp) 
      CALL MVMU2(fvcgnr,3,3,tnfcg,fvcgnd) 
      CALL MVMU2(fvacnr,3,3,tnfac,fvacnd) 

      fcnr=fcvnd+fvcgnr+fvacnr
      CALL VCP(mtcgnr,rcgnd,fvcgnr) 
      CALL VCP(mtacnr,racnd,fvacnr) 
      mtnr=mtvnd+mtcgnr+mtacnr!mtcgnr：公式推导论文式6.43Cw项；mtacnr：公式推导论文式6.43Cd项
    
      fcnr=fchbtgnd-fcnr 
      mtnr=mthbtgnd-mtnr
     

      trimv(1)=fcnr(3)
            !> vertical force.    
      trimv(2)=mtnr(2)
            !> pitch moment.          
      trimv(3)=mtnr(1)
            !> roll moment.  
      trimv(4)=fcnr(1)
            !> longitudinal force.   
end subroutine trimvct2


SUBROUTINE OUTPUT(THIS,IIRD,NUQ,NUQ2)
IMPLICIT NONE
CLASS(FTRIM) :: THIS
INTEGER,INTENT(IN) :: IIRD
REAL(RDT),INTENT(IN) :: NUQ,NUQ2
            write(*, *) iird 
            WRITE(*,*) 'RMS1:',NUQ
            WRITE(*,*) 'RMS2:',NUQ2
            write(*, 998) 'CT','TH0','TH1C','TH1S','THETAPA','DELTA','LAMIL'
            write(*, 997) THIS%FOR%AEROF%ct,&
                                THIS%FOR%th0,&
                                THIS%FOR%th1c,&
                                THIS%FOR%th1s,&
                                THIS%FOR%AEROF%thetapa,&
                                THIS%FOR%AEROF%delta,&
                                THIS%FOR%AEROF%LAMIL
			write(*, 997) THIS%FOR%AEROF%ct,&
                                THIS%FOR%th0/PI*180,&
                                THIS%FOR%th1c/PI*180,&
                                THIS%FOR%th1s/PI*180,&
                                THIS%FOR%AEROF%thetapa/PI*180,&
                                THIS%FOR%AEROF%delta/PI*180,&
                                THIS%FOR%AEROF%LAMIL
			WRITE(102,*) iird,NUQ,NUQ2
            
            if ( iird==THIS%FOR%nrd ) then
                write(103, '(I3,10f10.6)') iird,THIS%FOR%AEROF%ct,&
                                THIS%FOR%th0,&
                                THIS%FOR%th1c,&
                                THIS%FOR%th1s,&
                                THIS%FOR%AEROF%thetapa,&
                                THIS%FOR%AEROF%delta,&
                                THIS%FOR%AEROF%LAMIL						
            end if
								
								
997 format(10f10.6)
998 format(10a10)
END SUBROUTINE
END MODULE
