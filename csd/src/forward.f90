MODULE FORWARD_CLASS
USE EIGENBEAM_CLASS
USE AERO_CLASS
USE HTRIM_CLASS
USE SOLVER_CLASS
USE TFEM_CLASS
IMPLICIT NONE
!->+++++++++++++++++++++++++++++++++++++++++++++
TYPE,PUBLIC :: FORWARD

    TYPE(EIGENBEAM) :: BLADE
    TYPE(AERO) :: AEROF
    TYPE(HTRIM) :: MTRIM
    TYPE(SOLVER) :: SOLV
    TYPE(TFEM) :: TIMEFEM

    REAL(RDT) :: TH0,TH1C,TH1S,OMG
    REAL(RDT) :: PPSI

    INTEGER :: mktfe
    REAL(RDT) :: hdg
    INTEGER :: ndgop, nrd

    INTEGER :: NSPPR,NSTP,NHW

    INTEGER :: SOL=2

    REAL(RDT) :: HH,TSDT,TEND


    REAL(RDT),ALLOCATABLE:: UQSA(:,:),UQS1TA(:,:),UQS2TA(:,:)
    REAL(RDT),ALLOCATABLE:: FPSI(:,:,:),MPSI(:,:,:),AOAPSI(:,:)

    REAL(RDT),ALLOCATABLE:: BRFCINR(:,:),BRFCAER(:,:),BRFCR(:,:)
    REAL(RDT),ALLOCATABLE:: BRMTINR(:,:),BRMTAER(:,:),BRMTR(:,:)

    REAL(RDT), ALLOCATABLE:: HUBFCAER(:, :), HUBMTAER(:, :)
    REAL(RDT), ALLOCATABLE:: HUBFCINR(:, :), HUBMTINR(:, :)
    REAL(RDT), ALLOCATABLE:: HUBFCR(:, :), HUBMTR(:, :)

    REAL(RDT), ALLOCATABLE:: HUBFCAERHW(:, :), HUBMTAERHW(:, :)
    REAL(RDT), ALLOCATABLE:: HUBFCINRHW(:, :), HUBMTINRHW(:, :)
    REAL(RDT), ALLOCATABLE:: HUBFCRHW(:, :), HUBMTRHW(:, :)

    REAL(RDT) HUBFC0(3), HUBMT0(3)
    REAL(RDT) HUBFC1(3), HUBMT1(3),NODIMP1,NODIMP2
    
    CONTAINS

    PROCEDURE,PUBLIC :: SOLCASE
    PROCEDURE,PUBLIC :: CONSTRUCT_FORWARD
    PROCEDURE,PUBLIC ::  fwflginit
    PROCEDURE,PUBLIC :: fordriver
    PROCEDURE,PUBLIC ::  SOVLERSWRAP
    PROCEDURE,PUBLIC ::  OUTWVF
    PROCEDURE,PUBLIC :: HUBFORCE
    PROCEDURE,PRIVATE::  SETALL
    PROCEDURE,PUBLIC :: UPDATEINDUCEVF
    PROCEDURE,PUBLIC :: COPY
    PROCEDURE,PUBLIC :: fmlct
    PROCEDURE,PUBLIC :: OUTVMU
	
!	PROCEDURE,PUBLIC :: GETBPBM
!    PROCEDURE,PUBLIC :: GETBVQA
    PROCEDURE,PUBLIC :: GETBPABC
	
END TYPE FORWARD
CONTAINS
!->+++++++++++++++++++++++++++++++++++++++++++++
subroutine SOLCASE(THIS,INPU)    
implicit none
CLASS(FORWARD) :: THIS
TYPE(INPUT) :: INPU
REAL(RDT) :: DROP
INTEGER :: I,J,K
    WRITE(*,*) '+++++++++++++++++++++++++++++++++++++++++++++'
    WRITE(*,*) '构建计算模型'
    CALL THIS%CONSTRUCT_FORWARD(INPU)
!->++++++++++++++++++++++++++++++++++++++++++++
    WRITE(*,*) '前飞初始化'
    CALL THIS%fwflginit()
!->++++++++++++++++++++++++++++++++++++++++++++
    WRITE(*,*) '计算前飞特性'
    !OPEN(444,FILE='cl1_rotor.dat')
    !DO I=1,THIS%AEROF%NSTP
    !    DO J=1,NPTS
    !        READ(444,"(1000F)") DROP,DROP,(THIS%AEROF%CFDFPSI(K,J,I),K=1,3),&
	!							(THIS%AEROF%CFDMPSI(K,J,I),K=1,3)
    !            !IF(RLUDAUD(J).LT.0.5) THEN		
    !            !    THIS%AEROF%CFDFPSI(:,J,I)=0.0
    !            !    THIS%AEROF%CFDMPSI(:,J,I)=0.0			
	!			!END IF			
	!							
	!				!IF(DABS(THIS%AEROF%CFDFPSI(2,J,I)).GT.0.1) THEN
	!				!	THIS%AEROF%CFDFPSI(1:3,J,I)=0.0
	!				!END IF
	!				!IF(DABS(THIS%AEROF%CFDMPSI(1,J,I)).GT.0.1) THEN
	!				!	THIS%AEROF%CFDMPSI(1:3,J,I)=0.0
	!				!END IF					
	!							
	!							
    !    END DO
    !END DO
    !close(444)
    !THIS%AEROF%MKFM=6
    CALL THIS%fordriver()
!->++++++++++++++++++++++++++++++++++++++++++++
    WRITE(*,*) '分析力'
    CALL THIS%HUBFORCE()
    WRITE(*,*) THIS%HUBFC0,THIS%HUBMT0
!->++++++++++++++++++++++++++++++++++++++++++++
    WRITE(*,*) '输出前飞响应'
    CALL THIS%OUTWVF()
    CALL THIS%OUTVMU()
    
end subroutine SOLCASE 
!->+++++++++++++++++++++++++++++++++++++++++++++
subroutine CONSTRUCT_FORWARD(THIS,INPU)    
implicit none
CLASS(FORWARD) :: THIS
TYPE(INPUT) :: INPU

real(rdt), parameter:: oned=1.0D0

    THIS%mktfe=INPU%FWD%mktfe
    THIS%hdg=INPU%FWD%hdg
    THIS%ndgop=INPU%FWD%ndgop
    THIS%nrd=INPU%FWD%nrd

    THIS%TH0=INPU%IPT%TH0
    THIS%TH1C=INPU%IPT%TH1C
    THIS%TH1S=INPU%IPT%TH1S
    THIS%omg =INPU%MAT%OMG

    THIS%NODIMP1=INPU%GETMASS()*THIS%omg**2*INPU%MAT%R**2
    THIS%NODIMP2=INPU%GETMASS()*THIS%omg**2*INPU%MAT%R**3
    
    CALL THIS%BLADE%CONSTRUCT_EIGENBEAM( THIS%TH0, THIS%TH1C, THIS%TH1S,THIS%omg,INPU)    

    if(THIS%mktfe.eq.0) then
        CALL THIS%SOLV%CONSTRUCT_SOLVER(INPU)
        THIS%nsppr=IDNINT(360.0D0/THIS%hdg) 
    else
        CALL THIS%TIMEFEM%CONSTRUCT_TFEM(INPU)
        THIS%nsppr=THIS%TIMEFEM%npl*THIS%TIMEFEM%ntfe
    end if
          
    THIS%nstp=THIS%nsppr+1
    THIS%nhw=floor((THIS%nstp-2)/2.0) 

    CALL THIS%AEROF%CONSTRUCT_AERO(INPU,THIS%OMG,THIS%BLADE%R)
    CALL THIS%AEROF%CONSTRUCT_AERO2(THIS%NSTP,THIS%OMG,INPU)
    CALL THIS%MTRIM%CONSTRUCT_HTRIM(INPU)

    ALLOCATE(THIS%UQSA(THIS%BLADE%NDOFC,THIS%NSTP),&
                        THIS%UQS1TA(THIS%BLADE%NDOFC,THIS%NSTP),&
                        THIS%UQS2TA(THIS%BLADE%NDOFC,THIS%NSTP))

    ALLOCATE(THIS%FPSI(6,NPTS,THIS%NSTP),&
                        THIS%MPSI(6,NPTS,THIS%NSTP),&
                        THIS%AOAPSI(NPTS,THIS%NSTP))

    ALLOCATE(THIS%BRFCAER(3,THIS%NSTP),&
                        THIS%BRMTAER(3,THIS%NSTP),&
                        THIS%BRFCINR(3,THIS%NSTP),&
                        THIS%BRMTINR(3,THIS%NSTP),&
                        THIS%BRFCR(3,THIS%NSTP),&
                        THIS%BRMTR(3,THIS%NSTP))

    ALLOCATE(THIS%HUBFCAER(3,THIS%NSTP),&
                        THIS%HUBMTAER(3,THIS%NSTP),&
                        THIS%HUBFCINR(3,THIS%NSTP),&
                        THIS%HUBMTINR(3,THIS%NSTP),&
                        THIS%HUBFCR(3,THIS%NSTP),&
                        THIS%HUBMTR(3,THIS%NSTP))

    ALLOCATE(THIS%HUBFCAERHW(THIS%NSPPR,3),&
                        THIS%HUBMTAERHW(THIS%NSPPR,3),&
                        THIS%HUBFCINRHW(THIS%NSPPR,3),&
                        THIS%HUBMTINRHW(THIS%NSPPR,3),&
                        THIS%HUBFCRHW(THIS%NSPPR,3),&
                        THIS%HUBMTRHW(THIS%NSPPR,3))

    THIS%hh=THIS%hdg*ACOS(-ONED)/180.0D0/THIS%omg   

    if( this%mktfe .eq. 1) then

        call THIS%TIMEFEM%tfemini(THIS%BLADE%nhf,THIS%AEROF%nstp,THIS%omg,&
                                THIS%BLADE%q,THIS%BLADE%yst)
    end if      

end subroutine CONSTRUCT_FORWARD 
 !->+++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
subroutine fwflginit(THIS)
implicit none
CLASS(FORWARD) :: THIS

integer i,j

real(rdt) ladtpp, bt1c, bt1s


	open(111,file='wvf_tip.dat')
	CLOSE(111,STATUS="DELETE")
	open(111,file='stmpsi.dat')
	CLOSE(111,STATUS="DELETE")  
	open(111,file='fpsi.dat')
	CLOSE(111,STATUS="DELETE") 


    THIS%UQSA(:,:)=0.0D0
    THIS%UQS1TA(:,:)=0.0D0
    THIS%UQS2TA(:,:)=0.0D0
   
    THIS%FPSI=0.0D0
    THIS%MPSI=0.0D0
    THIS%AOAPSI=0.0D0

    THIS%AEROF%alp=THIS%AEROF%thetapa+THIS%AEROF%delta

    CALL THIS%UPDATEINDUCEVF(THIS%OMG,THIS%BLADE%R)

    THIS%BLADE%UQS=0.0D0
    THIS%BLADE%UQS1T=0.0D0
    THIS%BLADE%UQS2T=0.0D0
    
    THIS%PPSI=0.0D0
    
    CALL THIS%BLADE%UPDATE()
    CALL THIS%AEROF%MLEGGAUSS_AERO(THIS%OMG,THIS%BLADE,THIS%TH0,THIS%TH1C,THIS%TH1S,THIS%PPSI,THIS%SOL)
    CALL THIS%BLADE%MLEGGAUSS_ST_BEAM(THIS%TH0,THIS%TH1C,THIS%TH1S,THIS%PPSI,THIS%OMG,THIS%SOL)
    CALL THIS%BLADE%FORMULATE_ST()
    !CALL THIS%BLADE%OUTPUTMATRIX()

    CALL THIS%BLADE%INITIAL_UQSX()
    
    WRITE(*, *) 'CT:',THIS%AEROF%ct
    WRITE(*, *) 'TH0:',THIS%TH0
    WRITE(*, *) 'TH1C:',THIS%TH1C
    WRITE(*, *) 'TH1S:',THIS%TH1S
    WRITE(*, *) 'thetapa:',THIS%AEROF%thetapa
    WRITE(*, *) 'delta:',THIS%AEROF%delta
	open(141,file='run.dat')
    WRITE(141, *) 'CT:',THIS%AEROF%ct
    WRITE(141, *) 'TH0:',THIS%TH0
    WRITE(141, *) 'TH1C:',THIS%TH1C
    WRITE(141, *) 'TH1S:',THIS%TH1S
    WRITE(141, *) 'thetapa:',THIS%AEROF%thetapa
    WRITE(141, *) 'delta:',THIS%AEROF%delta
    WRITE(141, *) 'R:',THIS%BLADE%R
    WRITE(141, *) 'cchordrf:',THIS%AEROF%cchordrf
    close(141)
    write(*,*) 'vtip:',THIS%BLADE%uqs(1)/THIS%BLADE%R
    write(*,*) 'wtip:',THIS%BLADE%uqs(3)/THIS%BLADE%R
    write(*,*) 'ftip:',THIS%BLADE%uqs(5)

    !pause
    
    THIS%UQSA(:,1)=THIS%BLADE%UQS
    THIS%UQS1TA(:,1)=THIS%BLADE%UQS1T
    THIS%UQS2TA(:,1)=THIS%BLADE%UQS2T
    THIS%FPSI(:,:,1)=THIS%AEROF%FLCSTCOE
    THIS%MPSI(:,:,1)=THIS%AEROF%MLCSTCOE
    THIS%AOAPSI(:,1)=THIS%AEROF%AOALCST
            
!    if(mktrim .ne. -1 .or. mkim.eq.2) then
!        CALL BRFRMT(UQD,UQDC,SMD,UVC0,UVC2,UQDS,UNC,PPSI,SOL)  
!
!        BRFCINR(:, 1)=BRFCIN(:)/trf 
!        BRFCAER(:, 1)=BRFCAE(:)/trf 
!        BRFCR(:, 1)=BRFC(:)/trf             
!        BRMTINR(:, 1)=BRMTIN(:)/mrf
!        BRMTAER(:, 1)=BRMTAE(:)/mrf
!        BRMTR(:, 1)=BRMT(:)/mrf
!    end if
!      
!!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   
    THIS%tsdt=0.0D0
    THIS%tend=THIS%hh
!    
!    if(mkim.eq.2) then
!        ttinflow=0.0D0
!        ytinflow=0.0D0
!        do i=1,dm
!            ttinflow(i,i)=1.0D0
!        end do
!    end if
!
!!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!    WRITE(*,*) 'INITIALING OK'
!!    PAUSE 'PRESS ANY KEY TO CONTINUE'
!      return
!            
!995 format(1000f36.12)
end subroutine fwflginit
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE fordriver(THIS)
IMPLICIT NONE
CLASS(FORWARD) :: THIS
REAL(RDT), ALLOCATABLE:: UQS(:),UQSX(:)
REAL(RDT) NUQ
INTEGER IIRD
INTEGER :: is_converge
!external STATESPACE
REAL(RDT), EXTERNAL:: DNRM2

    ALLOCATE(UQS(THIS%BLADE%NDOFC),UQSX(THIS%BLADE%NDOFC))
	open(101,file='wake_rms.dat')
	
    do iird=1, THIS%nrd

        UQSX=THIS%UQSA(:, 1)
        call THIS%sovlerswrap()
        UQS=THIS%UQSA(:, 1)
        NUQ=DNRM2(THIS%BLADE%NDOFC,UQS-UQSX,1)
		
		write(*,*) iird,NUQ
		write(101,*) iird,NUQ

        if( nuq .lt. 1.0D-9.or.iird.eq.THIS%NRD) then
            exit
            is_converge=1 
        end if

        CALL THIS%UPDATEINDUCEVF(THIS%OMG,THIS%BLADE%R)
    END DO 
    
    DEALLOCATE(UQS,UQSX)
    close(101)
	
END SUBROUTINE
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine sovlerswrap(this)
implicit none
CLASS(FORWARD) :: THIS

!external statespace
integer :: ii
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
REAL(RDT), EXTERNAL:: DNRM2
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->        
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
integer iii, iij, i, j

real(rdt) tmp
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if( THIS%mktfe .eq. 1) then !使用时间有限元
        call THIS%TIMEFEM%timefmwp(THIS%BLADE,THIS%AEROF,&
                THIS%TH0,THIS%TH1C,THIS%TH1S,THIS%SOL,THIS%omg)
        call THIS%TIMEFEM%procuqstfm(THIS%BLADE,THIS%UQSA, THIS%UQS1TA,THIS%UQS2TA,THIS%AEROF%nstp)

    ELSE
    
        do 120 iii=1, THIS%nsppr 
                
            !II=II+1
     
            THIS%AEROF%FLCSTCOE=THIS%FPSI(:,:,iii+1)
            THIS%AEROF%MLCSTCOE=THIS%MPSI(:,:,iii+1)
            THIS%AEROF%AOALCST=THIS%AOAPSI(:,iii+1)

            CALL THIS%SOLV%SOLVERS(THIS%BLADE,THIS%AEROF,&
                    THIS%TH0,THIS%TH1C,THIS%TH1S,THIS%OMG,THIS%SOL,III,THIS%TEND,THIS%TSDT,THIS%HH)
            CALL THIS%SETALL(III)

            THIS%tsdt=THIS%tend
            THIS%tend=THIS%tend+THIS%hh
120     continue  

    END IF
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
995 format(400f26.12)
end subroutine sovlerswrap
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE SETALL(THIS,III)
IMPLICIT NONE
CLASS(FORWARD) :: THIS
INTEGER,INTENT(IN) :: III

    CALL THIS%BLADE%GETUQS()
            
    THIS%UQSA(:, iii+1)=THIS%BLADE%uqs
    THIS%UQS1TA(:, iii+1)=THIS%BLADE%uqs1T
    THIS%UQS2TA(:, iii+1)=THIS%BLADE%uqs2T
    THIS%FPSI(:,:,iii+1)=THIS%AEROF%FLCSTCOE
    THIS%MPSI(:,:,iii+1)=THIS%AEROF%MLCSTCOE
    THIS%AOAPSI(:,iii+1)=THIS%AEROF%AOALCST

    IF(III.EQ.THIS%nsppr) THEN
        THIS%UQSA(:, 1)=THIS%UQSA(:, THIS%NSTP)
        THIS%UQS1TA(:, 1)=THIS%UQS1TA(:, THIS%NSTP)
        THIS%UQS2TA(:, 1)=THIS%UQS2TA(:, THIS%NSTP)  
        THIS%FPSI(:,:,1)=THIS%FPSI(:,:,THIS%NSTP)
        THIS%MPSI(:,:,1)=THIS%MPSI(:,:,THIS%NSTP)
        THIS%AOAPSI(:,1)=THIS%AOAPSI(:,THIS%NSTP)
    END IF

END SUBROUTINE
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE OUTWVF(THIS)
IMPLICIT NONE
CLASS(FORWARD) :: THIS

INTEGER :: I,J
REAL(RDT) :: PPSI,tmp

INTEGER :: N
REAL(RDT),ALLOCATABLE :: RXUDA(:)
REAL(RDT) :: rxud
REAL(RDT),ALLOCATABLE :: PEX(:,:,:),QEX(:,:,:)
REAL(RDT),ALLOCATABLE :: VM(:,:),WM(:,:)
REAL(RDT),ALLOCATABLE :: UM(:,:),FM(:,:)
CHARACTER(100),EXTERNAL :: RtoChar
!!! TEMP
INTEGER :: NL,NI,K
LOGICAL :: THERE
REAL(RDT) :: SS
REAL(RDT),ALLOCATABLE :: RL(:),TOR(:)
REAL(RDT),EXTERNAL :: FWM

    ALLOCATE(VM(THIS%BLADE%NC,2*THIS%nstp),WM(THIS%BLADE%NC,2*THIS%nstp))
    ALLOCATE(UM(THIS%BLADE%NQ,THIS%nstp),FM(THIS%BLADE%NQ,THIS%nstp))
    CALL THIS%BLADE%DECOMUQ(THIS%UQSA,VM,WM,FM,UM,THIS%nstp)
    
    !open(89,FILE='input/fp.dat')
    !read(89,*) nl
    !ALLOCATE(RL(NL),TOR(NL))
    !do i=1,nl
    !    read(89,*) rl(i),tor(i)
    !end do
    !close(89)
    !rl=rl*PI/180
    !DO J=THIS%BLADE%NELE,0,-1
    !    DO I=1,THIS%nstp
    !        call rxlct2(NI, SS, nl, THIS%AEROF%PSIV(I), rl) 
    !        WM(J+1,2*I-1:2*I)=WM(J+1,2*I-1:2*I)/WM(1,2*I-1)*FWM(tor(NI),tor(NI+1),SS)
    !    END DO
    !END DO
    !DEALLOCATE(RL,TOR)
	
    !!open(89,FILE='input/lg.dat')
    !!read(89,*) nl
    !!ALLOCATE(RL(NL),TOR(NL))
    !!do i=1,nl
    !!    read(89,*) rl(i),tor(i)
    !!end do
    !!close(89)
    !!rl=rl*PI/180
    !!DO J=THIS%BLADE%NELE,0,-1
    !!    DO I=1,THIS%nstp
    !!        call rxlct2(NI, SS, nl, THIS%AEROF%PSIV(I), rl) 
    !!        VM(J+1,2*I-1:2*I)=VM(J+1,2*I-1:2*I)/VM(1,2*I-1)*FWM(tor(NI),tor(NI+1),SS)
    !!    END DO
    !!END DO
    !!DEALLOCATE(RL,TOR)
    !!!WM=0.0D0
    !!!DO J=THIS%BLADE%NELE,0,-1
    !!!    DO I=1,THIS%nstp
    !!!        WM(J+1,2*I-1:2*I)=WM(J+1,2*THIS%nstp-1:2*THIS%nstp)
    !!!    END DO
    !!!END DO
   
    !open(89,FILE='input/fm.dat')
    !read(89,*) nl
    !ALLOCATE(RL(NL),TOR(NL))
    !do i=1,nl
    !    read(89,*) rl(i),tor(i)
    !end do
    !close(89)
    !rl=rl*PI/180
    !DO J=THIS%BLADE%NQ,1,-1
    !    DO I=1,THIS%nstp
    !        call rxlct2(NI, SS, nl, THIS%AEROF%PSIV(I), rl) 
    !        FM(J,I)=FM(J,I)/FM(1,I)*FWM(tor(NI),tor(NI+1),SS)*PI/180
    !    END DO
    !END DO
	!fm=0
    
    CALL THIS%BLADE%COMUQ(THIS%UQSA,VM,WM,FM,UM,THIS%nstp)
    !CALL THIS%BLADE%DECOMUQ(THIS%UQSA,VM,WM,FM,UM,THIS%nstp)
    
	
    open(89,FILE='wvf_tip.dat',position='append')
    WRITE(89,'(3A,100(5A))') 'PSI',' flap',' lag',' torision',' flap',' lag',' torision'
    WRITE(89,*) 'deg',' NoDim',' NoDim',' NoDim',' NoDim',' NoDim',' NoDim'
    WRITE(89,'(3A,100(15A))') 'deg',' flap',' lag',' torision',' flap',' lag',' torision'
    DO I=1,THIS%nstp
        ppsi=THIS%AEROF%PSIV(I)
        tmp=THIS%th0+THIS%th1c*COS(PPSI)+THIS%th1s*SIN(PPSI)
        WRITE(89,995) PPSI*180/PI,(THIS%UQSA(1,I)*SIN(tmp)+THIS%UQSA(3,I)*COS(tmp))/THIS%BLADE%R*100,&
                                    -(THIS%UQSA(1,I)*COS(tmp)-THIS%UQSA(3,I)*SIN(tmp))/THIS%BLADE%R*100,&
                                    (THIS%UQSA(5,I))*100,THIS%UQSA(3,I)/THIS%BLADE%R*100,&
                                    THIS%UQSA(1,I)/THIS%BLADE%R*100,THIS%UQSA(5,I)*100
    END DO
    close(89)	
	
    open(89,FILE='wvf_beam.dat')

    DO J=THIS%BLADE%NELE,1,-1
        WRITE(89, 995) ((THIS%BLADE%ELEMEN(THIS%BLADE%NELE-J+1)%NOD(1)%RLUD,VM(J+1,2*I-1)/THIS%BLADE%R*100,WM(J+1,2*I-1)/THIS%BLADE%R*100,FM(2*J+1,I)/PI*180),I=1,1)
    END DO
    WRITE(89, 995) ((THIS%BLADE%ELEMEN(THIS%BLADE%NELE)%NOD(2)%RLUD,VM(1,2*I-1)/THIS%BLADE%R*100,WM(1,2*I-1)/THIS%BLADE%R*100,FM(1,I)/PI*180),I=1,1)

    close(89)
    INQUIRE( FILE='aim_slice.dat', EXIST=THERE ) 
    IF ( THERE ) THEN

        OPEN(10,FILE='aim_slice.dat')
        READ(10,*) N
        ALLOCATE(RXUDA(N),PEX(6,THIS%NSTP,N),QEX(6,THIS%NSTP,N))
        DO I=1,N
            READ(10,*) RXUDA(I)
        END DO
        CLOSE(10)

        open(135, file='fpsi.dat',position='append')
    
        DO I=1,THIS%NSTP
            DO J=1,N
                CALL THIS%fmlct(peX(:,I,J), qeX(:,I,J),RXUDA(J),I) 
            END DO
        END DO

        WRITE(135,'(3A,100(5A))') 'PSI',('  CL',I=1,N)
        WRITE(135,*) 'deg'
        WRITE(135,'(3A,100(15A))') 'deg',(' SLICE='//trim(ADJUSTL(RtoChar(RXUDA(I)))),I=1,N)
        DO I=1,THIS%NSTP
            write(135, 995) THIS%AEROF%PSIV(I)*180/PI, pex(3,I,1:N), pex(6,I,1:N)
        END DO

        CLOSE(135)   
    
        open(135, file='fpsi2.dat')
    
        DO J=1,NPTS
            write(135, 995) RLUDAUD(J), (THIS%FPSI(3,J,I),I=1,THIS%NSTP)
        END DO

        CLOSE(135)   
    
        READ(10,*) N
        IF(ALLOCATED(RXUDA)) DEALLOCATE(RXUDA,PEX,QEX)
        ALLOCATE(RXUDA(N),PEX(3,THIS%NSTP,N),QEX(3,THIS%NSTP,N))
        DO I=1,N
            READ(10,*) RXUDA(I)
        END DO
        CLOSE(10)
    
        open(135, file='stmpsi.dat',position='append')
    
        DO I=1,THIS%AEROF%nstp
            THIS%ppsi=THIS%AEROF%PSIV(I)
            THIS%BLADE%UQS=THIS%UQSA(:, i)
            THIS%BLADE%UQS1T=THIS%UQS1TA(:, i)
            THIS%BLADE%UQS2T=THIS%UQS2TA(:, i)

            CALL THIS%BLADE%UPDATE()
        
            DO J=1,N
                CALL THIS%BLADE%BDSTLD(RXUDA(J),THIS%TH0,THIS%TH1C,THIS%TH1S,THIS%PPSI,THIS%OMG,PEX(:,I,J))
            END DO
        END DO
  
        WRITE(135,'(3A,100(5A))') 'PSI',('  CL',I=1,N)
        WRITE(135,*) 'deg'
        WRITE(135,'(3A,100(15A))') 'deg',(' SLICE='//trim(ADJUSTL(RtoChar(RXUDA(I)))),I=1,N)
        DO I=1,THIS%NSTP
            write(135, 995) THIS%AEROF%PSIV(I)*180/PI, &
                        pex(1:3,I,1:N)/(THIS%AEROF%RHA*THIS%AEROF%nrbd*THIS%AEROF%cchordrf*THIS%BLADE%R**4*THIS%OMG**2)
        END DO

        CLOSE(135)   
   
        DEALLOCATE(RXUDA,PEX,QEX)
        
    END IF
   OPEN(444,FILE='cl3_rotor.dat')
    WRITE(444,*) 'ZONE I=',NPTS,',J=',THIS%AEROF%NSTP
    DO I=1,THIS%AEROF%NSTP
        DO J=1,NPTS
            WRITE(444,"(1000F)") RLUDAUD(J)*SIN(THIS%AEROF%PSIV(I)),-RLUDAUD(J)*COS(THIS%AEROF%PSIV(I)),&
                                 (THIS%FPSI(K,J,I),K=1,3), (THIS%MPSI(K,J,I),K=1,3)
        END DO
    END DO
    close(444)

995 format(400f26.12)

END SUBROUTINE
!->++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE OUTVMU(THIS)
IMPLICIT NONE
CLASS(FORWARD) :: THIS

INTEGER :: I
REAL(RDT) :: PPSI,tmp

    OPEN(122,FILE='vmu.dat',position='append')
    write(122, 995) THIS%AEROF%mu,THIS%th0,THIS%th1c,-THIS%th1s,&
                             THIS%AEROF%thetapa,THIS%AEROF%delta,&
                             THIS%AEROF%lamd+THIS%AEROF%lami
    close(122)

    OPEN(122,FILE='hubfm.dat',position='append')
    write(122, 995) THIS%AEROF%mu,THIS%HUBFC1(1:3)*THIS%AEROF%trf/THIS%NODIMP1,&
    THIS%HUBMT1(1:3)*THIS%AEROF%mrf/THIS%NODIMP2
    close(122)
    
995 format(400f26.12)
END SUBROUTINE
!->++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE HUBFORCE(THIS)
IMPLICIT NONE
CLASS(FORWARD) :: THIS

INTEGER :: I

    DO I=1,THIS%AEROF%nstp
        THIS%ppsi=THIS%AEROF%PSIV(I)
        THIS%BLADE%UQS=THIS%UQSA(:, i)
        THIS%BLADE%UQS1T=THIS%UQS1TA(:, i)
        THIS%BLADE%UQS2T=THIS%UQS2TA(:, i)

        CALL THIS%BLADE%UPDATE()
        CALL THIS%AEROF%MLEGGAUSS_AERO(THIS%OMG,THIS%BLADE,THIS%TH0,THIS%TH1C,THIS%TH1S,&
                                                                    THIS%PPSI,THIS%SOL)
        CALL THIS%BLADE%MLEGGAUSS_FR_BEAM(THIS%TH0,THIS%TH1C,THIS%TH1S,THIS%PPSI,THIS%OMG,THIS%SOL)
        CALL THIS%BLADE%FORMULATE_FR()

        THIS%BRFCINR(:, I)=THIS%BLADE%BRFCIN(:)/THIS%AEROF%trf 
        THIS%BRFCAER(:, I)=THIS%BLADE%BRFCAE(:)/THIS%AEROF%trf 
        THIS%BRFCR(:, I)=THIS%BLADE%BRFC(:)/THIS%AEROF%trf
        THIS%BRMTINR(:, I)=THIS%BLADE%BRMTIN(:)/THIS%AEROF%mrf
        THIS%BRMTAER(:, I)=THIS%BLADE%BRMTAE(:)/THIS%AEROF%mrf
        THIS%BRMTR(:, I)=THIS%BLADE%BRMT(:)/THIS%AEROF%mrf
    END DO

    call fmrbsgl2hub(THIS%HUBFCINR, THIS%HUBFCINRHW, 3, THIS%AEROF%nstp, THIS%AEROF%nrbd, THIS%nhw, THIS%BRFCINR) 
    call fmrbsgl2hub(THIS%HUBFCAER, THIS%HUBFCAERHW, 3, THIS%AEROF%nstp, THIS%AEROF%nrbd, THIS%nhw, THIS%BRFCAER) 
    call fmrbsgl2hub(THIS%HUBFCR, THIS%HUBFCRHW, 3, THIS%AEROF%nstp, THIS%AEROF%nrbd, THIS%nhw, THIS%BRFCR)

    call fmrbsgl2hub(THIS%HUBMTINR, THIS%HUBMTINRHW, 3, THIS%AEROF%nstp, THIS%AEROF%nrbd, THIS%nhw, THIS%BRMTINR) 
    call fmrbsgl2hub(THIS%HUBMTAER, THIS%HUBMTAERHW, 3, THIS%AEROF%nstp, THIS%AEROF%nrbd, THIS%nhw, THIS%BRMTAER) 
    call fmrbsgl2hub(THIS%HUBMTR, THIS%HUBMTRHW, 3, THIS%AEROF%nstp, THIS%AEROF%nrbd, THIS%nhw, THIS%BRMTR)

    THIS%HUBFC0(:)=THIS%HUBFCRHW(1, :)
    THIS%HUBMT0(:)=THIS%HUBMTRHW(1, :)
    
    DO I=1,3
        THIS%HUBFC1(:)=sqrt(THIS%HUBFCRHW(THIS%AEROF%nrbd*2, :)**2+THIS%HUBFCRHW(THIS%AEROF%nrbd*2+1, :)**2)
        THIS%HUBMT1(:)=sqrt(THIS%HUBMTRHW(THIS%AEROF%nrbd*2, :)**2+THIS%HUBMTRHW(THIS%AEROF%nrbd*2+1, :)**2)
    END DO
    
END SUBROUTINE
!->+++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE  UPDATEINDUCEVF(THIS,OMG,R)
USE MFWM
IMPLICIT NONE
CLASS(FORWARD) :: THIS
REAL(RDT),INTENT(IN) :: OMG,R

REAL(RDT) :: PSI,RB
REAL(RDT) :: HRAD,HDT

INTEGER :: I,J,K,ERR

    HRAD=THIS%HDG*PI/180.0D0
    HDT=HRAD/OMG
    
    THIS%AEROF%BVQA=0.0D0
    
    CALL THIS%AEROF%UNIFLOW()!计算单个旋翼的平均诱导速度
    
    
!!    write(*,*) lami,lamd,mu,ctc,ALP
    IMSEL:SELECT CASE(THIS%AEROF%MKIM)
        CASE(0)!前飞均匀入流
        
            DO I=1,NPTS!按气动插值点循环
                DO J=1,THIS%NSTP!按一圈站位数循环
                    THIS%AEROF%LAMIL=THIS%AEROF%LAMI
                    THIS%AEROF%BVQA(3,I,J)=-THIS%AEROF%LAMIL*OMG*R
                END DO
            END DO
            
       CASE(1)!前飞线性入流
       
            DO I=1,NPTS
                DO J=1,THIS%NSTP
                    PSI=HRAD*(J-1)
                    RB=RLUDAUD(I)
                    CALL THIS%AEROF%LINFLOW(rb,PSI)
                    THIS%AEROF%BVQA(3,I,J)=-THIS%AEROF%LAMIL*OMG*R
                END DO
            END DO      
            
!!!       CASE(2)
!!!
!!!            HUBFMINFLOW(1,:)=HUBFCAER(3,:)
!!!            HUBFMINFLOW(2,:)=-HUBMTAER(1,:)
!!!            HUBFMINFLOW(3,:)=HUBMTAER(2,:)
!!!
!!!            fminflowstb(1)=HUBFCAERHW(1, 3)
!!!            fminflowstb(2)=-HUBMTAERHW(1, 1)
!!!            fminflowstb(3)=HUBMTAERHW(1, 2)
!!!
!!!
!!!            call dynamicinflow(3, nstp, HDT, iird, mdyim)
!!!       
!!!       
!!!            DO I=1,NBPL
!!!                DO J=1,NSTP
!!!                    PSI=HRAD*(J-1)
!!!                    RB=RLUDAUD(I)
!!!                    call dynainflow(lamil, rb, psi, indv(:, J))
!!!                    BVQA(3,I,J)=-LAMIL*OMG*R
!!!                END DO
!!!            END DO     
!!!        
!!!
       CASE(3)   

		CALL THIS%GETBPABC()

        call VTX%vortexSOL(THIS%AEROF%BVQA,THIS%AEROF%alp,THIS%AEROF%mu,&
                                    THIS%AEROF%LAMI,THIS%AEROF%nrbd,THIS%BLADE%R,THIS%AEROF%cchordrf,&
                                    THIS%omg,THIS%AEROF%nstp,THIS%AEROF%ct)
        
		CALL THIS%GETBPABC()

		
       CASE(4)
		
!	    CALL THIS%GETBPBM()
!
!		CALL VPM_RE()
!		
!		CALL THIS%GETBVQA()

    END select IMSEL
	
	CALL THIS%AEROF%TECPLOT_BVQA()
    
END SUBROUTINE
!->+++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE COPY(THIS,FOR)
IMPLICIT NONE
CLASS(FORWARD) :: THIS,FOR

    THIS%BLADE%y=FOR%BLADE%Y
    THIS%BLADE%yP=FOR%BLADE%YP
    THIS%BLADE%yST=FOR%BLADE%YST

    THIS%TSDT=FOR%TSDT
    THIS%tend=FOR%tend

    THIS%UQSA=FOR%UQSA
    THIS%AEROF%BVQA=FOR%AEROF%BVQA


END SUBROUTINE
!->+++++++++++++++++++++++++++++++++++++++++++++
subroutine fmlct(THIS,pe, qe, rxud,I) 
implicit none
CLASS(FORWARD) :: THIS
INTEGER,INTENT(IN) :: I
real(rdt), intent(IN):: RXUD
real(rdt), intent(out):: pe(6), qe(6)

real(rdt) ss
INTEGER LBI,LBIP1,J

REAL(RDT), EXTERNAL:: FWM 

    call rxlct2(LBI, ss, NPTS, rxud, RLUDAUD) 

    LBIP1=LBI+1
    DO J=1,6
        pe(J)=FWM(THIS%FPSI(J,LBI,I),THIS%FPSI(J,LBIP1,I),SS)
        Qe(J)=FWM(THIS%MPSI(J,LBI,I),THIS%MPSI(J,LBIP1,I),SS)
    END DO
    
end subroutine fmlct
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!SUBROUTINE GETBPBM(THIS)
!USE MFWM
!USE SASSEMBLE
!IMPLICIT NONE
!CLASS(FORWARD) :: THIS
!INTEGER :: NI,NJ,I,J,K
!real(rdt) :: SSI,SSJ
!REAL(RDT),ALLOCATABLE :: VPM_PSI(:),RXX(:),rlud(:)
!
!	IF(VTX%IS_ELASTIC.EQ.1 .OR. CITER.EQ.0) THEN
!		CALL THIS%GETBPABC()
!		
!		ALLOCATE(VPM_PSI(MAINROTOR%NAZIM))
!		DO I=1,MAINROTOR%NAZIM
!			VPM_PSI(I)=2*PI/MAINROTOR%NAZIM*(I-1)
!		END DO
!		ALLOCATE(RXX(MAINROTOR%NSPAN))
!		DO I=1,MAINROTOR%NSPAN
!			RXX(I)=(MAINROTOR%RX(I)+MAINROTOR%RX(I+1))/2
!		END DO
!		ALLOCATE(RLUD(NPTS))
!		RLUD(1)=RLUDAUD(1)*THIS%BLADE%R
!		DO I=2,NPTS
!			RLUD(I)=(RLUDAUD(I)+RLUDAUD(I-1))*THIS%BLADE%R/2
!		END DO
!
!		DO I=1,MAINROTOR%NAZIM
!			call rxlct2(NI, SSI, THIS%AEROF%NSTP, VPM_PSI(I), THIS%AEROF%PSIV) 
!			DO J=1,MAINROTOR%NSPAN+1
!				call rxlct2(NJ, SSJ, NPTS, MAINROTOR%RX(J),RLUDAUD*THIS%BLADE%R) 
!
!				DO K=1,3
!					CALL INTP2D(MAINROTOR%bpam(K,J,I),SSI,SSJ,&
!								VTX%BPAM(K,NJ,NI),VTX%BPAM(K,NJ,NI+1),&
!								VTX%BPAM(K,NJ+1,NI),VTX%BPAM(K,NJ+1,NI+1))
!					CALL INTP2D(MAINROTOR%BPCM(K,J,I),SSI,SSJ,&
!								VTX%BPCM(K,NJ,NI),VTX%BPCM(K,NJ,NI+1),&
!								VTX%BPCM(K,NJ+1,NI),VTX%BPCM(K,NJ+1,NI+1))
!				END DO
!			END DO
!		END DO
!		
!		DO I=1,MAINROTOR%NAZIM
!			call rxlct2(NI, SSI, THIS%AEROF%NSTP, VPM_PSI(I), THIS%AEROF%PSIV) 
!			DO J=1,MAINROTOR%NSPAN
!				call rxlct2(NJ, SSJ, NPTS, RXX(J),RLUD) 
!
!				DO K=1,3
!					CALL INTP2D(MAINROTOR%BPBM(K,J,I),SSI,SSJ,&
!									VTX%BPBM(K,NJ,NI),VTX%BPBM(K,NJ,NI+1),&
!									VTX%BPBM(K,NJ+1,NI),VTX%BPBM(K,NJ+1,NI+1))
!						CALL INTP2D(MAINROTOR%BVBM(K,J,I),SSI,SSJ,&
!									VTX%BVBM(K,NJ,NI),VTX%BVBM(K,NJ,NI+1),&
!									VTX%BVBM(K,NJ+1,NI),VTX%BVBM(K,NJ+1,NI+1))
!				END DO
!			END DO
!		END DO
!		
!		DEALLOCATE(VPM_PSI,RXX,RLUD)
!	END IF
!
!END SUBROUTINE
!
!SUBROUTINE GETBVQA(THIS)
!USE SASSEMBLE
!IMPLICIT NONE
!CLASS(FORWARD) :: THIS
!INTEGER :: NI,NJ,I,J,K
!real(rdt) :: SSI,SSJ
!REAL(RDT),ALLOCATABLE :: VPM_PSI(:)
!
!	ALLOCATE(VPM_PSI(MAINROTOR%NAZIM))
!	DO I=1,MAINROTOR%NAZIM
!		VPM_PSI(I)=2*PI/MAINROTOR%NAZIM*(I-1)
!	END DO
!
!	DO I=1,THIS%AEROF%NSTP
!		call rxlct2(NI, SSI, MAINROTOR%NAZIM, THIS%AEROF%PSIV(I), VPM_PSI) 
!		DO J=1,NPTS
!			call rxlct2(NJ, SSJ, MAINROTOR%NSPAN+1, RLUDAUD(J)*THIS%BLADE%R,MAINROTOR%RX) 
!			DO K=1,3
!				CALL INTP2D(THIS%AEROF%BVQA(K,J,I),SSI,SSJ,&
!							MAINROTOR%bVAm(K,NJ,NI),MAINROTOR%bVAm(K,NJ,NI+1),&
!							MAINROTOR%bVAm(K,NJ+1,NI),MAINROTOR%bVAm(K,NJ+1,NI+1))
!			END DO
!		END DO
!	END DO
!	
!	DEALLOCATE(VPM_PSI)
!END SUBROUTINE

SUBROUTINE GETBPABC(THIS)
USE MFWM
IMPLICIT NONE
CLASS(FORWARD) :: THIS
INTEGER :: I,J,K

	DO I=1,THIS%nstp
		THIS%ppsi=THIS%AEROF%PSIV(I)
		THIS%BLADE%UQS=THIS%UQSA(:, i)
		THIS%BLADE%UQS1T=THIS%UQS1TA(:, i)
		THIS%BLADE%UQS2T=THIS%UQS2TA(:, i)

		CALL THIS%BLADE%UPDATE()
		call VTX%calbpabc(THIS%BLADE,THIS%OMG,THIS%PPSI,&
		THIS%TH0,THIS%TH1C,THIS%TH1S,THIS%AEROF%cchordrf,I)
	END DO
	
END SUBROUTINE    

END MODULE
