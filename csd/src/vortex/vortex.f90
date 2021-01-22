module VORTEX_CLASS
USE GlobalDataFun
USE MBEAMAERO
USE GlobalDataFun
USE MMATHS
IMPLICIT NONE
!->+++++++++++++++++++++++++++++++++++++++++++++
TYPE,PUBLIC :: VORTEX

   !----------------------------------------------------
    !�涨���ؽ�Ҷ1/4���ߵĳ�ʼλ��ΪX�ᣬ��ֱ��������ΪZ�ᣬ��������ΪԲ�ĵĵ�������ϵ
    !�涨��Ҷ��x��ת��y�ᣬ�ұ��Խ��Ľ�Ҷ������λ��Խ��
    !��ѧ����-------------------------------------------------
!    real(rdt),parameter::pi=3.1415926
	integer :: inp_cir=0
    !��Ҷ����-------------------------------------------------
    integer::nf,&                   !����β��Ȧ��
                  nb                        !��ҶƬ��
    real(rdt)::r,&	                        !��Ҷ�뾶
    		   ri0,&	                    !��������������� 
    		   omg			            !ת��
    !��Ҷ����-------------------------------------------------
    real(rdt)::ct                        !����ϵ��
       

    !��ʱ����------------------------------------------------
    real(rdt)::t_start,t_end                !��ʼ����ʱ��
    
    !�˶���---------------------------------------------------
    real(rdt)::u,&		                    !ǰ����
               vdx,&                        !Ԥ���յ��ٶ�
    	       alphs    		            !����ǰ���(̧ͷΪ��)
    
    !�������-------------------------------------------------
    integer ::   na,nat,nx,&	            !β��һ�ܵĽڵ���
    			       nw,&		        !β����Ȧ��
    			       ns,&		        !��Ҷ�ֶ���
    		           nseg		        !��Ҷβ���зֶ�����
    real(rdt),allocatable ::vor_tip(:)		    !ĳ��λ�������ڸ�����Ҷ����չ����ͷ�λ��
    real(rdt):: vor_root,rcct               !�������ͷ�λ��

    real(rdt) :: CHB
	INTEGER :: IS_ELASTIC
	REAL(RDT) :: RELAX
    !�������-------------------------------------------------
    integer::l_iter=0	                !��ǰ��������
    real(rdt)::rms_wake=1                 !β���仯��
    real(rdt),allocatable :: tmax(:),&		            !�����ֵ
               taoroot(:),&                !�����л���
	           tao(:),&	                !���н�Ҷĳһ��λ���¸��λ������ڶ�ά����Ϊ1ʱΪ��ǰ��λ��Ҷ����
	           bpq(:,:,:),&		        !��Ҷ������λ����1/4���ߵ�����(��λά��Ϊ1ʱΪ0��λ��)
	           bpb(:,:,:),&		        !��Ҷ������λ���º�Ե������(��λά��Ϊ1ʱΪ0��λ��)
	           bpc(:,:,:),&		        !��Ҷ������λ���¿��Ƶ�����(��λά��Ϊ1ʱΪ0��λ��)
	           grid(:,:,:,:),&     !��β��һȦ�ڵ�����
	           tip(:,:,:),&          !����������
	           vtip(:,:,:),&         !�������ٶ�
	           root(:,:,:),&         !����������
	           vgrid(:,:,:,:),&      !��β���ٶ�
	           vbpc(:,:,:),&           !���Ƶ��˶��ٶ�
	           vbpq(:,:,:)

        REAL(RDT), ALLOCATABLE:: bpam(:,:, :),bpbm(:,:, :),bpcm(:,:, :)
        REAL(RDT), ALLOCATABLE:: bvam(:,:, :),bvbm(:,:, :),bvcm(:,:, :)

	!-------------------------------------------------------------------
    CONTAINS
    PROCEDURE,PUBLIC :: sub_bwake
    PROCEDURE,PUBLIC :: stl
    PROCEDURE,PUBLIC :: GET_RC
    PROCEDURE,PUBLIC :: sub_rc
    PROCEDURE,PUBLIC :: solv_circul
    PROCEDURE,PUBLIC :: SET_circul
    PROCEDURE,PUBLIC :: sub_ve
    PROCEDURE,PUBLIC :: sub_tip
    PROCEDURE,PUBLIC :: sub_vi
    PROCEDURE,PUBLIC :: sub_newake
    PROCEDURE,PUBLIC :: sub_vwake
    PROCEDURE,PUBLIC :: sub_out
    PROCEDURE,PUBLIC :: CONSTRUCT_VOTEX
    PROCEDURE,PUBLIC :: vortexSOL
    PROCEDURE,PUBLIC :: calbpabc
    PROCEDURE,PUBLIC :: vor_allo
	
    PROCEDURE,PUBLIC :: GET_VEL
    PROCEDURE,PUBLIC :: WPTEST
	
END TYPE VORTEX

CONTAINS

subroutine CONSTRUCT_VOTEX(THIS,INPU,NSTP,NPTS)
use GlobalDataFun
USE INPUT_CLASS
implicit none
CLASS(VORTEX) :: THIS
TYPE(INPUT) :: INPU
INTEGER,INTENT(IN) :: NSTP,NPTS

INTEGER LIN,ISR
real(rdt) ::  vor_rootr,vor_tipr
character(len=100) :: copyfilenam
integer :: err,i

    THIS%ri0=INPU%vtx%ri0
    THIS%nf=INPU%vtx%nf
    THIS%nw=INPU%vtx%nw
    THIS%nseg=INPU%vtx%nseg
    THIS%nx=INPU%vtx%nx
    vor_rootr=INPU%vtx%vor_rootr
    vor_tipr=INPU%vtx%vor_tipr
    THIS%rcct=INPU%vtx%rcct
    THIS%ns=NPTS-1
    this%na=(NSTP-1)/THIS%nx
    this%IS_ELASTIC=INPU%vtx%IS_ELASTIC
    this%RELAX=INPU%vtx%RELAX
	
    ALLOCATE(THIS%vor_tip(this%na),THIS%tmax(this%na),THIS%taoroot(this%na),THIS%tao(THIS%ns*this%na), &
            THIS%bpq(3,THIS%ns+1,this%na),THIS%bpb(3,THIS%ns+1,this%na),&
            THIS%bpc(3,THIS%ns,this%na),this%grid(3,THIS%nseg+1,THIS%ns+1,this%na), &
            THIS%tip(3,this%na*THIS%nw+1,this%na),THIS%vtip(3,this%na*THIS%nw+1,this%na),&
            THIS%root(3,this%na*THIS%nw+1,this%na),THIS%vgrid(3,THIS%nseg+1,THIS%ns+1,this%na), &
            THIS%vbpc(3,THIS%ns,this%na),THIS%vbpq(3,THIS%ns+1,this%na),&
            STAT=ERR)
    
    IF(ERR.NE.0) THEN
        WRITE(*,998)
        STOP
    END IF  
    
    CALL THIS%vor_allo(NSTP)
    
    THIS%tmax=0.0
    THIS%vor_root=vor_rootr
    
    do i=1,this%na
        THIS%vor_tip(i)=vor_tipr
    end do
	       
998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
end subroutine


!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
subroutine vor_allo(THIS,nstp)
implicit none
CLASS(VORTEX) :: THIS
INTEGER,INTENT(IN) :: NSTP
integer err


      ALLOCATE(THIS%bpam(3,NPTS,nstp),THIS%bpbm(3,NPTS,nstp),&
                THIS%bpcm(3,NPTS,nstp),THIS%bvam(3,NPTS,nstp),&
                THIS%bvbm(3,NPTS,nstp),THIS%bvcm(3,NPTS,nstp),&
                  STAT=ERR)
                  
	IF(ERR.NE.0) THEN
		WRITE(*, 998)
		STOP
	END IF
	
998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
end subroutine vor_allo
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!initial��ʼ����--------------------------------------------------------------
!��ʱû�õ�
subroutine vortexSOL(THIS,bvqa,alp,mu,LAMI,nrbd,rb,CHB,omgb,nstp,ctx)
implicit none
CLASS(VORTEX) :: THIS
INTEGER,INTENT(IN)::nrbd,nstp
real(rdt),INTENT(IN):: alp,mu,rb,CHB,omgb,ctx,lami
real(rdt),INTENT(OUT):: bvqa(3,NPTS,nstp,2)
real(rdt)::tmp1,tmp2,tmp3
integer :: i,j,k,natt
REAL(KIND=RDT) :: PSI0=0

	THIS%bvbm(1,:,:)=-THIS%bvbm(1,:,:)+mu*omgb*rb
	THIS%bvbm(2,:,:)=-THIS%bvbm(2,:,:)
	THIS%bvbm(3,:,:)=-THIS%bvbm(3,:,:)-mu*omgb*rb*tan(alp)
	!THIS%bvbm(1,:,:)=mu*omgb*rb
	!THIS%bvbm(2,:,:)=0
	!THIS%bvbm(3,:,:)=-mu*omgb*rb*tan(alp)
    this%alphs=-alp
    THIS%u=mu
    THIS%nb=nrbd
    THIS%r=rb
    this%CHB=CHB
    THIS%omg=omgb
    natt=nstp-1

    THIS%ct=ctx
    THIS%vdx=-lami

!            if(THIS%l_iter .eq.0) then
!                call this%CONSTRUCT_VOTEX(inpu,natt)
!            end if


    if(THIS%l_iter .eq.0) then
		THIS%bpq=THIS%bpam(:,1:NPTS,1:natt:THIS%nx)
		THIS%bpb=THIS%bpcm(:,1:NPTS,1:natt:THIS%nx)
		THIS%bpc=THIS%bpbm(:,2:NPTS,1:natt:THIS%nx)
		THIS%vbpc=THIS%bvbm(:,2:NPTS,1:natt:THIS%nx)
        call THIS%sub_bwake()
        call this%sub_out(0)
    end if
          
    THIS%rms_wake=1.0D0

    do while(THIS%rms_wake>5.0e-4)
	    IF(this%IS_ELASTIC.EQ.1) THEN
			THIS%bpq=THIS%bpam(:,1:NPTS,1:natt:THIS%nx)
			THIS%bpb=THIS%bpcm(:,1:NPTS,1:natt:THIS%nx)
			THIS%bpc=THIS%bpbm(:,2:NPTS,1:natt:THIS%nx)
			THIS%vbpc=THIS%bvbm(:,2:NPTS,1:natt:THIS%nx)
		END IF
	
        THIS%l_iter=THIS%l_iter+1
		
		if(this%inp_cir.eq.0) then
			call this%solv_circul()
		end if
		
        call this%sub_tip()
        call this%sub_newake()

		call this%sub_out(2)
		call this%sub_out(22)
		call this%sub_out(3)
		call this%sub_out(5)
		call this%sub_out(8)
		call this%sub_out(1)
		call this%sub_out(4)
		!call this%WPTEST(PSI0)

        if(THIS%l_iter.gt.100) then
            if(mod(THIS%l_iter,10).eq.1) exit
        end if
    end do
    
    bvqa=0.0D0
    do i=1,this%na
        do j=1,THIS%ns+1
            call this%sub_vi(i,THIS%bpq(:,j,i),BVQA(:,j,1+(i-1)*THIS%nx,1))!����
        end do
    end do
     BVQA(:,:,nstp,1)=BVQA(:,:,1,1)!����
    do i=1,this%na
        do j=1,THIS%nx
            bvqa(:,:,j+(i-1)*THIS%nx,1)=bvqa(:,:,1+(i-1)*THIS%nx,1)+&!����
                                              (bvqa(:,:,1+i*THIS%nx,1)-bvqa(:,:,1+(i-1)*THIS%nx,1))/THIS%nx*(j-1.0)!����
        end do
    end do

end subroutine

subroutine vortexdeallo(THIS)
implicit none
CLASS(VORTEX) :: THIS
integer :: err

	
    DEALLOCATE(THIS%vor_tip,THIS%tmax,THIS%taoroot,THIS%tao, &
            THIS%bpq,THIS%bpb,THIS%bpc,this%grid, &
            THIS%tip,THIS%vtip,THIS%root,THIS%vgrid, &
            THIS%vbpc,THIS%vbpq,&
            STAT=ERR)

    IF(ERR.NE.0) THEN
        WRITE(*,999)
        STOP
    END IF  
    
998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
end subroutine

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
subroutine calbpabc(THIS,BLADE,OMG,PSI,TH0, TH1C, TH1S,cchordrf,ICT)
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
USE EIGENBEAM_CLASS
IMPLICIT NONE
CLASS(VORTEX) :: THIS
CLASS(BEAM) :: BLADE
REAL(RDT),INTENT(IN) :: OMG,PSI,cchordrf,TH0, TH1C, TH1S
INTEGER,INTENT(IN) :: ICT
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface:
!real(rdt),intent(out) :: bpa(3 ,nstp , nbpl),bpb(3 ,nstp , nbpl),bpc(3 ,nstp , nbpl)
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->        
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
INTEGER NBPLI
real(rdt) cdlct1(3,2),RPN1(3,2),RVN1(3,2)
real(rdt) cdlct2(3,1),RPN2(3,1),RVN2(3,1)
real(rdt) rxud
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> main body:
REAL(RDT),EXTERNAL :: triseries,fabsc

    CALL BLADE%GETTHP(TH0, TH1C, TH1S, PSI, OMG)

    DO NBPLI=1,NPTS
        cdlct1(:,1)=(/0.0D0,0.0D0,0.0D0/)
        cdlct1(:,2)=(/0.0D0,-3.0D0/4*cchordrf,0.0D0/)
        rxud=RLUDAUD(NBPLI)*BLADE%R

        call BLADE%rvndlct(RPN1, RVN1,2,RXUD,cdlct1,OMG,PSI,0,1)

        THIS%BPAM(:,NBPLI,ICT)=RPN1(:,1)
        THIS%BVAM(:,NBPLI,ICT)=RVN1(:,1)
        THIS%BPCM(:,NBPLI,ICT)=RPN1(:,2)
        THIS%BVCM(:,NBPLI,ICT)=RVN1(:,2)

        IF(NBPLI.GT.1) THEN
            rxud=(RLUDAUD(NBPLI)+RLUDAUD(NBPLI-1))/2.0D0*BLADE%R
        ELSEIF(NBPLI.EQ.1) THEN
            rxud=RLUDAUD(NBPLI)*BLADE%R
        ELSE 
            WRITE(*,*) "ERROR,IN SUBROUTINE BPABC"
            stop
        END IF
        
        cdlct2(:,1)=(/0.0D0,-1.0D0/2*cchordrf,0.0D0/)
        call BLADE%rvndlct(RPN2, RVN2,1,RXUD,cdlct2,OMG,PSI,0,1)
        THIS%BPBM(:,NBPLI,ICT)=RPN2(:,1)
        THIS%BVBM(:,NBPLI,ICT)=RVN2(:,1)
    END DO
         
end subroutine calbpabc
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!bwake��ʼβ��--------------------------------------------------------------------
subroutine sub_bwake(THIS)
!���this%grid(1:3,1:nseg+1,1:ns+1,1:this%na)
!     tip(1:3,1:this%na*nf+1,1:this%na)
!     root(1:3,1:this%na*nf+1,1:this%na)
implicit none
CLASS(VORTEX) :: THIS
real(rdt) tep,tmp,tmp1(3),tmp2(3),tmpr,fi,dpsi,psi0,psi1,psi,eta,v
integer i,j,n

    v=this%u*tan(this%alphs)+this%vdx
    dpsi=2.0*PI/real(this%na)	        
    do j=1,this%na
    	psi0=(j-1.0)*dpsi
!///////////////////////////��β��//////////////////////////////////////////
        do i=1,this%ns+1
            tmpr=SQRT(this%bpb(1,i,j)**2+this%bpb(2,i,j)**2)
            psi1=acos(SQRT(this%bpq(1,i,j)**2+this%bpq(2,i,j)**2)/tmpr) 
            psi=psi0-psi1
	        do n=1,this%nseg+1
	        	eta=(n-1.)*dpsi 
		        this%grid(1,n,i,j)=tmpr*COS(psi-eta)+this%u*eta*this%r
		        this%grid(2,n,i,j)=tmpr*SIN(psi-eta)
		        this%grid(3,n,i,j)=this%bpb(3,i,j)+v*eta*this%r
	        end do
        end do
!///////////////////////////������//////////////////////////////////////////                        
        tmp=(this%vor_tip(j)-this%ri0)/(1-this%ri0)
        tmp1(:)=(1-tmp)*this%bpb(:,1,j)+tmp*this%bpb(:,this%ns+1,j) 
        tmp2(:)=(1-tmp)*this%bpq(:,1,j)+tmp*this%bpq(:,this%ns+1,j)
        tmpr=SQRT(tmp1(1)**2+tmp1(2)**2)
        psi1=acos(SQRT(tmp2(1)**2+tmp2(2)**2)/tmpr)
        psi=psi0-psi1
        do i=1,this%na*this%nw+1
            eta=(i-1.)*dpsi
            this%tip(1,i,j)=tmpr*COS(psi-eta)+this%u*eta*this%r
            this%tip(2,i,j)=tmpr*SIN(psi-eta)
            this%tip(3,i,j)=tmp1(3)+v*eta*this%r
        end do
!///////////////////////������/////////////////////////////////////////
        tmp=(this%vor_root-this%ri0)/(1-this%ri0)
        tmp1(:)=(1-tmp)*this%bpb(:,1,j)+tmp*this%bpb(:,this%ns+1,j) 
        tmp2(:)=(1-tmp)*this%bpq(:,1,j)+tmp*this%bpq(:,this%ns+1,j)
        tmpr=SQRT(tmp1(1)**2+tmp1(2)**2)
        psi1=acos(SQRT(tmp2(1)**2+tmp2(2)**2)/tmpr)
        psi=psi0-psi1
        do i=1,this%na*this%nw+1
            eta=(i+this%nseg-1.)*dpsi
            this%root(1,i,j)=tmpr*COS(psi-eta)+this%u*eta*this%r
            this%root(2,i,j)=tmpr*SIN(psi-eta)
            this%root(3,i,j)=tmp1(3)+v*eta*this%r
        end do
!/////////////////////////////////////////////////////////////////////////
    end do
end subroutine

!ve����β���Խڵ���յ��ٶ�---------------------------------------------------------------
subroutine sub_ve(this,k,p,vb,vg,vt,vr)
!���vb(ns*this%na,3)
!    vg(ns*this%na,3)
!    vt(3,nb)
!    vr(3,nb)
implicit none
class(vortex) :: this
integer,intent(in) :: k
real(rdt),intent(in ) :: p(3)
real(rdt),intent(out) :: vb(this%ns*this%na,3),vg(this%ns*this%na,3),vt(3,this%nb),vr(3,this%nb)

integer i,j,ki,kg,bi,loop,tmp
real(rdt) a(3),b(3),c(3),d(3),rc,v1(3),v2(3),v3(3),v4(3),RC0,AVG
    

    vb=0.0
    vg=0.0
    vt=0.0
    vr=0.0
    rc0=0.075
    avg=sum(this%tmax)/this%na
    do bi=1,this%nb
        tmp=k+(bi-1)*this%na/this%nb
        ki=tmp-this%na*int((tmp-1)/this%na)
        do i=1,this%ns
            a(:)=this%bpq(:,i,ki)
            b(:)=this%bpq(:,i+1,ki)
            c(:)=this%bpb(:,i+1,ki)
            d(:)=this%bpb(:,i,ki)
            rc=rc0*this%chb
            call this%stl(p,d,a,rc,v1)
            call this%stl(p,a,b,rc,v2)
            call this%stl(p,b,c,rc,v3)
            loop=(ki-1)*this%ns+i
            vb(loop,:)=v1(:)+v2(:)+v3(:) 
            
            a(:)=this%grid(:,1,i,ki)
            b(:)=this%grid(:,1,i+1,ki)
            c(:)=this%grid(:,2,i+1,ki)
            d(:)=this%grid(:,2,i,ki)
            rc=rc0*this%chb
            call this%stl(p,b,c,rc,v1)
            call this%stl(p,d,a,rc,v3)
            rc=rc0*this%chb
            call this%stl(p,c,d,rc,v2)
            loop=(ki-1)*this%ns+i
            vg(loop,:)=v1(:)+v2(:)+v3(:)
        end do
        do j=2,this%nseg
            do i=1,this%ns
                a(:)=this%grid(:,j,i,ki)
                b(:)=this%grid(:,j,i+1,ki)
                c(:)=this%grid(:,j+1,i+1,ki)
                d(:)=this%grid(:,j+1,i,ki)
!/////////////////////////������/////////////////////////////
                RC=RC0*J*THIS%CHB
                call this%stl(p,a,b,rc,v1)
                call this%stl(p,c,d,rc,v3)
!////////////////////////β����//////////////////////////////
                RC=RC0*J*THIS%CHB
                call this%stl(p,b,c,rc,v2)
                call this%stl(p,d,a,rc,v4)
                kg=ki-j+1
                if(kg<=0) kg=kg+this%na
                loop=(kg-1)*this%ns+i
                vg(loop,:)=v1(:)+v2(:)+v3(:)+v4(:)
!////////////////////////////////////////////////////////////
            end do
        end do
!/////////////////////////���⽰����/////////////////////////////
        do i=1,this%na*this%nw
            a(:)=this%tip(:,i,ki)
            b(:)=this%tip(:,i+1,ki)
            c(:)=this%root(:,i+1,ki)
            d(:)=this%root(:,i,ki)
            IF(THIS%rcct.gt.0.0) THEN
                CALL THIS%sub_rc(i,rc,this%na)
            ELSE IF(THIS%rcct.LE.0.0) THEN
                call this%get_rc(i,AVG,rc)
            END IF
            call this%stl(p,a,b,rc,v1)
            call this%stl(p,c,d,rc,v2)
            vt(:,bi)=vt(:,bi)+v1(:)
            vr(:,bi)=vr(:,bi)+v2(:)         
        end do                
    end do
end subroutine

SUBROUTINE GET_RC(this,N,AVG,RC)
USE GlobalDataFun
implicit none
class(vortex) :: this
INTEGER,INTENT(IN) :: N
REAL(KIND=rdt) ETA,RC,TMP,AVG,MU,DET,A1

      MU =1.46D-4; A1=ABS(THIS%rcct)
      ETA=(N*2.0/THIS%NA+131.6/180.0)*PI
      DET=1+A1*AVG/MU
      TMP=4*1.25643*MU*DET
      RC =SQRT(TMP*ETA/this%OMG)

END SUBROUTINE GET_RC
!sub_rc����к˰뾶-----------------------------------------------------------------------------
subroutine sub_rc(this,n,rc,na)
USE GlobalDataFun
implicit none
class(vortex) :: this
integer,intent(in) :: n,na
real(rdt),intent(out) :: rc
real(rdt) eta,tmp
    eta=(n*2.0/na+131.6/180.0)*PI
    tmp=4*1.256*1.46e-4*this%rcct
    rc=SQRT(tmp*eta/this%omg)
end subroutine

!stl�ϳ������ֱ���ж��յ��ٶ�----------------------------------------------
subroutine stl(this,p,a,b,rc,v)
implicit none
class(vortex) :: this
real(rdt),intent(in) :: p(3),a(3),b(3),rc
real(rdt),intent(out) :: v(3)

real(rdt)  s(3),r(3),sm(3),tmp(4),rm

	s(1)=(p(2)-a(2))*(p(3)-b(3))-(p(3)-a(3))*(p(2)-b(2))
	s(2)=(p(3)-a(3))*(p(1)-b(1))-(p(1)-a(1))*(p(3)-b(3))
	s(3)=(p(1)-a(1))*(p(2)-b(2))-(p(2)-a(2))*(p(1)-b(1))
	r(1)=SQRT((p(1)-a(1))**2+(p(2)-a(2))**2+(p(3)-a(3))**2)	
	r(2)=SQRT((p(1)-b(1))**2+(p(2)-b(2))**2+(p(3)-b(3))**2)	
 	r(3)=(a(1)-b(1))**2+(a(2)-b(2))**2+(a(3)-b(3))**2
	tmp(1)=(p(1)-a(1))*(p(1)-b(1))+(p(2)-a(2))&
           *(p(2)-b(2))+(p(3)-a(3))*(p(3)-b(3))
	tmp(2)=r(1)*r(1)*r(2)*r(2)-tmp(1)*tmp(1)+&
           rc*rc*(r(1)*r(1)+r(2)*r(2)-2.*tmp(1))
	sm(1)=(p(1)-a(1))*(r(2)*r(2)-tmp(1))+&
          (p(1)-b(1))*(r(1)*r(1)-tmp(1))
	sm(2)=(p(2)-a(2))*(r(2)*r(2)-tmp(1))+&
          (p(2)-b(2))*(r(1)*r(1)-tmp(1))
	sm(3)=(p(3)-a(3))*(r(2)*r(2)-tmp(1))+&
          (p(3)-b(3))*(r(1)*r(1)-tmp(1))       
	rm=SQRT(sm(1)*sm(1)+sm(2)*sm(2)+sm(3)*sm(3))/r(3)
    if(rm<0.001*rc) then
 		v=0.0
	else 
 		tmp(3)=(r(1)+r(2))*(1.-tmp(1)/r(1)/r(2))
 		tmp(4)=1./4./3.14159/tmp(2)*tmp(3)
      	v(:)=s(:)*tmp(4)
  	end if

end subroutine	

!solv_circul��⻷��---------------------------------------------------------------------------
subroutine solv_circul(this)
implicit none
class(vortex) :: this
integer i,k,ki,m,ip1,ip2,ip3,flap,tmp 
real(rdt) psi0,psi1,beta,dbet,cp(3),tmpr
real(rdt) ht,hr
real(rdt),allocatable :: tmp1(:),tmp2(:),tmp3(:),tmp4(:),a1(:),a2(:),b1(:),b2(:),c(:),d(:)
real(rdt),allocatable :: vb(:,:),vg(:,:),vt(:,:),vr(:,:)
real(rdt),allocatable :: ma(:,:),mb(:),hg(:)
integer :: err


  ALLOCATE(tmp1(3),tmp2(3),tmp3(3),tmp4(3),a1(3),a2(3),b1(3),b2(3),c(3),d(3), &
   vb(this%ns*this%na,3),vg(this%ns*this%na,3),vt(3,this%nb),vr(3,this%nb), &
   ma(this%ns*this%na,this%ns*this%na),mb(this%ns*this%na),hg(this%na*this%ns), &
        STAT=ERR)

IF(ERR.NE.0) THEN
    WRITE(*,998)
    STOP
END IF  

ma=0.0
mb=0.0
do k=1,this%na
    psi0=(k-1.0)*2.0*PI/real(this%na) 
    do i=1,this%ns
        cp(:)=this%bpc(:,i,k)
!////////////////////////������Ƶ㷽������////////////////////////////////////
		tmp1(:)=this%bpq(:,i,k)
        tmp2(:)=this%bpq(:,i+1,k)
        tmp3(:)=this%bpb(:,i,k)
        tmp4(:)=this%bpb(:,i+1,k)
        call xl(cp,tmp2,a1)
        call xl(cp,tmp1,a2)
        call xl(cp,tmp3,b1)
        call xl(cp,tmp4,b2)
        call wi(a1,a2,tmp1)
        call wi(a2,b1,tmp2)
        call wi(b1,b2,tmp3)
        call wi(b2,a1,tmp4)
        call fang(tmp1,a1)
        call fang(tmp2,a2)
        call fang(tmp3,b1)
        call fang(tmp4,b2)
		d=(a1+a2+b1+b2)/4.0
!////////////////////////������յ��ٶ�////////////////////////////////////////
        tmp3(:)=this%vbpc(:,i,k)
	  !TMPR=SQRT(CP(1)**2+CP(2)**2)
	  !PSI1=ACOS(SQRT(TMP1(1)**2+TMP1(2)**2)/TMPR)
      !TMP3(1)= TMP3(1)+this%OMG*TMPR*SIN(PSI0-PSI1)
      !TMP3(2)= TMP3(2)-this%OMG*TMPR*COS(PSI0-PSI1)
      !TMP3(3)= TMP3(3)
        ip1=(k-1)*this%ns+i 
        mb(ip1)=-dot_product(tmp3,d) 
!////////////////////////�����յ��ٶ�//////////////////////////////////////////
        call this%sub_ve(k,cp,vb,vg,vt,vr)
        ma(ip1,:)=matmul(vb,d)+matmul(vg,d)
!//////////////////////////////////////////////////////////////////////////////
        do m=1,this%nb
            ht=dot_product(vt(:,m),d) 
            hr=dot_product(vr(:,m),d) 
            tmp=k+(m-1)*this%na/this%nb
            ki=tmp-this%na*int((tmp-1)/this%na)
            if (this%l_iter==1)then
                ip2=(ki-1)*this%ns+this%ns-1
                ip3=(ki-1)*this%ns+1
                ma(ip1,ip2)=ma(ip1,ip2)+ht
                ma(ip1,ip3)=ma(ip1,ip3)+hr
            else
                mb(ip1)=mb(ip1)-ht*this%tmax(ki)
                mb(ip1)=mb(ip1)-hr*this%taoroot(ki)
            end if
        end do
        !//////////////////////////////////////////////////////////////////////////////
    end do
end do

    call linsolver2(this%tao, 1, this%ns*this%na, ma, mb)

  DEALLOCATE(tmp1,tmp2,tmp3,tmp4,a1,a2,b1,b2,c,d, &
   vb,vg,vt,vr,&
   ma,mb,hg, &
        STAT=ERR)

IF(ERR.NE.0) THEN
    WRITE(*,999)
    STOP
END IF  

998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
end subroutine
!--------------------------------------------------------------
subroutine SET_circul(this,np,ns,psit,spant,tao)
implicit none
class(vortex) :: this
integer :: np,ns
real(rdt),intent(in) :: psit(np),spant(ns),tao(ns,np)

real(rdt) :: taop(ns,np)
integer :: j,k,ni,nj,idx
real(rdt) :: dpsi,ssi,ssj,PP,RR

	taop=tao*this%omg*this%r*this%CHB

    dpsi=2.0*PI/real(this%na)	        
    do j=1,this%na
		pp=(j-1.0)*dpsi
		do k=1,this%ns
			rr=SQRT(this%bpc(1,k,1)**2+this%bpc(2,k,1)**2)/this%r

			call rxlct2(NI, SSI, np, pp, psit) 
			call rxlct2(NJ, SSJ, ns, rr,spant) 
   
			idx=k+(j-1)*this%ns
			CALL INTP2D(this%tao(idx),SSI,SSJ,&
						taop(NJ,NI),taop(NJ,NI+1),&
						taop(NJ+1,NI),taop(NJ+1,NI+1))
		end do
	end do

end subroutine
!tip�󽰼��кͽ�������ǿ--------------------------------------------------------------
subroutine sub_tip(this)
implicit none
class(vortex) :: this
integer i,k,j,ipoint
real(rdt) tmp1,tmp2

    ipoint=this%ns
    do k=1,this%na
        j=1+(k-1)*this%ns
		this%tmax(k)=this%tao(j)
        this%taoroot(k)=this%tao(j+3)
	    do i=2,this%ns
            j=i+(k-1)*this%ns
		    if(this%tao(j)>this%tmax(k)) then
			    this%tmax(k)=this%tao(j)
			    ipoint=i
		    end if
	    end do
    end do
    
    !do k=1,this%na
    !    tmp2=0.0
    !    do i=ipoint+1,this%ns
    !        j=i+(k-1)*this%ns
    !        tmp1=SQRT((this%bpq(1,i+1,1)-this%bpq(1,i,1))**2+(this%bpq(2,i+1,1)-this%bpq(2,i,1))**2)
    !        tmp2=tmp2+this%tao(j)*tmp1
    !    end do
    !    tmp1=SQRT(this%bpq(1,ipoint,1)**2+this%bpq(2,ipoint,1)**2)               
    !    this%vor_tip(k)=(tmp1+tmp2/this%tmax(k))/this%r        
    !end do

end subroutine

!vi��ռ�����һ����ٶ�--------------------------------------------------------------
subroutine sub_vi(this,k,p,v)
implicit none
class(vortex) :: this
integer,intent(in) :: k
real(rdt),intent(in) ::p(3)
real(rdt),intent(out) ::v(3)
integer i,j,ki,tmp

real(rdt),allocatable :: vb(:,:),vg(:,:),vt(:,:),vr(:,:)
real(rdt),allocatable :: a(:),b(:),vbi(:),vgi(:),vfi(:)
integer :: err
    
    
ALLOCATE(vb(this%ns*this%na,3),vg(this%ns*this%na,3),vt(3,this%nb),vr(3,this%nb), &
    a(this%ns*this%na),b(this%ns*this%na),vbi(3),vgi(3),vfi(3), &
    STAT=ERR)

    IF(ERR.NE.0) THEN
        WRITE(*,998)
        STOP
    END IF  
    
    call this%sub_ve(k,p,vb,vg,vt,vr)
    vfi=0.0
    do i=1,3
        a(:)=vb(:,i)
        vbi(i)=dot_product(a,this%tao) 
		a(:)=vg(:,i)
	    vgi(i)=dot_product(a,this%tao) 
    end do
    do i=1,this%nb
        tmp=k+(i-1)*this%na/this%nb
        ki=tmp-this%na*int((tmp-1)/this%na)
        vfi(:)=vfi(:)+vt(:,i)*this%tmax(ki)
        vfi(:)=vfi(:)+vr(:,i)*this%taoroot(ki)
    end do
    v=vbi+vgi+vfi
    
DEALLOCATE(vb,vg,vt,vr,a,b,vbi,vgi,vfi, &
    STAT=ERR)
IF(ERR.NE.0) THEN
    WRITE(*,999)
    STOP
END IF   
    
998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
end subroutine


!vwake���㵱ǰ��λ��������β���ٶ�-------------------------------------------------------------
subroutine sub_vwake(this)
!���vgrid(1:3,1:this%nseg+1,1:this%ns+1,1:this%na)
!     vtip(1:3,1:this%na*nf+1,1:this%na)
implicit none
class(vortex) :: this
integer i,j,k
real(rdt) p(3),v(3)

    this%vtip=0.0D0
    do k=1,this%na
        do j=1,this%ns+1
            do i=1,this%nseg+1
                p(:)=this%grid(:,i,j,k)
                call this%sub_vi(k,p,v)
                this%vgrid(1,i,j,k)=v(1)+this%u*this%omg*this%r
                this%vgrid(2,i,j,k)=v(2)
                this%vgrid(3,i,j,k)=v(3)+this%u*this%omg*this%r*tan(this%alphs)
            end do
        end do
	    do i=1,this%na*this%nf+1
	        p(:)=this%tip(:,i,k)
	        call this%sub_vi(k,p,v)
            this%vtip(1,i,k)=v(1)+this%u*this%omg*this%r
            this%vtip(2,i,k)=v(2)
            this%vtip(3,i,k)=v(3)+this%u*this%omg*this%r*tan(this%alphs)
        end do
    end do

end subroutine	

!newwake����ȫ������β��-----------------------------------
subroutine sub_newake(this)
implicit none
class(vortex) :: this
integer i,j,k,k_1
real(rdt),allocatable::fgrid(:,:,:,:),fvgrid(:,:,:,:)
real(rdt),allocatable::ftip(:,:,:),fvtip(:,:,:)
real(rdt) t,tmp,tmpm(3)
integer :: err
    
    
ALLOCATE(fgrid(3,this%nseg+1,this%ns+1,this%na),fvgrid(3,this%nseg+1,this%ns+1,this%na), &
            ftip(3,this%na*this%nw+1,this%na),fvtip(3,this%na*this%nw+1,this%na), &
            STAT=ERR)

IF(ERR.NE.0) THEN
    WRITE(*,998)
    STOP
END IF  
    
    !if(this%u>0.04) relax=0.5
    !if(this%u>0.09) relax=0.7
    t=2*PI/real(this%na)/this%omg
    call this%sub_vwake()
!///////////////////////////////////////////////////////////////////////////////////
    fgrid=this%grid
    fvgrid=this%vgrid
    ftip=this%tip
    fvtip=this%vtip
!//////////////////////////////////Ԥ�ⲽ///////////////////////////////////////////
    do k=1,this%na
        k_1=k-1
        if(k==1) k_1=this%na
        do j=1,this%ns+1
            this%grid(:,1,j,k)=this%bpb(:,j,k)
            do i=2,this%nseg+1
                this%grid(:,i,j,k)=this%grid(:,i-1,j,k_1)+0.25*(this%vgrid(:,i-1,j,k)+this%vgrid(:,i,j,k)+&
                              this%vgrid(:,i-1,j,k_1)+this%vgrid(:,i,j,k_1))*t
            end do
        end do
        this%tip(:,1,k)=(this%bpb(:,this%ns,k)+this%bpb(:,this%ns+1,k))/2.0
        do i=2,this%na*this%nf+1
            this%tip(:,i,k)=this%tip(:,i-1,k_1)+0.25*(this%vtip(:,i-1,k)+this%vtip(:,i,k)+&
                       this%vtip(:,i-1,k_1)+this%vtip(:,i,k_1))*t
        end do
    end do
!//////////////////////////////////У����////////////////////////////////////////////
    call this%sub_vwake()
    this%grid=fgrid
    this%vgrid=0.5*(this%vgrid+fvgrid)
    this%tip=ftip
    this%vtip=0.5*(this%vtip+fvtip)
    do k=1,this%na
        k_1=k-1
        if(k==1) k_1=this%na
        do j=1,this%ns+1
            this%grid(:,1,j,k)=this%bpb(:,j,k)
            do i=2,this%nseg+1
                this%grid(:,i,j,k)=this%grid(:,i-1,j,k_1)+0.25*(this%vgrid(:,i-1,j,k)+this%vgrid(:,i,j,k)+&
                              this%vgrid(:,i-1,j,k_1)+this%vgrid(:,i,j,k_1))*t
                this%grid(:,i,j,k)=THIS%relax*this%grid(:,i,j,k)+(1-THIS%relax)*fgrid(:,i,j,k)
            end do
        end do
        this%tip(:,1,k)=(this%bpb(:,this%ns,k)+this%bpb(:,this%ns+1,k))/2.0
        do i=2,this%na*this%nf+1
            this%tip(:,i,k)=this%tip(:,i-1,k_1)+0.25*(this%vtip(:,i-1,k)+this%vtip(:,i,k)+&
                       this%vtip(:,i-1,k_1)+this%vtip(:,i,k_1))*t
            this%tip(:,i,k)=THIS%relax*this%tip(:,i,k)+(1-THIS%relax)*ftip(:,i,k)
        end do
        do i=this%na*this%nf+2,this%na*this%nw+1
            this%tip(1,i,k)=this%tip(1,i-this%na,k)+this%u*this%omg*this%r*this%na*t
            this%tip(2,i,k)=this%tip(2,i-this%na,k)
            this%tip(3,i,k)=this%tip(3,i-this%na,k)+(this%vdx+this%u*this%alphs)*this%omg*this%r*this%na*t
        end do
    end do
    !//////////////////////////////����RMS////////////////////////////////////////////
    this%rms_wake=0
    do k=1,this%na
        tmp=0
        do i=1,this%na*this%nw+1
            tmpm(:)=ftip(:,i,k)-this%tip(:,i,k)
            tmp=tmp+tmpm(1)**2+tmpm(2)**2+tmpm(3)**2
        end do
        this%rms_wake=this%rms_wake+tmp/this%na/(this%na*this%nw+1)
    end do
    !/////////////////////////////////////////////////////////////////////////////////
    
    
DEALLOCATE(fgrid,fvgrid, ftip,fvtip, &
            STAT=ERR)

IF(ERR.NE.0) THEN
    WRITE(*,999)
    STOP
END IF  
998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
end subroutine		
!���������ʾ-----------------------------------
subroutine sub_out(this,opmode)
implicit none
class(vortex) :: this
integer opmode,k,i,j
            real(rdt)::v(3),tmp(3)
            character(len=80)   filename,form1,form2
            select case(opmode)
            case(0) !----------------------------------------------��ʼ����
                call cpu_time(this%t_start)
                open(10,file='wake.dat')
                close(10,status="delete")
				open(10,file='wake2.dat')
                close(10,status="delete")
                open(10,file='tao.dat')
                close(10,status="delete")
                open(10,file='tao_tip.dat')
                close(10,status="delete")
                open(10,file='vbpc.dat')
                close(10,status="delete")
                open(10,file='vbpq.dat')
                close(10,status="delete")
                open(10,file='ct.dat')
                close(10,status="delete")
                open(10,file='a0s.dat')
                close(10,status="delete")
                open(10,file='a1s.dat')
                close(10,status="delete")
                open(10,file='b1s.dat')
                close(10,status="delete")
                open(10,file='fi7.dat')
                close(10,status="delete")
                open(10,file='cit1.dat')
                close(10,status="delete")
                open(10,file='cit2.dat')
                close(10,status="delete")
                open(10,file='load_blade.dat')
                close(10,status="delete")
                open(10,file='load_circle.dat')
                close(10,status="delete")
                open(10,file='vh.dat')
                close(10,status="delete")
                open(10,file='vorp.dat')
                close(10,status="delete")
                open(10,file='bladeshape1.dat')
                close(10,status="delete")
                open(10,file='bladeshape2.dat')
                close(10,status="delete")
                open(10,file='freewake_rms.dat')
                close(10,status="delete")            
            case(1) !----------------------------------------------��Ļ���
                write(form1,'(i4)') this%l_iter
                write(form2,'(e10.3)') this%rms_wake
                write(*,*) 'Iteration:',trim(form1),'    RMS:',trim(form2)
                open(11,file="freewake_rms.dat",position='append',action='write')
                write(11,*) trim(form1),trim(form2)
                close(11)
            case(2) !----------------------------------------------���������
                open(11,file="wake.dat",position='append',action='write')
                write(form1,'(i4.4)') this%l_iter
                DO J=1,THIS%nb
                    write(11,*)'variables = "x","y","z","u","v","w"'
                    write(11,*)'zone  t="1"','f=point i=',(this%nw*this%na+1),'j=1 k=1'
                    write(11,*)'STRANDID=1 ', 'SOLUTIONTIME=',trim(form1)
                    do i=1,this%na*this%nw+1
                        write(11,993) this%tip(:,i,1+this%na/THIS%NB*(J-1))/this%r ,this%vtip(:,i,1)
                    end do
                END DO
            case(22) !----------------------------------------------���������
                open(11,file="wake2.dat",position='append',action='write')
                write(form1,'(i4.4)') this%l_iter
                DO J=1,1
                    write(11,*)'variables = "x","y","z"'
                    write(11,*)'zone  t="1"','f=point i=',THIS%ns+1,'j=',THIS%nseg+1,' k=1'
                    write(11,*)'STRANDID=1 ', 'SOLUTIONTIME=',trim(form1)
					do k=1,THIS%nseg+1
						do i=1,THIS%ns+1
							write(11,993) this%grid(:,k,i,1+this%na/THIS%NB*(J-1))/this%r
						end do
					end do
                END DO
            case(3) !----------------------------------------------�������
            	open(10,file='tao.dat')
            	open(11,file='tao_tip.dat')
            	!write(10,'(2(I5))') this%l_iter
            	!write(form1,*) "(",this%ns,"F9.3)"
            	do j=1,this%ns
                    write(10,993) SQRT(this%bpc(1,j,1)**2+this%bpc(2,j,1)**2)/this%r,&
                    (this%tao(j+(i-1)*this%ns)/this%r**2/this%omg,i=1,this%na)
                end do    
            	do i=1,this%na
            	    write(11,993) this%tmax(i)
                end do
            	close(10)
            	close(11) 
            case(4) !----------------------------------------------�������
                open(11,file="vorp.dat")
                do i=1,this%na*this%nw+1
                    write(11,993) 360.0/this%na*(i-1),sqrt(this%tip(1,i,1)**2+this%tip(2,i,1)**2)/this%r,&
                                    this%tip(3,i,1)/this%r,sqrt(this%tip(1,i,this%na/2)**2+this%tip(2,i,this%na/2)**2)/this%r,&
                                    this%tip(3,i,this%na/2)/this%r
                end do
                close(11)
            case(5) !----------------------------------------------�����Ҷ����
                open(11,file="bladeshape1.dat")
                write(form1,'(i4.4)') this%l_iter
                do k=1,this%na
                    write(11,*)'variables = "x","y","z"'
                    write(11,*)'zone  t="',k,'" DATAPACKING=BLOCK i=',(this%ns+1),'j=2 k=1' 
                    !write(11,*)'STRANDID=1', 'SOLUTIONTIME=',trim(form1)
                    do j=1,3
                        write(11,993) (THIS%bpb(j,i,k)/this%r,i=1,this%ns+1),&
                         (THIS%bpq(j,i,k)/this%r,i=1,this%ns+1)
                    end do
                end do
                close(11)         
                open(11,file="bladeshape2.dat")
                write(form1,'(i4.4)') this%l_iter
                do k=1,this%na
                    write(11,*)'variables = "x","y","z"'
                    write(11,*)'zone  t="',k,'" DATAPACKING=BLOCK i=',(this%ns),'j=2 k=1' 
                    !write(11,*)'STRANDID=1', 'SOLUTIONTIME=',trim(form1)
                    do j=1,3
                        write(11,993) (THIS%bpc(j,i,k)/this%r,i=1,this%ns),&
                         ((THIS%bpb(j,i,k)+THIS%bpb(j,i+1,k))/2/this%r,i=1,this%ns)
                    end do
                end do
                close(11)   
            case(8) !---------------------------------------------�����Ҷ���յ��ٶ�
            	open(10,file='vbpq.dat',position='append')
            	open(11,file='vbpc.dat',position='append')
                write(form1,'(i4.4)') this%l_iter 
                write(10,*)'variables = "x","y"'
                write(10,*)'zone  t="',trim(form1),'" f=point i=',this%ns,'j=1 k=1'
                write(11,*)'variables = "x","y"'
                write(11,*)'zone  t="',trim(form1),'" f=point i=',this%ns,'j=1 k=1'
                k=5  !--------------------------------------------����յ��ٶȵ�λ��
                do i=1,this%ns
                    tmp(:)=(this%bpq(:,i,k)+this%bpq(:,i+1,k))/2
                    call this%sub_vi(k,tmp,v)
                    write(10,993) SQRT(tmp(1)*tmp(1)+tmp(2)*tmp(2))/this%r,v(3)/this%omg/this%r
                    tmp(:)=this%bpc(:,i,k)
                    call this%sub_vi(k,tmp,v)
                    write(11,993) SQRT(tmp(1)*tmp(1)+tmp(2)*tmp(2))/this%r,v(3)/this%omg/this%r
                end do
                close(10)
                close(11)
           end select
        993 format(1000f26.12)
        end subroutine	
!����ٶ�ֵ-----------------------------------
SUBROUTINE GET_VEL(THIS,PSI0,NP,XP,VP)
IMPLICIT NONE
CLASS(VORTEX) :: THIS
INTEGER,INTENT(IN)            :: NP
REAL(KIND=RDT),INTENT(IN)     :: XP(3,NP),PSI0
REAL(KIND=RDT),INTENT(OUT)    :: VP(3,NP)


REAL(KIND=RDT)                :: WORKWAKE(3,THIS%NA*THIS%NW+1,THIS%NB),WORKTAO(THIS%NB)
REAL(KIND=RDT)                :: PSI(THIS%NB)
REAL(KIND=RDT)                :: P(3),A(3),B(3),V(3),V1(3),AVG,RC
INTEGER                       :: I,J,KB,IP

PSI(1)  = PSI0

DO KB   = 2,THIS%NB
  PSI(KB)   = PSI(1) + (KB-1)*2.0*PI/THIS%NB
END DO

DO J  = 1,THIS%NA*THIS%NW+1
  DO I  = 1,3
!    CALL FFT_INTERPOLATION(THIS%NA,THIS%TIP(I,J,:),THIS%NB,PSI,WORKWAKE(I,J,:))
	CALL FFTINP(THIS%NA,THIS%NB,3*(THIS%NA*THIS%NW+1),THIS%TIP(:,:,:),WORKWAKE(:,J,:))
  END DO
END DO
!CALL FFT_INTERPOLATION(THIS%NA,THIS%TMAX,THIS%NB,PSI,WORKTAO)
CALL FFTINP(THIS%NA,THIS%NB,3*(THIS%NA*THIS%NW+1),THIS%TIP(:,:,:),WORKWAKE(:,J,:))

AVG   = SUM(THIS%TMAX)/THIS%NA
VP    = 0
!$OMP PARALLEL DO PRIVATE(IP,KB,I,P,A,B,RC,V,V1)
DO IP   = 1,NP
  P   = XP(:,IP)
  DO KB=1,THIS%NB
    V   = 0
    DO I  = 1,THIS%NA*THIS%NW
      A(:)=WORKWAKE(:,I,KB)
      B(:)=WORKWAKE(:,I+1,KB)
      CALL THIS%GET_RC(I,AVG,RC)
      CALL THIS%STL(P,A,B,RC,V1)
      V(:)=V(:)+V1(:)
    END DO
    VP(:,IP) = VP(:,IP) + V*WORKTAO(KB)
  END DO
END DO
!$OMP END PARALLEL DO
END SUBROUTINE
!����-----------------------------------
SUBROUTINE WPTEST(THIS,PSI0)
IMPLICIT NONE
CLASS(VORTEX) :: THIS
REAL(KIND=RDT),INTENT(IN)   :: PSI0
INTEGER,PARAMETER           :: NNN=101
REAL(KIND=RDT)              :: X(3,NNN),V(3,NNN)
INTEGER                     :: I,IU

DO I=0,NNN-1
  X(1,I+1) = 1.2*THIS%R*INT(I-NNN/2)/INT(NNN-1-NNN/2)
  X(2,I+1) = 0
  X(3,I+1) =-0.1*THIS%R
END DO
CALL THIS%GET_VEL(PSI0,NNN,X,V)
OPEN(IU,FILE="DOWNWASH-0.1.DAT")
DO I=1,NNN
  WRITE(IU,*) X(1,I)/THIS%R,V(3,I) !/(-VDX*OMG*R)
END DO
!DO I=1,NNN
!  CALL GET_VELOCITY(1,X(:,I),V(:,I))
!  WRITE(IU,*) X(1,I)/R,V(3,I) !/(-VDX*OMG*R)
!END DO
CLOSE(IU)

DO I=0,NNN-1
  X(1,I+1) = 1.2*THIS%R*INT(I-NNN/2)/INT(NNN-1-NNN/2)
  X(2,I+1) = 0
  X(3,I+1) =-0.4*THIS%R
END DO
CALL THIS%GET_VEL(PSI0,NNN,X,V)
OPEN(IU,FILE="DOWNWASH-0.4.DAT")
DO I=1,NNN
  WRITE(IU,*) X(1,I)/THIS%R,V(3,I) !/(-VDX*OMG*R)
END DO
CLOSE(IU)

END SUBROUTINE WPTEST
end module

MODULE MFWM
USE VORTEX_CLASS

    TYPE(VORTEX) :: VTX

END MODULE MFWM




!
