!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
subroutine ptix2p(y,m,c,n,cho)
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
implicit none
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
integer, intent(in):: n,cho
real(RDT), intent(in):: m(n,n),c(n,n)
real(RDT), intent(inout):: y(2*n)
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->        
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
    integer :: n2,err
    real(rdt),allocatable :: tp1(:),tp2(:)
    n2=2*n
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> main body: 
    if(cho.eq.0) return 
    ALLOCATE(tp1(n),tp2(n),&
            STAT=ERR)

	IF(ERR.NE.0) THEN
		WRITE(*,998)
		STOP
	END IF  
    

    tp1=matmul(m,y(n+1:n2))
    tp2=matmul(c,y(1:n))
    y(1:n)=y(1:n)
    y(n+1:n2)=tp1+1.0D0/2*tp2
    
    
	DEALLOCATE(tp1,tp2,&
            STAT=ERR)

	IF(ERR.NE.0) THEN
		WRITE(*,999)
		STOP
	END IF     
	 
998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
end subroutine


subroutine ptip2x(y,m,c,n,cho)
use GlobalDataFun
implicit none
integer, intent(in):: n,cho
real(RDT), intent(in):: m(n,n),c(n,n)
real(RDT), intent(inout):: y(2*n)

    integer :: n2,err
    real(rdt),allocatable :: invm(:,:),tp1(:),tp2(:)
    n2=2*n
    
    if(cho.eq.0) return
    ALLOCATE(tp1(n),tp2(n),invm(n,n),&
            STAT=ERR)

	IF(ERR.NE.0) THEN
		WRITE(*,998)
		STOP
	END IF  
	

    call mtxinv2(invm, 1, n, m)
    tp1=matmul(c,y(1:n))
    tp2=y(n+1:n2)-1.0/2*tp1 
 
    y(1:n)=y(1:n)
    y(n+1:n2)=matmul(invm,tp2)
    
    
	DEALLOCATE(tp1,tp2,invm,&
            STAT=ERR)

	IF(ERR.NE.0) THEN
		WRITE(*,999)
		STOP
	END IF     
	 
998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
end subroutine
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!->+++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> First-Order ode solver, dy/dt=f(y,t), explicit Runge-Kutta methods, one step.
!-> ch=1: runge-kutta-gill method.
!-> ch=0: classical runge-kutta method.
!->       coded by XiaoYu, 05162009.
!->       edited by rotor, 05192009.
subroutine ptisol(yout,y1,n,ta,ha,fa,y,h)
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
integer, intent(in):: n
real(RDT), intent(in):: ta(n,n),ha(n,n),fa(n),y(n),h
real(RDT), intent(out):: yout(n),y1(n)
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->      
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
    integer :: err
    real(rdt),allocatable :: yt(:),invha(:,:),tp1(:),tp2(:),tp3(:,:)

   
    ALLOCATE(tp1(n),tp2(n),tp3(n,n),&
                yt(n),invha(n,n),&
            STAT=ERR)

	IF(ERR.NE.0) THEN
		WRITE(*,998)
		STOP
	END IF

    call mtxinv2(invha, 1, n, ha)

    tp1=matmul(invha,fa)
    tp2=matmul(invha,fa)

    yt=matmul(ta,y)
    y1=matmul(ta,tp1)
    y1=y1-tp2
!    yt=matmul(ta,y+tp1)
!    yout=yt-tp2
    yout=yt+y1
    
!    tp1=matmul(ta,fa)
!    tp2=tp1+1
!
!    yt=matmul(ta,y)
!    y1=tp2/2*h
!!    y1=y1-tp2
!!    yt=matmul(ta,y+tp1)
!!    yout=yt-tp2
!    yout=yt+y1
    
    
    
	DEALLOCATE(tp1,tp2,yt,invha,&
            STAT=ERR)

	IF(ERR.NE.0) THEN
		WRITE(*,999)
		STOP
	END IF     
	 
998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
end subroutine
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
subroutine calpticoe(ha,fa,m,c,k,f,n,cho)
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
implicit none
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
integer, intent(in):: n,cho
real(RDT), intent(in):: m(n,n),c(n,n),k(n,n),f(n)
real(RDT), intent(out):: ha(2*n,2*n),fa(2*n)
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->        
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
    integer :: n2,err,i
    real(rdt),allocatable :: tp1(:,:),tp2(:,:),tp3(:,:)
    real(rdt),allocatable :: ma(:,:),mb(:,:),mc(:,:),md(:,:),invm(:,:)
   
    n2=2*n
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> main body:
    ALLOCATE(tp1(n,n),tp2(n,n),tp3(n,n), &
                ma(n,n),mb(n,n),mc(n,n),md(n,n),invm(n,n),&
            STAT=ERR)

	IF(ERR.NE.0) THEN
		WRITE(*,998)
		STOP
	END IF
	
    call mtxinv2(invm, 1, n, m)
	if(cho.eq.1) then
        tp1=matmul(invm,c)
        tp2=matmul(c,tp1)
        tp3=matmul(c,invm)
        ma=-1.0/2*tp1
        mb=1.0/4*tp2-k
        mc=-1.0/2*tp3
        md=invm
        
        ha(1:n,1:n)=ma(:,:)
        ha(1:n,n+1:n2)=md(:,:)
        ha(n+1:n2,1:n)=mb(:,:)
        ha(n+1:n2,n+1:n2)=mc(:,:)
        
        fa(1:n)=0.0
        fa(n+1:n2)=f(:)
    else if(cho.eq.0) then
        tp1=matmul(invm,c)
        tp2=matmul(invm,k)
        tp3=0.0
        do i=1,n
            tp3(i,i)=1.0
        end do
        ha(1:n,1:n)=0
        ha(1:n,n+1:n2)=tp3
        ha(n+1:n2,1:n)=-tp2
        ha(n+1:n2,n+1:n2)=-tp1
        
        fa(1:n)=0.0
        fa(n+1:n2)=matmul(invm,f)
    end if
        
    
	DEALLOCATE(tp1,tp2,tp3, &
                ma,mb,mc,md,invm,&
            STAT=ERR)

	IF(ERR.NE.0) THEN
		WRITE(*,999)
		STOP
	END IF     
	 
998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
end subroutine
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
subroutine calptitran(t,ha,n,h,ptic)
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
implicit none
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
integer, intent(in):: n,ptic
real(RDT), intent(in):: ha(n,n)
real(RDT), intent(in):: h
real(RDT), intent(out):: t(n,n)
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->        
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
    integer :: i,err
    real(rdt) :: pticp,tao
    real(rdt),allocatable :: tp1(:,:),tp2(:,:),tp3(:,:),tp4(:,:),tp5(:,:)
    real(rdt),allocatable :: ha1(:,:) ,ha2(:,:),ta(:,:),tt(:,:),ttt(:,:)
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> main body:  
    ALLOCATE(tp1(n,n),tp2(n,n),ha1(n,n) ,ha2(n,n),&
                ta(n,n),tt(n,n),ttt(n,n),&
            STAT=ERR)

	IF(ERR.NE.0) THEN
		WRITE(*,998)
		STOP
	END IF
       
    pticp=2.0D0**ptic
    tao=h/pticp
    ha1=ha*tao
    ha2=matmul(ha1,ha1)

    tp1=ha2/12.0+ha1/3.0
    do 20 i=1,n
        tp1(i,i)=tp1(i,i)+1.0
20  end do
    
    tp2=matmul(ha2,tp1)
    
    ta=ha1+1.0D0/2*tp2
!    write(*,*) ta(1,1)
!    pause
    do 30 i=1,ptic
        ttt=matmul(ta,ta)
        tt=ta*2.0+ttt
        ta=tt
30  end do
    t=ta
    do 40 i=1,n
        t(i,i)=t(i,i)+1.0
40  end do


	DEALLOCATE(tp1,tp2,ha1,ha2,&
                ta,tt,ttt,&
            STAT=ERR)

	IF(ERR.NE.0) THEN
		WRITE(*,999)
		STOP
	END IF     
	 
998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
 end subroutine
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  



!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
subroutine accupti(tn,yn,t,y,n)
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
implicit none
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
integer, intent(in):: n
real(RDT), intent(in):: t(n,n),y(n)
real(RDT), intent(inout):: tn(n,n),yn(n)
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->        
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> main body:  
    tn=matmul(t,tn)
    yn=matmul(t,yn)
    yn=yn+y


end subroutine accupti
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
