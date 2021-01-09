subroutine remove0(BRMTR,DM,NSTP)
use GlobalDataFun
implicit none
INTEGER,INTENT(IN) :: DM,NSTP
REAL(RDT),INTENT(INOUT) :: BRMTR(DM,NSTP)

character(80) filenam
integer err
integer m, n
integer lensav
real(rdt) omega
real(rdt), allocatable :: x(:, :), xot(:, :), t(:)
real(rdt), allocatable :: wsave(:)

      m=NSTP-1
      N=DM
      lensav=2*m+15
      
      ALLOCATE( x(m, n), xot(m,n), t(m), wsave(lensav), &
            STAT=ERR)
                  
	IF(ERR.NE.0) THEN
		WRITE(*, 998)
		STOP
	END IF
	
	X=TRANSPOSE(BRMTR(:,1:NSTP-1))

     call fftf(m, n, x, xot, lensav, wsave)
      
      x=xot
      
!      call drop(m,n,x,xot,0,0)
!      call extr(m,n,x,xot,1,10)
      
      call fftb(m, n, xot, x, lensav, wsave)
      
    BRMTR(:,1:NSTP-1)=TRANSPOSE(X)
    BRMTR(:,NSTP)=BRMTR(:,1)
    

      DEALLOCATE( x, xot, t, wsave, &
            STAT=ERR)

	IF(ERR.NE.0) THEN
		WRITE(*, 999)
		STOP
	END IF

998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
end subroutine remove0

subroutine fftwrapper()
use GlobalDataFun
implicit none

character(80) filenam
integer err
integer m, n,m2,j
integer lensav
real(rdt) omega
real(rdt), allocatable :: x(:, :), xot(:, :),  x2(:, :), xot2(:, :), t(:), t2(:)
real(rdt), allocatable :: wsave(:)

      m=100
      n=3
      m2=200
      lensav=2*m2+15
      
      ALLOCATE( x(m, n), xot(m,n),x2(m2, n), xot2(m2,n), t(m), t2(m2), wsave(lensav), &
            STAT=ERR)
                  
	IF(ERR.NE.0) THEN
		WRITE(*, 998)
		STOP
	END IF
	
	
      omega=100

	call tdfft(x,m,n,omega,t)
    do j=1,m
    write(112,993) t(j),x(j,1:n)
    end do
      call fftf(m, n, x, xot, lensav, wsave)
      x2=0.0D0
      x2(1:m,1:n)=xot(1:m,1:n)
      do j=1,m2
        t2(j)=(t(m)-t(1))/(m2-1)*(j-1)+t(1)
      end do
      call fftb(m2, n, x2, xot2, lensav, wsave)
    do j=1,m2
    write(113,993) t2(j),xot2(j,1:n)
    end do

!      filenam='fftout'
!      call outmatx(filenam, m, n, xot) 


!      DEALLOCATE( x, xot,x2,xot2, t, wsave, &
!            STAT=ERR)
!
!	IF(ERR.NE.0) THEN
!		WRITE(*, 999)
!		STOP
!	END IF
    stop
993 format(1000f26.12)
998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
end subroutine fftwrapper

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> to provide FFT test data.
!-> n+1 stations, and n devisions, provide the data of the first n stations.
!->       coded by XiaoYu, 05112009.
!->       edited by rotor, 05122009.
subroutine tdfft(x, n, m, omega, t)
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
integer, intent(in) :: n, m
real(rdt), intent(in) :: omega
real(rdt), intent(inout) :: x(n, m), t(n)

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->       dffti, dfftf       
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
integer i,j
real(rdt) f0, dt
	
      dt=datan(1.0D0)*8.0D0/omega/n
      !    f0=1.0D0/(atan(1.0D0)*8.0D0/omega)

      do i=1, n
            t(i)=(i-1)*dt
      end do
      do 5 i=1, m
        do j=1,n
            x(j, i)=1.2-i+(4.9-i)*COS(omega*1.0*t(j))-(30.0-i)*SIN(omega*1.0*t(j)) & 
                  -(0.0)*COS(omega*2.0*t(j))+(6.0-i)*SIN(omega*2.0*t(j)) &
                  -(4.0-i)*COS(omega*6.0*t(j))+(-3.0-i)*SIN(omega*6.0*t(j)) 
         end do
5    continue

end subroutine tdfft
!->+++++++++++++++++++++++++++++++++++++++++++++



!->+++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> FFT wrapper.
!-> n+1 stations, and n devisions, input the data of the first n stations.
!->       coded by XiaoYu, 05112009.
!->       edited by rotor, 05122009.
subroutine fftf(n, m, xin, xou, lensav, wsave)
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
integer :: n, m, lensav
real(rdt), intent(out):: xou(n, m)
real(rdt), intent(inout):: xin(n, m), wsave(lensav)

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->       dffti, dfftf       
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
integer i

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

	call dffti(n, wsave)
	
      do 10 i=1, m	
	      call dfftf(n, xin(:, i), wsave)
10  continue      
	xou=xin*2.0/n
	xou(1, :)=xou(1, :)/2.0D0
	xou(3:n:2, :)=-xou(3:n:2, :)
	
end subroutine fftf
!->+++++++++++++++++++++++++++++++++++++++++++++


subroutine fftb(n, m, xin, xou, lensav, wsave)
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
integer :: n, m, lensav
real(rdt), intent(out):: xou(n, m)
real(rdt), intent(inout):: xin(n, m), wsave(lensav)

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->       dffti, dfftf       
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
integer i

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

	call dffti(n, wsave)
	xin(1, :)=xin(1, :)*2.0D0
	xin(3:n:2, :)=-xin(3:n:2, :)
	xou=xin/2.0D0
      do 10 i=1, m	
	      call dfftb(n, xou(:, i), wsave)
10  continue      

end subroutine fftb
!->+++++++++++++++++++++++++++++++++++++++++++++

subroutine drop(m,n,x,xot,a,b)
use GlobalDataFun
implicit none
integer,intent(in) :: m,n,a,b
real(rdt),intent(in) :: x(m,n)
real(rdt),intent(out) :: xot(m,n)
 
    xot=x
    xot(a*2+1:b*2+1,:)=0.0D0
 
end
 !->+++++++++++++++++++++++++++++++++++++++++++++
subroutine extr(m,n,x,xot,a,b)
use GlobalDataFun
implicit none
integer,intent(in) :: m,n,a,b
real(rdt),intent(in) :: x(m,n)
real(rdt),intent(out) :: xot(m,n)
 
    xot=0.0D0
    xot(a*2+1:b*2+1,:)=x(a*2+1:b*2+1,:)
 
end    
!->+++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE FFTINP(n1,n2,m,xin,XOUT)
use GlobalDataFun
IMPLICIT NONE
INTEGER,INTENT(IN) :: M,N1,N2
REAL(RDT),INTENT(OUT) :: XOUT(M,N2)
REAL(RDT),INTENT(IN) :: xin(M,N1)
integer lensav,nst,nst2,dm,dmt
real(rdt) omega
real(rdt), allocatable :: x(:, :), x2(:,:),xot(:, :), xot2(:,:)
real(rdt), allocatable :: wsave(:)

    nst=N1
    nst2=N2
    dm=M
    lensav=2*(nst2-1)+15

    ALLOCATE( x(nst-1, dm), xot(nst-1, dm), x2(nst2-1, dm),xot2(nst2-1, dm), wsave(lensav))
              
    x=transpose(XIN(:, 1:nst-1))
    dmt=dm
    call fftf(nst-1, dmt, x, xot, lensav, wsave)
    
    x2=0.0D0
    !x2(1:nst-1,1:dm)=xot(1:nst-1,1:dm)
    !call extr(nst-1,dm,xot,x2,0,5)
    x2(0*2+1:5*2+1,1:dm)=xot(0*2+1:5*2+1,1:dm)
    call fftb(nst2-1, dmt, x2, xot2, lensav, wsave)
    XOUT(:,1:nst2-1)=transpose(xot2)
    XOUT(:,nst2)=XOUT(:,1)

    deALLOCATE( x,xot, x2,xot2, wsave)
END SUBROUTINE
!->+++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE FFTINP2(n1,n2,m,xin,XOUT)
use GlobalDataFun
IMPLICIT NONE
INTEGER,INTENT(IN) :: M,N1,N2
REAL(RDT),INTENT(OUT) :: XOUT(M,N2)
REAL(RDT),INTENT(IN) :: xin(M,N1)
integer lensav,nst,nst2,dm,dmt
real(rdt) omega
real(rdt), allocatable :: x(:, :), x2(:,:),xot(:, :), xot2(:,:)
real(rdt), allocatable :: wsave(:)

    nst=N1
    nst2=N2
    dm=M
    lensav=2*(nst2-1)+15

    ALLOCATE( x(nst-1, dm), xot(nst-1, dm), x2(nst2-1, dm),xot2(nst2-1, dm), wsave(lensav))
              
    x=transpose(XIN(:, 1:nst-1))
    dmt=dm
    call fftf(nst-1, dmt, x, xot, lensav, wsave)
    
    x2=0.0D0
    !x2(1:nst-1,1:dm)=xot(1:nst-1,1:dm)
    call extr(nst-1,dm,xot,x2,0,10)
    !x2(0*2+1:5*2+1,1:dm)=xot(0*2+1:5*2+1,1:dm)
    call fftb(nst2-1, dmt, x2, xot2, lensav, wsave)
    XOUT(:,1:nst2-1)=transpose(xot2)
    XOUT(:,nst2)=XOUT(:,1)

    deALLOCATE( x,xot, x2,xot2, wsave)
END SUBROUTINE

!FFT≤Â÷µ-----------------------------------
!SUBROUTINE FFT_INTERPOLATION(NSRC,SRC,NTAR,PSI,TAR)
!USE FFTW3
!use GlobalDataFun
!IMPLICIT NONE
!INTEGER,INTENT(IN)                :: NSRC,NTAR
!REAL(KIND=RDT),INTENT(IN)         :: SRC(NSRC),PSI(NTAR)
!REAL(KIND=RDT),INTENT(OUT)        :: TAR(NTAR)
!REAL(KIND=C_DOUBLE)               :: X(NSRC)
!COMPLEX(KIND=C_DOUBLE_COMPLEX)    :: Y(NSRC/2+1)
!TYPE(C_PTR)                       :: PLAN
!INTEGER                           :: I,J
!
!X = SRC
!
!PLAN = FFTW_PLAN_DFT_R2C_1D(INT(NSRC,KIND=C_INT),X,Y,FFTW_ESTIMATE)
!CALL FFTW_EXECUTE_DFT_R2C(PLAN,X,Y)
!CALL FFTW_DESTROY_PLAN(PLAN)
!
!Y = Y*2.0/NSRC
!
!DO I = 1,NTAR
!  TAR(I) = REAL(Y(1)) / 2.0
!  DO J=1,NSRC/2
!    TAR(I) = TAR(I) + REAL(Y(J+1))*COS(J*PSI(I)) - AIMAG(Y(J+1))*SIN(J*PSI(I))
!  END DO
!END DO
!
!END SUBROUTINE FFT_INTERPOLATION

