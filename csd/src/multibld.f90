
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
subroutine fmrbsgl2hub(fmhub, fmhubhw, dm, nst, nb, nhw, fmsgl) 
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
integer, intent(in):: dm, nst, nb, nhw
real(rdt), intent(in):: fmsgl(dm, nst)
real(rdt), intent(out):: fmhub(dm, nst), fmhubhw(nst-1, dm)

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->        
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
character(80) filenam
integer err
integer lensav
real(rdt), parameter:: oned=1.0D0
real(rdt) psit, pai
real(rdt), allocatable :: fmmul(:, :, :)
real(rdt), allocatable :: x(:, :)
real(rdt), allocatable :: wsave(:), xot(:, :)

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

      pai=ACOS(-oned)
      psit=2.0D0*pai

      lensav=2*(nst-1)+15
      
      ALLOCATE( fmmul(dm, nb, nst), &
            x(nst-1, dm), &
            wsave(lensav), xot(nst-1, dm), & 
            STAT=ERR)
                  
	IF(ERR.NE.0) THEN
		WRITE(*, 998)
		STOP
	END IF
	

      call fmsgl2mul(fmmul, dm, nb, nst, nhw, psit, fmsgl) 
     
      call fmrb2hb(fmhub, dm, nb, nst, fmmul)  
 
      x=transpose(fmhub(:, 1:nst-1))
      call fftf(nst-1, dm, x, xot, lensav, wsave)

        fmhubhw=xot
      
      DEALLOCATE( fmmul, &
            x, &
            STAT=ERR)

	IF(ERR.NE.0) THEN
		WRITE(*, 999)
		STOP
	END IF

998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
end subroutine fmrbsgl2hub
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> psit = omga * T = 2 * Pi.
subroutine fmsgl2mul(fmmul, dm, nb, nst, nhw, psit, fmsgl) 
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
integer, intent(in):: dm, nb, nst, nhw
real(rdt) psit
real(rdt), intent(in):: fmsgl(dm, nst)
real(rdt), intent(out):: fmmul(dm, nb, nst)

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->        
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
character(80) filenam
integer err
integer lensav
integer ist, ibd, ihw
real(rdt) pai, dpibd, dpist
real(rdt) psibd, psi0
real(rdt) psbd(nb)
real(rdt), allocatable :: x(:, :), xot(:, :), fmsglfft(:, :)
real(rdt), allocatable :: wsave(:)
integer dmt
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

      lensav=2*(nst-1)+15
      
      ALLOCATE( x(nst-1, dm), xot(nst-1, dm), wsave(lensav), &
            fmsglfft(dm, nst-1), & 
            STAT=ERR)
                  
	IF(ERR.NE.0) THEN
		WRITE(*, 998)
		STOP
	END IF
	
      dpist=psit/(nst-1) 	
      dpibd=psit/nb
      
	x=transpose(fmsgl(:, 1:nst-1))
	dmt=dm
	
      call fftf(nst-1, dmt, x, xot, lensav, wsave)
      fmsglfft=transpose(xot)

      
      do 5 ibd=1, nb
            psbd(ibd)=dble(ibd-1)*dpibd 
5    continue

      do 10 ist=1, nst
            psi0=(ist-1)*dpist  
            do 20 ibd=1, nb
                  psibd=psi0+psbd(ibd)
                  fmmul(:, ibd, ist)=fmsglfft(:, 1)
                  do 30 ihw=1, nhw 
                        fmmul(:, ibd, ist)=fmmul(:, ibd, ist)+fmsglfft(:, 2*ihw)*cos(ihw*psibd)&      
                              +fmsglfft(:, 2*ihw+1)*sin(ihw*psibd)
30              continue
20        continue
10  continue


      DEALLOCATE( x, xot, wsave, &
            fmsglfft, & 
            STAT=ERR)

	IF(ERR.NE.0) THEN
		WRITE(*, 999)
		STOP
	END IF

998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
end subroutine fmsgl2mul
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> to transfer the forces and moments from rotating coordinates system to no-rotating
!-> hub coordinates systems.
subroutine fmrb2hb(fmhb, dm, nb, nst, fmrb) 
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
integer, intent(in):: dm, nb, nst
real(rdt), intent(in):: fmrb(dm, nb, nst)
real(rdt), intent(out):: fmhb(dm, nst)

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->        
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
integer ist, ibd
real(rdt), parameter:: oned=1.0D0
real(rdt) pai, dpibd, dpist
real(rdt) psibd, psi0
real(rdt) psbd(nb)
real(rdt) ttmt(dm, dm), vthb(dm), vtrb(dm)

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

      pai=ACOS(-oned)
      dpibd=2.0D0*pai/nb
      dpist=2.0D0*pai/(nst-1) 

      do 5 ibd=1, nb
            psbd(ibd)=dble(ibd-1)*dpibd 
5    continue
  
      !>   dm==3
      fmhb=0.0D0
      do 10 ist=1, nst
            psi0=(ist-1)*dpist  
            do 20 ibd=1, nb
                  psibd=psbd(ibd)+psi0
                  CALL TXYZ(ttmt, 3, -psibd)
                  vtrb(:)=fmrb(:, ibd, ist)
                  CALL MVMU2(vthb, dm, dm, ttmt, vtrb)
                  fmhb(:, ist)=fmhb(:, ist)+vthb(:)
20        continue             
10  continue             

end subroutine fmrb2hb
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
