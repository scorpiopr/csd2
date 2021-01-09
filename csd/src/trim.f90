MODULE HTRIM_CLASS
USE INPUT_CLASS
IMPLICIT NONE
!->+++++++++++++++++++++++++++++++++++++++++++++
TYPE,PUBLIC :: HTRIM

    INTEGER :: MKTRIM
    INTEGER :: isinput
    REAL(RDT) :: rlxr,rlxv
    REAL(RDT) :: fcdfpdka
    REAL(RDT) :: rcgb(3),racb(3), fchbtg(3), mthbtg(3)

    REAL(RDT) :: CW
    REAL(RDT) :: lockn,THTW

    CONTAINS

    PROCEDURE,PUBLIC :: CONSTRUCT_HTRIM
    PROCEDURE,PUBLIC :: trimest

END TYPE HTRIM
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++PUBLIC++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
CONTAINS
SUBROUTINE CONSTRUCT_HTRIM(THIS,INPU)
IMPLICIT NONE
CLASS(HTRIM) :: THIS
TYPE(INPUT) :: INPU

    THIS%THTW=0.0
    THIS%lockn=INPU%IPT%lockn
    THIS%cw=INPU%TRM%cw
    THIS%MKTRIM=INPU%IPT%MKTRIM
    THIS%isinput=INPU%TRM%isinput
    THIS%rlxr=INPU%TRM%rlxr
    THIS%rlxv=INPU%TRM%rlxv
    THIS%fcdfpdka=INPU%TRM%fcdfpdka
    THIS%rcgb=INPU%TRM%rcgb
    THIS%racb=INPU%TRM%racb
    THIS%fchbtg=INPU%TRM%fchbtg
    THIS%mthbtg=INPU%TRM%mthbtg

END SUBROUTINE
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> to estimate the controlling pitch angles.
!-> W. J., H. T., pp. 193.
subroutine trimest(THIS,th75, th1s, th1c, mu, ct, ladtpp, bt1c, bt1s, siga, alcr)
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
implicit none
CLASS(HTRIM) :: THIS
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
real(RDT), intent(in):: mu, ct, ladtpp, bt1c, bt1s, siga, alcr
real(RDT), intent(out):: th75, th1s, th1c

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
real(RDT) x2, x4, x6, xs, xc, xd, bt0

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

      !>bt1s=0.0D0
      !>bt1c=0.0D0 
      !>ladtpp=lad 
        
      x2=mu*mu
      x4=x2*x2
      x6=x4*x2
      xs=1.0D0-x2+9.0D0*x4/4.0D0
      xc=6.0D0*ct/siga/alcr 
      xd=xc+0.375d0*x2*this%thtw
       
      bt0=(1.0D0-19.0D0*x2/18.0D0+1.5d0*x4)*xc
      bt0=bt0+(0.05d0+29.0D0*x2/120.0D0-0.2d0*x4+0.375d0*x6)*this%thtw
      bt0=bt0+(1.0D0/6.0D0-7.0D0*x2/12.0D0+0.25d0*x4)*ladtpp
      bt0=bt0*this%lockn/xs/8.0D0
      th1c=bt1s+4.0D0*mu*bt0/(1.0D0+0.5d0*x2)/3.0D0  
      th1s=-bt1c-(8.0D0*mu*xd/3.0D0+2.0D0*mu*ladtpp*(1.0D0-1.5d0*x2))/xs
      th75=((1.0D0+1.5d0*x2)*xd+1.5d0*ladtpp*(1.0D0-0.5d0*x2))/xs

end  subroutine trimest
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
END MODULE
