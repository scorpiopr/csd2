MODULE INFLOW_CLASS
USE INPUT_CLASS
IMPLICIT NONE
!->+++++++++++++++++++++++++++++++++++++++++++++
TYPE,PUBLIC :: INFLOW

    INTEGER  UNITYPE,LINTYPE
    REAL(RDT) KHLMDA,KFLMDA
    real(rdt) FXLMDA,FYLMDA
    
    CONTAINS

    PROCEDURE,PUBLIC ::  CONSTRUCT_INFLOW
    PROCEDURE,PUBLIC :: VEQ3
    PROCEDURE,PUBLIC ::linearinflow

    PROCEDURE,PRIVATE :: VEQ32
    PROCEDURE,PRIVATE :: jonsoninflow


END TYPE INFLOW
!->+++++++++++++++++++++++++++++++++++++++++++++
CONTAINS
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++PUBLIC++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE CONSTRUCT_INFLOW(THIS,INPU)
IMPLICIT NONE
CLASS(INFLOW) :: THIS
CLASS(INPUT) :: INPU

    THIS%UNITYPE=INPU%ifw%UNITYPE
    THIS%KHLMDA=INPU%ifw%KHLMDA
    THIS%KFLMDA=INPU%ifw%KFLMDA

    THIS%LINTYPE=INPU%ifw%LINTYPE
    THIS%FXLMDA=INPU%ifw%FXLMDA
    THIS%FYLMDA=INPU%ifw%FYLMDA

END SUBROUTINE

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> to chose the inflow model.
!-> coded by xiaoyu, 0603.2009.
!-> edited by rotor, 0616.2009. 
subroutine inflowmodel(lamil, lami, lamd, mu, psi, rb, vid, cho)
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
real(rdt), intent(out):: lamil
real(rdt), intent(in):: lami, lamd, mu, psi, rb
real(rdt), intent(in):: vid(3)
integer, intent(in):: cho

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
    
      IF(CHO.EQ.0) THEN
    
            lamil=lami
            
      ELSE IF(CHO.EQ.1) THEN
!            CALL dressinflow(SS,LL,HE,EL1,R,PSI)

            call dressinflow(lamil, lami, lamd, mu, rb, psi,1)
            
      ELSE IF(CHO.EQ.2) THEN
      
            call dynainflow(lamil, rb, psi, vid)       
      
      ELSE 
            write(*, *) ' error, inflow model.', 900011    
      END IF
 
end subroutine inflowmodel
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


subroutine linearinflow(THIS,lamids, lami, lamd, mu, rb, psi,J)
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
implicit none
CLASS(INFLOW) :: THIS
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
real(rdt), intent(out):: lamids
real(rdt), intent(in):: lami, lamd, mu, rb, psi
INTEGER :: J

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

    IF(THIS%LINTYPE.EQ.0) THEN
        CALL dressinflow(lamids, lami, lamd, mu, rb, psi,J)
    ELSE IF(THIS%LINTYPE.EQ.1) THEN
        CALL THIS%jonsoninflow(lamids, lami, lamd, mu, rb, psi)
    END IF


end subroutine linearinflow
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> the Dress inflow model.
!-> coded by xiaoyu, 0603.2009.
!-> edited by rotor, 0610.2009. 
subroutine dressinflow(lamids, lami, lamd, mu, rb, psi,J)
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
real(rdt), intent(out):: lamids
real(rdt), intent(in):: lami, lamd, mu, rb, psi

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->        
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
real(rdt) tp1, tp2, tp3
real(rdt) csx, snx, kx, ky
INTEGER :: J!上下旋翼判断:1上；-1下

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> main body:
    
      csx=COS(psi)
      snx=SIN(psi) 
     
      tp1=1.0D0-1.8d0*mu**2.0D0
      tp2=(lami+lamd)/mu
      tp3=SQRT(1.0+tp2**2.0)
      kx=4.0D0/3.0D0*(tp1*tp3-tp2)
      ky=-2.0D0*J*mu

      tp1=kx*csx
      tp2=ky*snx
      tp3=rb*(tp1+tp2)
      lamids=lami*(1.0D0+tp3)

end subroutine dressinflow
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> the Dress inflow model.
!-> coded by xiaoyu, 0603.2009.
!-> edited by rotor, 0610.2009. 
subroutine jonsoninflow(THIS,lamids, lami, lamd, mu, rb, psi)
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
implicit none
CLASS(INFLOW) :: THIS
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
real(rdt), intent(out):: lamids
real(rdt), intent(in):: lami, lamd, mu, rb, psi

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->        
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
real(rdt) tp1, tp2, tp3
real(rdt) csx, snx, kx, ky

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

      csx=COS(psi)
      snx=SIN(psi) 
     
      tp1=lami+lamd
      tp2=SQRT(TP1**2+MU**2)
      TP3=TP2+ABS(TP1)
      kx=THIS%FXLMDA*MU/TP3-THIS%FYLMDA*2*0.0
      ky=THIS%FXLMDA*0.0/TP3-THIS%FYLMDA*2*MU

      tp1=kx*csx
      tp2=ky*snx
      tp3=rb*(tp1+tp2)
      lamids=lami*(1.0D0+tp3)

end subroutine jonsoninflow
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> the dynmic inflow model.
!-> coded by xiaoyu, 0603.2009.
!-> edited by rotor, 0616.2009. 
subroutine dynainflow(lamids, rb, psi, vid)   
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
real(rdt), intent(out):: lamids
real(rdt), intent(in):: rb, psi
real(rdt), intent(in):: vid(3)

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->        
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
real(rdt) tp1, tp2, tp3
real(rdt) csx, snx, kx, ky

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> main body:
    
      csx=COS(psi)
      snx=SIN(psi) 
     
      kx=vid(3)
      ky=vid(2)

      tp1=kx*csx
      tp2=ky*snx
      tp3=rb*(tp1+tp2)
      lamids=vid(1)+tp3

end subroutine  dynainflow
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


SUBROUTINE VEQ3(THIS,AS,CTC,MU,LAMI,LAMD)
implicit none
CLASS(INFLOW) :: THIS

REAL(RDT), INTENT(INOUT):: MU,AS,CTC
REAL(RDT), INTENT(OUT):: LAMI,LAMD

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++

IF(THIS%UNITYPE.EQ.0) THEN
        CALL VEQ31(AS,CTC,MU,LAMI,LAMD)
    ELSE IF(THIS%UNITYPE.EQ.1) THEN
        CALL THIS%VEQ32(AS,CTC,MU,LAMI,LAMD)
    END IF
    
END SUBROUTINE 

SUBROUTINE VEQ32(THIS,AS,CTC,MU,LAMI,LAMD)
implicit none
CLASS(INFLOW) :: THIS

REAL(RDT), INTENT(IN):: MU,AS,CTC
REAL(RDT), INTENT(OUT):: LAMI,LAMD

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++
REAL(RDT) CT,LAMII,EPS
INTEGER NI
REAL(RDT) COE,LAM,LAMH

    EPS=1.0D-10
    CT=CTC/2
    LAMD=MU*SIN(AS)/COS(AS)
    LAMH=SQRT(CT/2)
    NI=1

    IF( MU**2+(2*LAMD+3*LAMH)**2.LT.LAMH**2 ) THEN
        LAMI=THIS%KHLMDA*LAMD*((0.373*LAMD**2+0.598*MU**2)/LAMH**2-1.991)
        WRITE(*,*) LAMI,'ENTER VORTEX RING'
        PAUSE
    ELSE 
        LAMI=CT/2/SQRT(CT/2/THIS%KHLMDA**2+(MU/THIS%KFLMDA)**2)           
        DO 5
            LAM=LAMI+LAMD
            COE=(LAM**2/THIS%KHLMDA**4+MU**2/THIS%KFLMDA**2)**1.5
            
            LAMII=CT/2*(LAM*(LAM+LAMI)/THIS%KHLMDA**4+(MU/THIS%KFLMDA)**2)/ &
                  (COE+CT/2*LAM/THIS%KHLMDA**4)
            IF( ABS(LAMII-LAMI) .LT. EPS ) THEN
                EXIT
            END IF
            LAMI=LAMII
            NI=NI+1
            IF(NI.GT.1000) THEN
                WRITE(*,*) NI,CTC,'ERROR IN VEQ3'
                STOP
            END IF
5		    CONTINUE
    END IF
    
END SUBROUTINE 


SUBROUTINE VEQ31(AS,CTC,MU,LAMI,LAMD)
implicit none

REAL(RDT), INTENT(IN):: MU,AS,CTC
REAL(RDT), INTENT(OUT):: LAMI,LAMD

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++
REAL(RDT) CT,LAMII,EPS
INTEGER NI

    EPS=1.0D-10
    CT=CTC/2
    LAMD=MU*SIN(AS)/COS(AS)
    NI=1

    IF( MU .LT. 1.0D-6 ) THEN!悬停状态
        LAMI=1.0*SQRT(CT/2)
    ELSE !前飞状态
        LAMI=CT/2/MU
        LAMI=CT/2/SQRT(MU**2+CT/2)           
        DO 5
            LAMII=CT/2/SQRT(MU**2+(LAMD+LAMI)**2)
            !LAMII=(LAMD+CT*(MU**2+2*(LAMD+LAMI)**2)/2/(MU**2+(LAMD+LAMI)**2)**1.5)/(1+CT*(LAMD+LAMI)/2/(MU**2+(LAMD+LAMI)**2)**1.5)!祁浩天论文前飞均匀入流模型
            IF( ABS(LAMII-LAMI) .LT. EPS ) THEN
                EXIT
            END IF
            LAMI=LAMII
            NI=NI+1
            IF(NI.GT.1000) THEN
                WRITE(*,*) NI,CTC,'ERROR IN VEQ3'
                EXIT
            END IF
5		    CONTINUE
    END IF
END SUBROUTINE 

END
