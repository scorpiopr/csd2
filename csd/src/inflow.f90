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
    ELSE IF(THIS%UNITYPE.EQ.2) THEN
        CALL VEQ33(AS,CTC,MU,LAMI,LAMD)
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

!简单定点迭代法求诱导入流比
SUBROUTINE VEQ31(AS,CTC,MU,LAMI,LAMD)
implicit none

REAL(RDT), INTENT(IN):: MU,AS,CTC
REAL(RDT), INTENT(OUT):: LAMI,LAMD

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++
REAL(RDT) :: EPSMIN=100.0
REAL(RDT) CT,LAMII,EPS
INTEGER NI

    EPS=1.0D-10
    CT=CTC/2
    LAMD=MU*SIN(AS)/COS(AS)
    NI=1

    IF( MU .LT. 1.0D-6 ) THEN!悬停状态
        LAMI=1.0*SQRT(CT/2)
    ELSE !前飞状态
        LAMI=CT/2/SQRT(MU**2+CT/2)           
        DO 5
            LAMII=CT/2/SQRT(MU**2+(LAMD+LAMI)**2)!简单定点迭代法求诱导入流比，公式参考Principles of Helicopter Aerodynamics,J.Gordon.Leishman,2ed,P133(95)中的式2.132
            !LAMII=(LAMD+CT*(MU**2+2*(LAMD+LAMI)**2)/2/(MU**2+(LAMD+LAMI)**2)**1.5)/(1+CT*(LAMD+LAMI)/2/(MU**2+(LAMD+LAMI)**2)**1.5)!祁浩天论文前飞均匀入流模型
            
            IF(EPSMIN.GT.ABS(LAMII-LAMI))THEN
                EPSMIN=ABS(LAMII-LAMI)
            END IF
            
            IF( ABS(LAMII-LAMI) .LT. EPS ) THEN
                EXIT
            END IF
            LAMI=LAMII
            NI=NI+1
            IF(NI.GT.1000) THEN
                WRITE(*,*) EPSMIN
                WRITE(*,*) NI,CTC,'ERROR IN VEQ3'
                EXIT
            END IF
5		    CONTINUE
    END IF
END SUBROUTINE 

!牛顿迭代法求桨盘均匀诱导入流比
SUBROUTINE VEQ33(AS,CTC,MIU_XY,TEMP2,MIU_Z)
    IMPLICIT NONE
    
    REAL(RDT), INTENT(IN):: MIU_XY,AS,CTC
    REAL(RDT), INTENT(OUT):: TEMP2,MIU_Z
    
    REAL(RDT) :: CT,TEMP1,FUNC_LAMBDAi,DFUNC_LAMBDAi
    INTEGER :: J=1
    
    CT=CTC/2
    MIU_Z=MIU_XY*SIN(AS)/COS(AS)
    
    IF( MIU_XY .LT. 1.0D-6 ) THEN!悬停状态
        TEMP2=DSQRT(CT/2)
    ELSE !前飞状态
        TEMP1=DSQRT(CT/2)
        
        !迭代进一步求初始均匀诱导速度
        DO J = 1, 10000
            IF(DABS(TEMP1) < 1.D-12) TEMP1 = 1.D-12
            !FUNC_LAMBDA = DSQRT((MIU_Z + TEMP1)**2 + MIU_XY**2) + CT * 5.D-1 / TEMP1       !公式参考Principles of Helicopter Aerodynamics,J.Gordon.Leishman,2ed,P133(95)起的式2.125、2.126、2.135、2.136
            !DFUNC_LAMBDA = (MIU_Z + TEMP1) / DSQRT((MIU_Z + TEMP1)**2 + MIU_XY**2) - CT * 5.D-1 / TEMP1**2  !但要注意：1)书中公式为桨盘总无量纲入流的函数，而本程序中为诱导速度无量纲入流的函数;2)书中公式λ默认从桨盘上方流向下方为正，在本程序直接与程序坐标系结合，为沿局部坐标系Z轴负方向流动。因此函数中诱导速度项符号与书中相反。
            FUNC_LAMBDAi = DSQRT((MIU_Z + TEMP1)**2 + MIU_XY**2) - CT * 5.D-1 / TEMP1!桨盘诱导速度向下为正
            DFUNC_LAMBDAi = (MIU_Z + TEMP1) / DSQRT((MIU_Z + TEMP1)**2 + MIU_XY**2) + CT * 5.D-1 / TEMP1**2
            IF(DABS(DFUNC_LAMBDAi) < 1.D-12) THEN
                WRITE(*,*) "诱导速度初始化导数项为零!"
                PAUSE
                PAUSE
                STOP
            END IF
            TEMP2 = TEMP1 - FUNC_LAMBDAi / DFUNC_LAMBDAi               !牛顿迭代法  x(n+1) = x(n) - [f(x(n))/f'(x(n))]
            IF(DABS(TEMP2 - TEMP1) < 1.D-8) THEN
                EXIT
            ELSE IF(J == 9999) THEN
                WRITE(*,*) "诱导速度初始化未收敛!"
                PAUSE
                PAUSE
                STOP
            END IF
            TEMP1 =TEMP2
        END DO
   END IF
END SUBROUTINE

END    
