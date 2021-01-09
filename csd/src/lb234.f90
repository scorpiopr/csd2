MODULE LB_CLASS
USE LB1_CLASS
USE LB234_CLASS
IMPLICIT NONE
!->+++++++++++++++++++++++++++++++++++++++++++++
TYPE,PUBLIC :: LB

    TYPE(LB1) :: LB1
    TYPE(LB2) :: LB2
    TYPE(LB3) :: LB3
    TYPE(LB4) :: LB4

    CONTAINS

    PROCEDURE,PUBLIC :: GETPARA
    PROCEDURE,PUBLIC :: LBPRE
    PROCEDURE,PUBLIC :: LBDELOCA
    PROCEDURE,PUBLIC :: LBINI
    PROCEDURE,PUBLIC :: LBMODEL
    PROCEDURE,PUBLIC :: TESTLBMODEL3
    PROCEDURE,PUBLIC :: TESTLBMODEL4
    
END TYPE LB
CONTAINS

SUBROUTINE LBPRE(THIS,J,CHO)
IMPLICIT NONE
CLASS(LB) :: THIS
INTEGER :: J,CHO
    IF(CHO.EQ.0) THEN
        THIS%LB1%AFASTD=THIS%LB2%AFASTDT(J)
        THIS%LB1%QSTD=THIS%LB2%QSTDT(J)
        THIS%LB1%CNCFSTD=THIS%LB2%CNCFSTDT(J)
        THIS%LB1%FPMSTD=THIS%LB2%FPMSTDT(J)
        THIS%LB1%FPNSTD=THIS%LB2%FPNSTDT(J)
        THIS%LB1%CVSTD=THIS%LB2%CVSTDT(J)
        
        THIS%LB1%XC=THIS%LB2%XCT(J)
        THIS%LB1%YC=THIS%LB2%YCT(J)
        THIS%LB1%ZC=THIS%LB2%ZCT(J)
        
        THIS%LB1%D1=THIS%LB2%D1T(J)
        THIS%LB1%D2=THIS%LB2%D2T(J)
        THIS%LB1%D3=THIS%LB2%D3T(J)
        THIS%LB1%D4=THIS%LB2%D4T(J)
        THIS%LB1%D5=THIS%LB2%D5T(J)
        THIS%LB1%D6=THIS%LB2%D6T(J)
        THIS%LB1%D7M=THIS%LB2%D7MT(J)
        THIS%LB1%D7N=THIS%LB2%D7NT(J)
            
        THIS%LB1%DELTA_XX1STD=THIS%LB2%DELTA_XX1STDT(J)
        THIS%LB1%DELTA_XX2STD=THIS%LB2%DELTA_XX2STDT(J)
        THIS%LB1%DELTA_XX3STD=THIS%LB2%DELTA_XX3STDT(J)
        THIS%LB1%DELTA_XX4STD=THIS%LB2%DELTA_XX4STDT(J)
        
        THIS%LB1%CNV=THIS%LB2%CNVT(J)
        THIS%LB1%TAOV=THIS%LB2%TAOVT(J)
    ELSE IF(CHO.EQ.1) THEN
        THIS%LB2%AFASTDT(J)=THIS%LB1%AFASTD
        THIS%LB2%QSTDT(J)=THIS%LB1%QSTD
        THIS%LB2%CNCFSTDT(J)=THIS%LB1%CNCFSTD
        THIS%LB2%FPMSTDT(J)=THIS%LB1%FPMSTD
        THIS%LB2%FPNSTDT(J)=THIS%LB1%FPNSTD
        THIS%LB2%CVSTDT(J)=THIS%LB1%CVSTD
        
        THIS%LB2%XCT(J)=THIS%LB1%XC
        THIS%LB2%YCT(J)=THIS%LB1%YC
        THIS%LB2%ZCT(J)=THIS%LB1%ZC
        
        THIS%LB2%D1T(J)=THIS%LB1%D1
        THIS%LB2%D2T(J)=THIS%LB1%D2
        THIS%LB2%D3T(J)=THIS%LB1%D3
        THIS%LB2%D4T(J)=THIS%LB1%D4
        THIS%LB2%D5T(J)=THIS%LB1%D5
        THIS%LB2%D6T(J)=THIS%LB1%D6
        THIS%LB2%D7MT(J)=THIS%LB1%D7M
        THIS%LB2%D7NT(J)=THIS%LB1%D7N
          
        THIS%LB2%DELTA_XX1STDT(J)=THIS%LB1%DELTA_XX1STD
        THIS%LB2%DELTA_XX2STDT(J)=THIS%LB1%DELTA_XX2STD
        THIS%LB2%DELTA_XX3STDT(J)=THIS%LB1%DELTA_XX3STD
        THIS%LB2%DELTA_XX4STDT(J)=THIS%LB1%DELTA_XX4STD
        
        THIS%LB2%CNVT(J)=THIS%LB1%CNV
        THIS%LB2%TAOVT(J)=THIS%LB1%TAOV
     ELSE IF(CHO.EQ.3) THEN
        THIS%LB1%AFASTD=THIS%LB4%AFASTDO
        THIS%LB1%QSTD=THIS%LB4%QSTDO
        THIS%LB1%CNCFSTD=THIS%LB4%CNCFSTDO
        THIS%LB1%FPMSTD=THIS%LB4%FPMSTDO
        THIS%LB1%FPNSTD=THIS%LB4%FPNSTDO
        THIS%LB1%CVSTD=THIS%LB4%CVSTDO
        
        THIS%LB1%XC=THIS%LB4%XCO
        THIS%LB1%YC=THIS%LB4%YCO
        THIS%LB1%ZC=THIS%LB4%ZCO
        
        THIS%LB1%D1=THIS%LB4%D1O
        THIS%LB1%D2=THIS%LB4%D2O
        THIS%LB1%D3=THIS%LB4%D3O
        THIS%LB1%D4=THIS%LB4%D4O
        THIS%LB1%D5=THIS%LB4%D5O
        THIS%LB1%D6=THIS%LB4%D6O
        THIS%LB1%D7M=THIS%LB4%D7MO
        THIS%LB1%D7N=THIS%LB4%D7NO
            
        THIS%LB1%DELTA_XX1STD=THIS%LB4%DELTA_XX1STDO
        THIS%LB1%DELTA_XX2STD=THIS%LB4%DELTA_XX2STDO
        THIS%LB1%DELTA_XX3STD=THIS%LB4%DELTA_XX3STDO
        THIS%LB1%DELTA_XX4STD=THIS%LB4%DELTA_XX4STDO
        
        THIS%LB1%CNV=THIS%LB4%CNVO
        THIS%LB1%TAOV=THIS%LB4%TAOVO
    ELSE IF(CHO.EQ.4) THEN
        THIS%LB4%AFASTDO=THIS%LB1%AFASTD
        THIS%LB4%QSTDO=THIS%LB1%QSTD
        THIS%LB4%CNCFSTDO=THIS%LB1%CNCFSTD
        THIS%LB4%FPMSTDO=THIS%LB1%FPMSTD
        THIS%LB4%FPNSTDO=THIS%LB1%FPNSTD
        THIS%LB4%CVSTDO=THIS%LB1%CVSTD
        
        THIS%LB4%XCO=THIS%LB1%XC
        THIS%LB4%YCO=THIS%LB1%YC
        THIS%LB4%ZCO=THIS%LB1%ZC
        
        THIS%LB4%D1O=THIS%LB1%D1
        THIS%LB4%D2O=THIS%LB1%D2
        THIS%LB4%D3O=THIS%LB1%D3
        THIS%LB4%D4O=THIS%LB1%D4
        THIS%LB4%D5O=THIS%LB1%D5
        THIS%LB4%D6O=THIS%LB1%D6
        THIS%LB4%D7MO=THIS%LB1%D7M
        THIS%LB4%D7NO=THIS%LB1%D7N
          
        THIS%LB4%DELTA_XX1STDO=THIS%LB1%DELTA_XX1STD
        THIS%LB4%DELTA_XX2STDO=THIS%LB1%DELTA_XX2STD
        THIS%LB4%DELTA_XX3STDO=THIS%LB1%DELTA_XX3STD
        THIS%LB4%DELTA_XX4STDO=THIS%LB1%DELTA_XX4STD
        
        THIS%LB4%CNVO=THIS%LB1%CNV
        THIS%LB4%TAOVO=THIS%LB1%TAOV
    END IF
END SUBROUTINE


SUBROUTINE LBDELOCA(THIS)
IMPLICIT NONE
CLASS(LB) :: THIS

    CALL THIS%LB2%LBDE1()
    CALL THIS%LB3%LBDE2()

END SUBROUTINE

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> THE PURPOSE:
SUBROUTINE GETPARA(THIS,ALCSRF,AFA1,DAFA1,SS1,SS2,K0,K1,K2,CD0,DF,CN1,TP,TF,TV,TAOV1,AFA0,CM0,MA)
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> MODULES USED:
IMPLICIT NONE
CLASS(LB) :: THIS
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> DUMMY ARGUMENTS / INTERFACE: 
real(rdt),INTENT(IN) ::MA
real(rdt),INTENT(OUT) :: ALCSRF,AFA1,DAFA1,SS1,SS2,K0,K1,K2,CD0,DF,CN1,TP,TF,TV,TAOV1,AFA0,CM0
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> DUMMY ARGUMENTS / INTERFACE: 
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->       
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->       READDATA1, READDATA2
!->       EIGANS, HOVER
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
real(rdt) :: VMA(7)
real(rdt) :: VALCSRF(7),VAFA1(7),DVAFA1(7),VSS1(7),VSS2(7),VK0(7),VK1(7),VK2(7),VCD0(7)
real(rdt) :: VDF(7),VCN1(7),VTP(7),VTF(7),VTV(7),VTAOV1(7)
real(rdt) :: PI,SS
INTEGER :: I,J,IJ,imk

    PI=ACOS(0.0D0)*2
!for NACA0012
    DATA VMA/0.3,0.4,0.5,0.6,0.7,0.75,0.8/
    DATA VALCSRF /0.108,0.113,0.117,0.127,0.154,0.175,0.216/
    DATA VAFA1 /15.25,12.5,10.5,8.5,5.6,3.5,0.7/
    DATA DVAFA1 /2.1,2.0,1.45,1.0,0.8,0.2,0.1/
    DATA VSS2 /3.0,3.25,3.5,4.0,4.5,3.5,0.7/
    DATA VSS1 /2.3,1.6,1.2,0.7,0.5,0.8,0.18/
    DATA VK0 /0.0025,0.006,0.02,0.038,0.030,0.001,-0.01/
    DATA VK1 /-0.135,-0.135,-0.125,-0.12,-0.09,-0.13,0.02/
    DATA VK2 /0.04,0.05,0.04,0.04,0.05,-0.02,-0.01/
    DATA VCD0 /0.0085,0.008,0.0077,0.0078,0.0078,0.0079,0.0114/
    DATA VDF /8.0,7.75,6.2,6.0,5.9,5.5,4.0/
    DATA VCN1 /1.45,1.2,1.05,0.92,0.68,0.5,0.18/
    
    DATA VTP /1.7,1.8,2.0,2.5,3.0,3.3,4.3/
    DATA VTF /3.0,2.5,2.2,2.0,2.0,2.0,2.0/
    DATA VTV /6.0,6.0,6.0,6.0,6.0,6.0,4.0/
    DATA VTAOV1 /7.0,9.0,9.0,9.0,9.0,9.0,9.0/
    CM0=0.0D0
    AFA0=0.0D0
!----------------------------------------------------分割线----------------------------
!for SC1095
    DATA VMA/0.3,0.4,0.5,0.6,0.7,0.75,0.8/
    DATA VALCSRF /0.13,0.1203, 0.1287, 0.1397, 0.144, 0.154, 0.175/
    DATA VAFA1 /15.25,12.5,10.5,8.5,5.6,3.5,0.7/
    DATA DVAFA1 /2.1,2.0,1.45,1.0,0.8,0.2,0.1/
    DATA VSS2 /3.0,3.25,3.5,4.0,4.5,3.5,0.7/
    DATA VSS1 /2.3,1.6,1.2,0.7,0.5,0.8,0.18/
    DATA VK0 /0.0025,0.006,0.02,0.038,0.030,0.001,-0.01/
    DATA VK1 /-0.135,-0.135,-0.125,-0.12,-0.09,-0.13,0.02/
    DATA VK2 /0.04,0.05,0.04,0.04,0.05,-0.02,-0.01/
    DATA VCD0 /0.0085,0.008,0.0077,0.0078,0.0078,0.0079,0.0114/
    DATA VDF /3.0,7.75,6.2,6.0,5.9,5.5,4.0/
    DATA VCN1 /1.45,1.2,1.05,0.92,0.68,0.5,0.18/
    
    DATA VTP /1.7,1.8,2.0,2.5,3.0,3.3,4.3/
    DATA VTF /3.0,2.5,2.2,2.0,2.0,2.0,2.0/
    DATA VTV /6.0,6.0,6.0,6.0,6.0,6.0,4.0/
    DATA VTAOV1 /7.0,9.0,9.0,9.0,9.0,9.0,9.0/
    CM0=-0.03
    AFA0=0.0D0
!----------------------------------------------------addddjussst begin----------------------------
!!!for sc1095
!    VALCSRF=0.13
!!    AFA0=-0.2711/180*3.14
!!    VAFA1(1)=12.9
!!    VSS1(1)=0.9896
!!    VSS2(1)=3.435
!    CM0=-0.03
!!    VK0(1)=-0.0001
!!    VK1(1)=-0.1277
!!    VK2(1)=0.0306 
!    VCN1=1.4
!    VDF=3
!for sc1095
    !VALCSRF=0.13
!    AFA0=-0.2711/180*3.14
!    VAFA1(1)=12.9
!    VSS1(1)=0.9896
!    VSS2(1)=3.435
    !CM0=-0.03
!    VK0(1)=-0.0001
!    VK1(1)=-0.1277
!    VK2(1)=0.0306 
    !VCN1=1.4
    !VDF=3
!----------------------------------------------------addddjussst end----------------------------
    IF((MA-0.3).LE.1e-13) THEN
        ALCSRF=VALCSRF(1)
        ALCSRF=ALCSRF/PI*180.0
        AFA1=VAFA1(1)
        AFA1=AFA1*PI/180.0
        DAFA1=DVAFA1(1)
        DAFA1=DAFA1*PI/180.0 
        SS1=VSS1(1)
        SS2=VSS2(1)
        K0=VK0(1)
        K1=VK1(1)
        K2=VK2(1)
        CD0=VCD0(1)
        DF=VDF(1)   
        CN1=VCN1(1)
        TP=VTP(1)
        TF=VTF(1)
        TV=VTV(1)
        TAOV1=VTAOV1(1)
    ELSE IF((MA-0.8).GE.1e-13) THEN
        ALCSRF=VALCSRF(7)
        ALCSRF=ALCSRF/PI*180.0
        AFA1=VAFA1(7)
        AFA1=AFA1*PI/180.0
        DAFA1=DVAFA1(7)
        DAFA1=DAFA1*PI/180.0 
        SS1=VSS1(7)
        SS2=VSS2(7)
        K0=VK0(7)
        K1=VK1(7)
        K2=VK2(7)
        CD0=VCD0(7)
        DF=VDF(7)   
        CN1=VCN1(7)
        TP=VTP(7)
        TF=VTF(7)
        TV=VTV(7)
        TAOV1=VTAOV1(7)
    ELSE
         call rxlct2(imk, ss, 7, MA, VMA) 
        IJ=IMK
        J=IMK+1
        ALCSRF=SS*VALCSRF(J)+(1.0-SS)*VALCSRF(IJ)
        ALCSRF=ALCSRF/PI*180.0
        AFA1=SS*VAFA1(J)+(1.0-SS)*VAFA1(IJ)
        AFA1=AFA1*PI/180.0
        DAFA1=SS*DVAFA1(J)+(1.0-SS)*DVAFA1(IJ)
        DAFA1=DAFA1*PI/180.0
        SS1=SS*VSS1(J)+(1.0-SS)*VSS1(IJ)
        SS2=SS*VSS2(J)+(1.0-SS)*VSS2(IJ)
        K0=SS*VK0(J)+(1.0-SS)*VK0(IJ)
        K1=SS*VK1(J)+(1.0-SS)*VK1(IJ)
        K2=SS*VK2(J)+(1.0-SS)*VK2(IJ)
        CD0=SS*VCD0(J)+(1.0-SS)*VCD0(IJ)         
        DF=SS*VDF(J)+(1.0-SS)*VDF(IJ)
        CN1=SS*VCN1(J)+(1.0-SS)*VCN1(IJ)
        TP=SS*VTP(J)+(1.0-SS)*VTP(IJ)
        TF=SS*VTF(J)+(1.0-SS)*VTF(IJ)
        TV=SS*VTV(J)+(1.0-SS)*VTV(IJ)
        TAOV1=SS*VTAOV1(J)+(1.0-SS)*VTAOV1(IJ)  
   END IF
   
!----------------------------------------------------addddjussst begin----------------------------
!!for sc1095
!    ALCSRF=ALCSRF*1.1
!!    AFA0=-0.2711/180*3.14
!!    VAFA1(1)=12.9
!!    VSS1(1)=0.9896
!!    VSS2(1)=3.435
!    CM0=-0.03
!!    VK0(1)=-0.0001
!!    VK1(1)=-0.1277
!!    VK2(1)=0.0306 
!    CN1=1.4
!    DF=3
!----------------------------------------------------addddjussst end----------------------------
   
END SUBROUTINE

subroutine LBINI(THIS,delfa,OMG,nstpt)
IMPLICIT NONE
CLASS(LB) :: THIS
real(rdt),intent(in) :: delfa,omg
integer,intent(in) :: nstpt
    call THIS%LB1%lbini1(delfa,OMG)
    call THIS%LB2%lbini2(nstpt)
    call THIS%LB3%lbini3(nstpt)
end subroutine




!CHO1=0    输出附着流计算结果  CHO1=1  输出分离及动态失速流计算结果
!CHO2=0    单涡脱落                   CHO2=1 多重涡脱落
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> THE PURPOSE:
SUBROUTINE LBMODEL(THIS,CNTAL,CMTAL,CCTAL,AFAI,QI,CCHORDRF,MA,CHO1,CHO2)
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> MODULES USED:
IMPLICIT NONE
CLASS(LB) :: THIS
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> DUMMY ARGUMENTS / INTERFACE: 
real(rdt),INTENT(IN) :: CCHORDRF
real(rdt),INTENT(IN) :: AFAI,QI,MA
real(rdt),INTENT(OUT) :: CNTAL,CMTAL,CCTAL
INTEGER,INTENT(IN) :: CHO1,CHO2
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> DUMMY ARGUMENTS / INTERFACE: 
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->       
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->       READDATA1, READDATA2
!->       EIGANS, HOVER
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:

!INTEGER ERR
real(rdt) :: AFA,Q
real(rdt) :: DELTA_S,DELTA_X
real(rdt) :: BETA,TI,ASPEED,PI
real(rdt) :: UTAL
real(rdt) :: CNF,CNCF,CNIF,CMF,CMCF,CMIF,CCF

real(rdt) :: CNFP,CNCFF,CMCFF
real(rdt) :: FPM,FPN,FPPM,FPPN
real(rdt) :: CV,CMV 
real(rdt) :: AFA00
real(rdt) :: MFEN

real(rdt) :: ALCSRF,AFA1,DAFA1,SS1,SS2,K0,K1,K2,CD0,DF,CN1,TP,TF,TV,TAOV1,AFA0,CM0
real(rdt) :: AFA1M,AFA1N
real(rdt) :: ETA
     

    PI=ACOS(0.0D0)*2
    ASPEED=340
    ETA=0.95

    AFA=AFAI
    UTAL=MA*ASPEED
    Q=QI*CCHORDRF/UTAL

    AFA00=0.0D0
    MFEN=2.0

    TI=CCHORDRF/ASPEED
    DELTA_S=2*UTAL*THIS%LB1%DELTA_T/CCHORDRF
    BETA=SQRT(1.0-MA**2)
    
    CALL THIS%GETPARA(ALCSRF,AFA1,DAFA1,SS1,SS2,K0,K1,K2,CD0,DF,CN1,TP,TF,TV,TAOV1,AFA0,CM0,MA)
    
    AFA=AFA-AFA0
    
    CALL THIS%LB1%LBATTACH(CNCF,CNIF,CMCF,CMIF,CCF,AFA,Q,K0,ALCSRF,MA,TI,ETA,DELTA_S)
    CNF=CNCF+CNIF
    CMF=CMCF+CMIF+CM0
!*************************OUTPUT1****************************
    IF(CHO1.EQ.0) THEN

        THIS%LB1%AFASTD=AFA
        THIS%LB1%QSTD=Q

        CCF=CNF*SIN(AFA)-CCF*COS(AFA)+CD0

        CNTAL=CNF
        CMTAL=CMF
        CCTAL=CCF
        RETURN  
    END IF
    
!*************************PRESSURE CORRECTION****************************

    CALL THIS%LB1%PRCORRECT(CNFP,AFA1M,AFA1N,TF,TV,AFA1,TAOV1,DAFA1,CN1,AFA,Q,CNCF,CNIF,TP,DELTA_S)
    
!*************************T.E. SEPERATE****************************

    CALL LBSEPERATE(FPPM,THIS%LB1%FPMSTD,THIS%LB1%D7M,CNFP,ALCSRF,AFA00,AFA1M,SS1,SS2,DELTA_S,TF)
    CALL LBSEPERATE(FPPN,THIS%LB1%FPNSTD,THIS%LB1%D7N,CNFP,ALCSRF,AFA00,AFA1N,SS1,SS2,DELTA_S,TF)
    CNCFF=((1+SQRT(FPPN))/2)**2*CNCF
    CMCFF=((1+SQRT(FPPM))/2)**2*CNCF*(K0+K1*(1-FPPM)+K2*SIN(FPPM**MFEN*PI))
    CCF=ETA*CNCF*SIN(AFA+Q/2-THIS%LB1%XC-THIS%LB1%YC)*SQRT(FPPN)

!*************************DYNAMIC STALL**************************** 
    IF(CNFP.GT.CN1) THEN
        CALL THIS%LB1%DYSTALL(CMV,CCF,CNFP,CN1,TAOV1,DF,CNCF,FPPN,DELTA_S,AFA,TV,CHO2)
        CNF=CNCFF+CNIF+THIS%LB1%CNV
        CMF=CMCFF+CMIF+CMV+CM0
        CCF=CCF
    ELSE
        THIS%LB1%CNV=0.0D0
        CV=0.0D0
        THIS%LB1%CVSTD=0.0D0
        CNF=CNCFF+CNIF
        CMF=CMCFF+CMIF+CM0
        CCF=CCF
    END IF

!*************************OUTPUT2****************************
    IF(CHO1.EQ.1) THEN
    
        THIS%LB1%AFASTD=AFA
        THIS%LB1%QSTD=Q
        
        THIS%LB1%CNCFSTD=CNCF
        
        CCF=CNF*SIN(AFA)-CCF*COS(AFA)+CD0
        CNF=CNF

        CNTAL=CNF
        CMTAL=CMF
        CCTAL=CCF
        RETURN
    ELSE
        WRITE(*,*) "ERROR"
        STOP
    END IF
!*************************END****************************



996 FORMAT(I2,10F10.5) 
997 FORMAT(2F10.5,2X,A50) 
998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
END SUBROUTINE
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> THE PURPOSE:
SUBROUTINE TESTLBMODEL3(THIS)
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> MODULES USED:
IMPLICIT NONE
CLASS(LB) :: THIS
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> DUMMY ARGUMENTS / INTERFACE: 
INTEGER NSPPR,N,I,J
INTEGER :: ERR
REAL(RDT) :: CCHORDRF
REAL(RDT),ALLOCATABLE :: AFA(:,:),Q(:,:),CNTAL(:,:),CMTAL(:,:),CCTAL(:,:),MA(:,:)
REAL(RDT),ALLOCATABLE :: KFRE(:),AFA0(:),AFAA(:),OMG(:)
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> DUMMY ARGUMENTS / INTERFACE: 
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->       
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->       READDATA1, READDATA2
!->       EIGANS, HOVER
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
INTEGER :: NSPPRT,NSTPT,WIDTH
REAL(RDT) :: V,DELFA,PI
character( len = 10000 ) :: copyfilenam

	NSPPR=360
	N=100
	NSPPRT=NSPPR*N
	NSTPT=4
    ALLOCATE(AFA(NSTPT,NSPPRT),Q(NSTPT,NSPPRT),MA(NSTPT,NSPPRT),&
        CNTAL(NSTPT,NSPPRT),CMTAL(NSTPT,NSPPRT),CCTAL(NSTPT,NSPPRT),&
        KFRE(NSTPT),AFA0(NSTPT),AFAA(NSTPT),OMG(NSTPT),&
        STAT=ERR)
    
	IF(ERR.NE.0) THEN
		WRITE(*,998)
		STOP
	END IF
	
	PI=2*ACOS(0.0)
	
	
	MA(1,:)=0.4
	MA(2,:)=0.5
	MA(3,:)=0.6
	MA(4,:)=0.8

    KFRE=40000.0
    WIDTH=10
    CCHORDRF=0.1
    
	DO J=1,NSTPT
	    OMG(J)=KFRE(J)
	    DO I=1,1
            Q(J,I)=0.0
            AFA(J,I)=0.0
        END DO
        DO I=2,WIDTH
            Q(J,I)=KFRE(J)
            AFA(J,I)=5.0/(WIDTH-1)*(I-1)
        END DO
	    DO I=WIDTH+1,NSPPRT
            Q(J,I)=0.0
            AFA(J,I)=5.0
        END DO
    END DO

    DELFA=5.0/(WIDTH-1)*PI/180
    
    AFA=AFA*PI/180.0
    Q=Q*PI/180.0
    DO J=1,NSTPT
        CALL THIS%LBINI(DELFA,OMG(J),NSTPT)
        DO I=1,NSPPRT
            CALL THIS%LBMODEL(CNTAL(J,I),CMTAL(J,I),CCTAL(J,I),AFA(J,I),Q(J,I),CCHORDRF,MA(J,I),0,0)
        END DO
        CALL THIS%LBDELOCA()
    END DO

    DO I=1,NSTPT
        write( copyfilenam,'(i10000)' ) int(I)
        open ( I*10, file = './LBTEST' // trim(adjustl( copyfilenam )) // '.DAT' )
    END DO

    DO I=WIDTH+1,NSPPRT
        DO J=1,NSTPT
            V=MA(J,J)*340
            WRITE(J*10,997) V*DELFA/OMG(J)*2/CCHORDRF*(I-WIDTH),CNTAL(J,I)/(5.0*PI/180),CMTAL(J,I)/(5.0*PI/180),CCTAL(J,I)/(5.0*PI/180)
        END DO
    END DO
    
    
    DO I=1,NSTPT
        CLOSE(I*10)
    END DO
    
    DEALLOCATE(AFA,Q,MA,CNTAL,CMTAL,CCTAL,KFRE,AFA0,AFAA,&
        STAT=ERR)
    
	IF(ERR.NE.0) THEN
		WRITE(*,999)
		STOP
	END IF
    STOP
997 FORMAT(4F26.6) 
998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
END SUBROUTINE
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> THE PURPOSE:
SUBROUTINE TESTLBMODEL4(THIS)
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> MODULES USED:
IMPLICIT NONE
CLASS(LB) :: THIS
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> DUMMY ARGUMENTS / INTERFACE: 
INTEGER NSPPR,N,I,J,K
INTEGER :: ERR
real(rdt) :: CCHORDRF
real(rdt),ALLOCATABLE :: AFA(:,:),Q(:,:),CNTAL(:,:),CMTAL(:,:),CCTAL(:,:),MA(:,:)
real(rdt),ALLOCATABLE :: KFRE(:),AFA0(:),AFAA(:),OMG(:)
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> DUMMY ARGUMENTS / INTERFACE: 
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->       
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->       READDATA1, READDATA2
!->       EIGANS, HOVER
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
INTEGER :: NSPPRT,NSTPT,XMAX1,XMAX2
real(rdt) :: V,DELFA,PI,DY
CHARACTER(LEN=10000) :: COPYFILENAM
REAL(RDT),EXTERNAL :: FIND_MAX
	NSPPR=360
	N=10
	NSPPRT=NSPPR*N
	NSTPT=2
    ALLOCATE(AFA(NSTPT,NSPPRT),Q(NSTPT,NSPPRT),MA(NSTPT,NSPPRT),&
        CNTAL(NSTPT,NSPPRT),CMTAL(NSTPT,NSPPRT),CCTAL(NSTPT,NSPPRT),&
        KFRE(NSTPT),AFA0(NSTPT),AFAA(NSTPT),OMG(NSTPT),&
        STAT=ERR)
    
	IF(ERR.NE.0) THEN
		WRITE(*,998)
		STOP
	END IF
!	pause
	PI=2*ACOS(0.0)
	DELFA=2*PI/NSPPR
	
	AFA0(1)=0.0
	AFAA(1)=4.0
	MA(1,:)=0.4
	KFRE(1)=0.074
	
	AFA0(2)=0.0
	AFAA(2)=4.0
	MA(2,:)=0.7
	KFRE(2)=0.074
	
    DO I=1,NSTPT
        write( copyfilenam,'(i10000)' ) int(I)
        open ( I*10, file = './LBTEST' // trim(adjustl( copyfilenam )) // '.DAT' )
    END DO

    CCHORDRF=0.4
	
	DO K=1,20
	
        KFRE=0.005*K*8
    	
	    DO J=1,NSTPT
            V=MA(J,J)*340
	        OMG(J)=KFRE(J)*2*V/CCHORDRF

            DO I=1,NSPPRT
                AFA(J,I)=AFA0(J)+AFAA(J)*SIN(DELFA*(I-1))
                Q(J,I)=OMG(J)*AFAA(J)*COS(DELFA*(I-1))
            END DO
        END DO
        
        
        AFA=AFA*PI/180.0
        Q=Q*PI/180.0
        DO J=1,NSTPT
            CALL THIS%LBINI(	DELFA,OMG(J),NSTPT)
            DO I=1,NSPPRT
                CALL THIS%LBMODEL(CNTAL(J,I),CMTAL(J,I),CCTAL(J,I),AFA(J,I),Q(J,I),CCHORDRF,MA(J,I),1,0)
            END DO
            CALL THIS%LBDELOCA()
        END DO
    
!	    DO I=1,NSPPRT
            DO J=1,NSTPT
                DY=FIND_MAX(NSPPR,AFA(J,NSPPR*(N-1)+1:NSPPRT),XMAX1)
                DY=FIND_MAX(NSPPR,CNTAL(J,NSPPR*(N-1)+1:NSPPRT),XMAX2)
                WRITE(J*10,997) KFRE(J),DY/(AFAA(J)*PI/180.0),(XMAX1-XMAX2)*DELFA/PI*180.0
!                WRITE(J*10,997) KFRE(J),AFA(J,I)/PI*180,CNTAL(J,I),CMTAL(J,I),CCTAL(J,I)
            END DO
!        END DO
        
    END DO    
    
    DO I=1,NSTPT
        CLOSE(I*10)
    END DO
    
    DEALLOCATE(AFA,Q,MA,CNTAL,CMTAL,CCTAL,KFRE,AFA0,AFAA,&
        STAT=ERR)

	IF(ERR.NE.0) THEN
		WRITE(*,999)
		STOP
	END IF
997 FORMAT(4F26.6) 
998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
END SUBROUTINE
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
END MODULE
