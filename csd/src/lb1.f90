MODULE LB1_CLASS
USE GlobalDataFun
IMPLICIT NONE

TYPE,PUBLIC :: LB1

    real(rdt) :: delta_t
    
    real(rdt) :: afastd,qstd
    real(rdt) :: xc,yc,zc
    real(rdt) :: d1,d2,d3,d4,d5
    real(rdt) :: delta_xx1std,delta_xx2std,delta_xx3std,delta_xx4std 
    real(rdt) :: cncfstd,fpmstd,fpnstd
    real(rdt) :: d6,d7n,d7m
    real(rdt) :: cnv,taov
    real(rdt) :: cvstd

    CONTAINS

    PROCEDURE,PUBLIC :: LBINI1
    PROCEDURE,PUBLIC :: PRELAG
    PROCEDURE,PUBLIC :: LBCORR
    PROCEDURE,PUBLIC :: DYSTALL
    PROCEDURE,PUBLIC :: PRCORRECT
    PROCEDURE,PUBLIC :: LBATTACH

end TYPE LB1
    CONTAINS

    SUBROUTINE LBATTACH(THIS,CNCF,CNIF,CMCF,CMIF,CCF,AFA,Q,K0,ALCSRF,MA,TI,ETA,DELTA_S)
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> MODULES USED:
IMPLICIT NONE
CLASS(LB1) :: THIS
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> DUMMY ARGUMENTS / INTERFACE: 
real(rdt),INTENT(IN) :: AFA,Q,ALCSRF,MA,TI,DELTA_S,K0,ETA
real(rdt),INTENT(OUT) :: CNCF,CNIF,CMCF,CMIF,CCF
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
real(rdt) :: A1,A2,A3,A4,A5,B1,B2,B3,B4,B5
real(rdt) :: BETA,CM0,PI
real(rdt) :: ED,EQ,KI,KQ
real(rdt) :: DELTA_X,DELTA_XX1,DELTA_XX2,DELTA_XX3,DELTA_XX4
real(rdt) :: CNQF,CMQF,ONE=1.0

	A1=0.3
	A2=0.7
	B1=0.14
	B2=0.53
!   A1=0.636
!   A2=0.364
!   B1=0.339
!   B2=0.249  
!   A1=0.482
!   A2=0.518
! 	B1=0.684
!   B2=0.235  
	A3=1.5
	B3=0.25
	A4=-0.5
	B4=0.1
	A5=1.0
	B5=0.5
    
	CM0=0.0D0 

	BETA=SQRT(1.0-MA**2)
	PI=ACOS(0.0D0)*2
    
!*************************CNCF****************************

	DELTA_X=(AFA+Q/2)-(THIS%AFASTD+THIS%QSTD/2)
	CALL LBLAG(THIS%XC,EXP(-B1*BETA**2*DELTA_S),A1,DELTA_X)
	CALL LBLAG(THIS%YC,EXP(-B2*BETA**2*DELTA_S),A2,DELTA_X)
	CNCF=ALCSRF*(AFA+Q/2-THIS%XC-THIS%YC)
	CCF=ETA*CNCF*TAN(AFA+Q/2-THIS%XC-THIS%YC)


!*************************CNIF****************************

   
!   KI=0.75/(1.0-MA+PI*BETA*MA**2*(A1*B1+A2*B2))
    KI=0.75/(1.0-MA+ALCSRF/2*BETA**2*MA**2*(A1*B1+A2*B2))
	ED=EXP(-THIS%DELTA_T/KI/TI)
	DELTA_XX1=(AFA-THIS%AFASTD)/THIS%DELTA_T

	CALL LBLAG(THIS%D1,ED,ONE,DELTA_XX1-THIS%DELTA_XX1STD)

	CNIF=4.*(KI*TI)/MA*(DELTA_XX1-THIS%D1)
    

!*************************CNQF****************************

!   KQ=0.75/(1.0-MA+2*PI*BETA*MA**2*(A1*B1+A2*B2)) 
	KQ=0.75/(1.0-MA+ALCSRF*BETA**2*MA**2*(A1*B1+A2*B2)) 
	EQ=EXP(-THIS%DELTA_T/KQ/TI)
	DELTA_XX2=(Q-THIS%QSTD)/THIS%DELTA_T

	CALL LBLAG(THIS%D2,EQ,ONE,DELTA_XX2-THIS%DELTA_XX2STD)  

	CNQF=(KQ*TI)/MA*(DELTA_XX2-THIS%D2)
    
!*************************CMCF****************************

	DELTA_X=Q-THIS%QSTD
	CALL LBLAG(THIS%ZC,EXP(-B5*BETA**2*DELTA_S),A5,DELTA_X)  
	CMCF=-ALCSRF/16.0*(Q-THIS%ZC)

	CMCF=CM0+K0*CNCF+CMCF
    
!*************************CMIF****************************

    KI=0.8*(A3*B4+A4*B3)/(B3*B4*(1-MA)) 
!   KI=(A3*B4+A4*B3)/(B3*B4*(1-MA)) 
	ED=EXP(-THIS%DELTA_T/KI/TI/B3)
	EQ=EXP(-THIS%DELTA_T/KI/TI/B4)
	DELTA_XX3=(AFA-THIS%AFASTD)/THIS%DELTA_T

	CALL LBLAG(THIS%D3,ED,ONE,DELTA_XX3-THIS%DELTA_XX3STD)  
	CALL LBLAG(THIS%D4,EQ,ONE,DELTA_XX3-THIS%DELTA_XX3STD) 

	CMIF=-(KI*TI)/MA*(A3*B3*(DELTA_XX3-THIS%D3)+A4*B4*(DELTA_XX3-THIS%D4))
        
!*************************CMQF****************************

    KQ=5.6/(15*(1-MA)+3*ALCSRF/2*BETA**2*MA**2*A5*B5)
!   KQ=7.0/(15*(1-MA)+3*PI*BETA*MA**2*A5*B5)
	EQ=EXP(-THIS%DELTA_T/KQ/TI)
	DELTA_XX4=(Q-THIS%QSTD)/THIS%DELTA_T

	CALL LBLAG(THIS%D5,EQ,ONE,DELTA_XX4-THIS%DELTA_XX4STD) 

	CMQF=-7.0/12.0*(KQ*TI)/MA*(DELTA_XX4-THIS%D5)


	CNCF=CNCF
	CNIF=CNIF+CNQF

	CMCF=CMCF
	CMIF=CMIF+CMQF

	CCF=CCF

	THIS%DELTA_XX1STD=DELTA_XX1
	THIS%DELTA_XX2STD=DELTA_XX2
	THIS%DELTA_XX3STD=DELTA_XX3 
	THIS%DELTA_XX4STD=DELTA_XX4
    
    
END SUBROUTINE




SUBROUTINE PRCORRECT(THIS,CNFP,AFA1M,AFA1N,TF,TV,AFA1,TAOV1,DAFA1,CN1,AFA,Q,CNCF,CNIF,TP,DELTA_S)
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> MODULES USED:
IMPLICIT NONE
CLASS(LB1) :: THIS
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> DUMMY ARGUMENTS / INTERFACE: 
real(rdt),INTENT(IN) :: AFA,Q,CNCF,CNIF,TP,DELTA_S,TAOV1,DAFA1,CN1,AFA1
real(rdt),INTENT(OUT) :: CNFP,AFA1M,AFA1N
real(rdt),INTENT(INOUT) :: TF,TV
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

    CALL THIS%PRELAG(CNFP,CNCF,CNIF,TP,DELTA_S)
    
    IF(CNFP.GT.CN1) THEN
        THIS%TAOV=THIS%TAOV+DELTA_S
    ELSE
        THIS%TAOV=-DELTA_S
    END IF
    
    CALL THIS%LBCORR(CNFP,AFA1M,AFA1N,TF,TV,AFA1,TAOV1,DAFA1,CN1,AFA,Q,DELTA_S)
    
END SUBROUTINE





SUBROUTINE PRELAG(THIS,CNFP,CNCF,CNIF,TP,DELTA_S)
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> MODULES USED:
IMPLICIT NONE
CLASS(LB1) :: THIS
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> DUMMY ARGUMENTS / INTERFACE: 
real(rdt),INTENT(IN) :: CNCF,CNIF,TP,DELTA_S
real(rdt),INTENT(OUT) :: CNFP
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
real(rdt) :: EP
real(rdt) :: DELTA_X,ONE=1.0

    EP=EXP(-DELTA_S/TP)
    DELTA_X=CNCF-THIS%CNCFSTD

    CALL LBLAG(THIS%D6,EP,ONE,DELTA_X) 
    
    CNFP=CNCF+CNIF-THIS%D6
    
END SUBROUTINE
    


SUBROUTINE LBCORR(THIS,CNFP,AFA1M,AFA1N,TF,TV,AFA1,TAOV1,DAFA1,CN1,AFA,Q,DELTA_S)
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> MODULES USED:
IMPLICIT NONE
CLASS(LB1) :: THIS
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> DUMMY ARGUMENTS / INTERFACE: 
real(rdt),INTENT(IN) :: AFA,Q,TAOV1,DAFA1,CN1,DELTA_S,AFA1
real(rdt),INTENT(IN) :: CNFP
real(rdt),INTENT(INOUT) :: TF,TV,AFA1N,AFA1M
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



!+++++++++++++++++++++++CORRECTION1++++++++++++++++++++++++++++++

	AFA1M=AFA1
	AFA1N=AFA1
	IF((CNFP.GT.CN1).AND.(AFA.LT.THIS%AFASTD)) THEN
		AFA1N=AFA1N-(1-THIS%FPNSTD)**0.25*DAFA1
    END IF


!+++++++++++++++++++++++CORRECTION2++++++++++++++++++++++++++++++

    IF((THIS%TAOV.GT.0).OR.(THIS%FPNSTD.LE.0.7)) THEN 
!    IF(CNFP.LT.CN1)THEN 
        TF=TF/2.0
	END IF   

!+++++++++++++++++++++++CORRECTION3++++++++++++++++++++++++++++++

    IF((THIS%TAOV.GT.0).AND.(AFA.GT.THIS%AFASTD)) THEN 
        TF=TF/2.0
! forsc1095 begin 
	    TV=TV*8.0
! forsc1095 end 
	END IF   

!+++++++++++++++++++++++CORRECTION4++++++++++++++++++++++++++++++

    	IF((CNFP.GT.CN1).AND.(AFA.LT.THIS%AFASTD)) THEN
!    IF(AFA.LT.AFASTD) THEN
	    TV=TV/2.0
! forsc1095 begin 
        TF=TF*4.0
! forsc1095 end 
    	END IF

!!	
!	IF((TAOV.GT.0).AND.(TAOV.LT.TAOV1).AND.(QSTD*Q.LT.0)) THEN
!	    TF=TF/2.0
!	    TV=TV/2.0
!	END IF
!!	
!	IF((TAOV.GT.TAOV1)) THEN
!	    TF=TF*4.0
!	END IF
!	
!


!    IF((CNFP.LT.CN1).AND.(AFA.LT.AFASTD)) THEN

    
!    IF((CNFP.GT.CN1).AND.(AFA.GT.AFASTD)) THEN
!!    IF(AFA.LT.AFASTD) THEN
!	    AFA1=AFA1-FPSTD**0.25*DAFA1
!    END IF
    
END SUBROUTINE


SUBROUTINE DYSTALL(THIS,CMV,CCF,CNFP,CN1,TAOV1,DF,CNCF,FPP,DELTA_S,AFA,TV,CHO)
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> MODULES USED:
IMPLICIT NONE
CLASS(LB1) :: THIS
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> DUMMY ARGUMENTS / INTERFACE: 
real(rdt),INTENT(IN) :: CNFP,CN1,TAOV1,DF,CNCF,FPP,DELTA_S,AFA,TV
INTEGER,INTENT(IN) :: CHO
real(rdt),INTENT(OUT) :: CMV
real(rdt),INTENT(INOUT) :: CCF
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
real(rdt) :: EV,KV,CV,DELTA_X,CPV,BIGFA,PI,TST

    PI=ACOS(0.0D0)*2
   
    IF((THIS%TAOV.GT.0).AND.(THIS%TAOV.LT.TAOV1)) THEN
        CALL VLBU(THIS%CNV,CV,THIS%CVSTD,CNCF,FPP,DELTA_S,TV,1)
        CPV=0.2*(1-COS(PI*THIS%TAOV/TAOV1))
        CMV=-CPV*THIS%CNV
    ELSE IF ((THIS%TAOV.GT.TAOV1).AND.(THIS%TAOV.LT.2*TAOV1)) THEN
        CALL VLBU(THIS%CNV,CV,THIS%CVSTD,CNCF,FPP,DELTA_S,TV,0)
        CPV=0.2*(1-COS(PI*THIS%TAOV/TAOV1))
        CMV=-CPV*THIS%CNV
    ELSE 
        IF((CHO.EQ.1).AND.(AFA.GT.THIS%AFASTD).AND.(THIS%TAOV.GT.TAOV1)) THEN
            TST=2*(1-FPP)/0.20
            THIS%TAOV=-TST+TAOV1
        END IF
        CALL VLBU(THIS%CNV,CV,THIS%CVSTD,CNCF,FPP,DELTA_S,TV,0)
        CMV=0.0
    END IF
        
    
!    IF((TAOV.GT.0).AND.(TAOV.LT.2*TAOV1)) THEN
    BIGFA=FPP**(DF*(CNFP-CN1))
    CCF=CCF*BIGFA

    
    THIS%CVSTD=CV
    
END SUBROUTINE



SUBROUTINE VLBU(CNVT,CV,CVSTD,CNCF,FPP,DELTA_S,TV,CHO)
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> MODULES USED:
IMPLICIT NONE
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> DUMMY ARGUMENTS / INTERFACE: 
real(rdt),INTENT(IN) :: CNCF,FPP,DELTA_S,TV,CVSTD
INTEGER,INTENT(IN) :: CHO
real(rdt),INTENT(OUT) :: CNVT,CV
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
real(rdt) :: EV,KV,DELTA_X,PI,ONE=1.0

    PI=ACOS(0.0D0)*2
    
    EV=EXP(-DELTA_S/TV)

    IF(CHO.EQ.1) THEN
        KV=(1.0+SQRT(FPP))**2.0/4.0
        CV=CNCF*(1-KV)
        DELTA_X=CV-CVSTD
        
        CALL LBLAG(CNVT,EV,ONE,DELTA_X)
    ELSE IF(CHO.EQ.0) THEN
        KV=(1.0+SQRT(FPP))**2.0/4.0
        CV=CNCF*(1-KV)
        
        CNVT=CNVT*EV
    ELSE
        WRITE(*,*) 'VLBU ERROR'
    END IF
 
END SUBROUTINE



SUBROUTINE LBLAG(DI,EX,COE,DELTA)
IMPLICIT NONE
real(rdt),INTENT(INOUT) :: DI
real(rdt),INTENT(IN) :: EX,COE,DELTA

	    DI=DI*EX+COE*DELTA*EX**0.5
	    
END SUBROUTINE

SUBROUTINE LBSEPERATE(FPP,FPSTD,D7,CNFP,ALCSRF,AFA00,AFA1,SS1,SS2,DELTA_S,TF)
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> MODULES USED:
IMPLICIT NONE
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> DUMMY ARGUMENTS / INTERFACE: 
real(rdt),INTENT(IN) :: CNFP,ALCSRF,AFA00,AFA1,SS1,SS2,DELTA_S,TF
real(rdt),INTENT(INOUT) :: FPSTD,D7
real(rdt),INTENT(OUT) :: FPP
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
real(rdt) :: AFAEF,PI,FP,DELTA_X,EF,ONE=1.0

    
	AFAEF=CNFP/ALCSRF
    PI=ACOS(0.0D0)*2
    
    IF(AFAEF.LE.AFA1) THEN
        FP=1-0.3*EXP((AFAEF-AFA00-AFA1)/PI*180/SS1)
    ELSE
        FP=0.04+0.66*EXP((AFA1-(AFAEF-AFA00))/PI*180/SS2)
    END IF 

    
    EF=EXP(-DELTA_S/TF)
    DELTA_X=FP-FPSTD
    
    CALL LBLAG(D7,EF,ONE,DELTA_X) 
	
    FPP=FP-D7

    FPSTD=FP
  
    
END SUBROUTINE

SUBROUTINE LBINI1(THIS,DELFA,OMG)
IMPLICIT NONE
CLASS(LB1) :: THIS
REAL(RDT),INTENT(IN) :: DELFA,OMG
    THIS%DELTA_T=DELFA/OMG
	
    THIS%AFASTD=0.0D0
    THIS%QSTD=0.0D0
    THIS%CNCFSTD=0.0D0
    THIS%FPMSTD=0.0D0
    THIS%FPNSTD=0.0D0
    THIS%CVSTD=0.0D0

    
    THIS%XC=0.0D0
    THIS%YC=0.0D0
    THIS%ZC=0.0D0
    
    THIS%D1=0.0D0
    THIS%D2=0.0D0
    THIS%D3=0.0D0
    THIS%D4=0.0D0
    THIS%D5=0.0D0
    THIS%D6=0.0D0
    THIS%D7M=0.0D0
    THIS%D7N=0.0D0
       
    THIS%DELTA_XX1STD=0.0D0
    THIS%DELTA_XX2STD=0.0D0
    THIS%DELTA_XX3STD=0.0D0
    THIS%DELTA_XX4STD=0.0D0
    
    THIS%CNV=0.0D0
    THIS%TAOV=0.0D0
    

END SUBROUTINE
END MODULE