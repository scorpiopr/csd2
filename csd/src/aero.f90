MODULE AERO_CLASS
USE GlobalDataFun
USE INFLOW_CLASS
USE SECTIONF_CLASS
USE BEAM_CLASS
USE MBEAMAERO
USE VORTEX_CLASS
USE LB_CLASS
IMPLICIT NONE
!->+++++++++++++++++++++++++++++++++++++++++++++
TYPE,PUBLIC :: AERO

    REAL(RDT),ALLOCATABLE :: FLCSTCOE(:,:),MLCSTCOE(:,:),AOALCST(:)

    REAL(RDT),ALLOCATABLE :: BVQA(:,:,:,:)!共轴诱导速度
    REAL(RDT),ALLOCATABLE ::PSIV(:)
    
    real(RDT), allocatable:: CFDFNR(:, :),CFDMNR(:,:)
    real(RDT), allocatable:: CFDFPSI(:, :, :),CFDMPSI(:,:, :)
    
    INTEGER :: NSTP

    INTEGER :: MKIM

    REAL(RDT) :: ALC,CD0,XAD,BSC,cchordrf

    REAL(RDT) :: MU,CTC,CT
    REAL(RDT) :: LAMIL,LAMI,LAMD
    REAL(RDT) :: alp,thetapa,delta
    INTEGER :: nrbd,MKFM
    REAL(RDT) :: RHA
    REAL(RDT) :: ALCSRF,sgslid

    REAL(RDT) :: trf,mrf

    TYPE(INFLOW) :: INF
    TYPE(SECTIONF) :: SEF
    TYPE(LB) :: LB

    CONTAINS
    PROCEDURE,PUBLIC :: CONSTRUCT_AERO
    PROCEDURE,PUBLIC :: CONSTRUCT_AERO2
    PROCEDURE,PUBLIC :: MLEGGAUSS_AERO
    PROCEDURE,PUBLIC :: GETTH0
    PROCEDURE,PUBLIC :: GETTH02
    PROCEDURE,PUBLIC :: GETCT
    PROCEDURE,PUBLIC :: UNIFLOW
    PROCEDURE,PUBLIC :: LINFLOW
    PROCEDURE,PUBLIC :: FMEL
    PROCEDURE,PUBLIC :: LADIDHV
    PROCEDURE,PUBLIC :: beamaero

    PROCEDURE,PRIVATE :: AEROMODEL
    PROCEDURE,PRIVATE :: GETLINDV
    PROCEDURE,PRIVATE :: GREENBERGST
    PROCEDURE,PRIVATE :: GREENBERG
    PROCEDURE,PRIVATE :: AEROST
    PROCEDURE,PUBLIC :: TECPLOT_BVQA
	
END TYPE AERO
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++PUBLIC++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
CONTAINS
!->+++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE CONSTRUCT_AERO(THIS,INPU,OMG,R)
IMPLICIT NONE
CLASS(AERO) :: THIS
CLASS(INPUT) :: INPU
REAL(RDT),INTENT(IN) :: OMG,R
INTEGER :: I

    CALL THIS%beamaero(INPU)

    ALLOCATE(THIS%FLCSTCOE(6,NPTS),&
                        THIS%MLCSTCOE(6,NPTS),&
                        THIS%AOALCST(NPTS))
    ALLOCATE(THIS%CFDFNR(3,NPTS),THIS%CFDMNR(3,NPTS))

    THIS%ALC=INPU%IPT%ALCSRF
    THIS%CD0=0.01
    THIS%XAD=0.0D0
    THIS%cchordrf=INPU%IPT%cchordrf
    THIS%BSC=THIS%cchordrf/2
    THIS%thetapa=INPU%IPT%thetapa
    THIS%delta=INPU%IPT%delta
    THIS%alp=THIS%thetapa+THIS%delta

    THIS%MU=INPU%IPT%MU
    THIS%CTC=INPU%IPT%CTC
    THIS%CT=THIS%CTC/2

    THIS%nrbd=INPU%IPT%nrbd
    THIS%MKFM=INPU%IPT%MKFM
    THIS%RHA=INPU%IPT%RHA
    THIS%ALCSRF=INPU%IPT%ALCSRF

    THIS%MKIM=INPU%IPT%MKIM

    CALL THIS%INF%CONSTRUCT_INFLOW(INPU)
    CALL THIS%SEF%CONSTRUCT_SECTIONF(INPU)

    THIS%sgslid=THIS%nrbd*THIS%cchordrf/R/PI

    THIS%trf=THIS%RHA*PI*OMG**2*R**4
    THIS%MRF=THIS%trf*R

END SUBROUTINE
!->+++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE CONSTRUCT_AERO2(THIS,NSTP,OMG,INPU)
USE MFWM
IMPLICIT NONE
CLASS(AERO) :: THIS
TYPE(INPUT) :: INPU
INTEGER,INTENT(IN) :: NSTP
REAL(RDT),INTENT(IN) :: OMG
INTEGER :: I

    THIS%NSTP=NSTP
    ALLOCATE(THIS%BVQA(3,NPTS,nstp,2),THIS%PSIV(NSTP))
    ALLOCATE(THIS%CFDFPSI(3,NPTS, nstp),THIS%CFDMPSI(3,NPTS, NSTP))

    DO I=1,NSTP
        THIS%PSIV(I)=2.0*PI/(NSTP-1)*(I-1)
    END DO

    IF(THIS%MKIM.EQ.3.OR.THIS%MKIM.EQ.4) THEN
        CALL VTX%CONSTRUCT_VOTEX(INPU,THIS%NSTP,NPTS)
    END IF

    SELECT CASE(THIS%MKFM)
        case(3)
            call THIS%LB%LBini(360.0/(NSTP-1)*PI/180,OMG,NPTS)
!        CASE(4)
!            call cfdini()
    end select

END SUBROUTINE
!->+++++++++++++++++++++++++++++++++++++++++++++
subroutine beamaero(THIS,INPU)    
USE MBEAMAERO
implicit none
CLASS(AERO) :: THIS
TYPE(INPUT) :: INPU

INTEGER :: I
REAL(RDT) :: C
    IF(.NOT.ALLOCATED(RLUDAUD)) THEN
        NPTS=INPU%IPT%NPTS
        ALLOCATE(RLUDAUD(NPTS))
        ALLOCATE(FLCST_ST(3,NPTS),MLCST_ST(3,NPTS))
        ALLOCATE(FLCST_FR(3,NPTS),MLCST_FR(3,NPTS))

        NNOD=INPU%MAT%NBPL
        ALLOCATE(rlud(NNOD))

        DO I=1,NNOD!按展向截面循环
            RLUD(I)=INPU%MAT%RL(I)
        END DO

        DO I=1,NPTS
            
            RLUDAUD(I)=(RLUD(NNOD)-RLUD(1))/(NPTS-1)*(I-1)+RLUD(1)!确定展向气动点相对位置：均匀分布
        END DO
		
		!C=0.9
        !DO I=1,NPTS
        !    RLUDAUD(I)=(RLUD(NNOD)-RLUD(1))*sin(C*PI/2*(i-1)/(NPTS-1))/sin(C*PI/2)+RLUD(1)
        !END DO


    END IF

end subroutine beamaero 
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE GETTH0(THIS,TH0)
IMPLICIT NONE
CLASS(AERO) :: THIS
REAL(RDT),INTENT(OUT) :: TH0
REAL(RDT) :: TH75

    TH75=3.0*THIS%CTC/THIS%alcsrf/THIS%sgslid+1.5*THIS%LAMI
    TH0=0.8*TH75 

END SUBROUTINE
!->+++++++++++++++++++++++++++++++++++++++++++++
REAL(RDT) FUNCTION GETTH02(THIS)
IMPLICIT NONE
CLASS(AERO) :: THIS

    GETTH02=6.0*THIS%CT/THIS%alcsrf/THIS%sgslid+1.5*THIS%LAMI

END FUNCTION
!->+++++++++++++++++++++++++++++++++++++++++++++
REAL(RDT) FUNCTION GETCT(THIS,TH75)
IMPLICIT NONE
CLASS(AERO) :: THIS
REAL(RDT),INTENT(IN) :: TH75
    GETCT=THIS%ALCSRF*THIS%SGSLID/2.0*(TH75/3.0-THIS%LAMI/2.0)!SGSLID:实度。公式取自Principles of Helicopter Aerodynamics Leishman中式4.111,其中前进比为0，无线性扭转

END FUNCTION
!->+++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE UNIFLOW(THIS)
IMPLICIT NONE
CLASS(AERO) :: THIS

    CALL THIS%INF%VEQ3(THIS%ALP,THIS%CTC,THIS%MU,THIS%LAMI,THIS%LAMD)
    THIS%LAMIL=THIS%LAMI

END SUBROUTINE
!->+++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE LADIDHV(THIS,TH75)
implicit none
CLASS(AERO) :: THIS
REAL(RDT), INTENT(IN):: TH75
REAL(RDT) XX

    XX=THIS%alcsrf*THIS%sgslid
    
    THIS%LAMI=(SQRT(1.0+24.0*TH75/XX)-1.0)*XX/16.0
    THIS%LAMIL=THIS%LAMI

END SUBROUTINE
!->+++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE LINFLOW(THIS,rb,PSI,J)
IMPLICIT NONE
CLASS(AERO) :: THIS
REAL(RDT),INTENT(IN) :: RB,PSI
INTEGER :: J
 
    CALL THIS%INF%linearinflow(THIS%lamil, THIS%lami, THIS%lamd, THIS%mu, rb, psi,J)

END SUBROUTINE
!->+++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE MLEGGAUSS_AERO(THIS,OMG,BLADE,TH0,TH1C,TH1S,PSI,SOL)
IMPLICIT NONE
CLASS(AERO) :: THIS
INTEGER,INTENT(IN) :: SOL
REAL(RDT),INTENT(IN)::TH0,TH1C,TH1S,PSI,OMG
CLASS(BEAM) :: BLADE
INTEGER :: NBPLIA,FCHO

REAL(RDT),EXTERNAL:: FABSC,triseries

    CALL BLADE%GETTHP(TH0, TH1C, TH1S, PSI, OMG)
    
    FCHO=0
    IF(SOL.EQ.1.or.SOL.EQ.2) THEN
        DO nbplia=1,npts
            CALL THIS%FMEL(NBPLIa,OMG,BLADE,PSI,FCHO,SOL)
        END DO
    END IF

END SUBROUTINE
!->+++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE FMEL(THIS,NBPLIa,OMG,BLADE,PSI,FCHO,SOL)
USE MBEAMAERO
IMPLICIT NONE
CLASS(AERO) :: THIS
CLASS(BEAM) :: BLADE
REAL(RDT),INTENT(IN) :: OMG,PSI
INTEGER,INTENT(IN)::NBPLIa,FCHO,SOL

REAL(RDT) SS,LL,RXUD,he

INTEGER NBPLI,NBPLIP1 

real(rdt) rb

REAL(RDT) TM1(3,3)
REAL(RDT) TER(3,3)
REAL(RDT) TDE(3,3),TDET(3,3),TDED1T(3,3)

REAL(RDT) REAE(3),VEAE(3),VEAED1T(3),VAE(3),VAED1T(3),UDP(3),UDPD1T(3)
REAL(RDT) UF(3),TRN(3,3),TBRBP(3,3),TBRTP(3,3),TRND1T(3,3),TBRTPD1T(3,3)

REAL(RDT) TVCP1(3),TVCP2(3),TVCP3(3),TVCP4(3)
REAL(RDT) TMV1(3),TMV2(3)

REAL(RDT) TAF,TAFD1T,TAFD2T,TAFG,TAFGD1T,TAFGD2T 

REAL(RDT) TX1

REAL(RDT) INDVL(3)

INTEGER UDPCHO

REAL(RDT) LDST,LDAT
 
    rxud=RLUDAUD(NBPLIa)
    call rxlct(nbpli, ss, ll, NNOD, rxud, rlud) 

    ll=ll*BLADE%R
    HE=BLADE%ELEMEN(NBPLI)%NOD(1)%RLUD*BLADE%R-BLADE%EL1

    NBPLIP1=NBPLI+1
    CALL BLADE%ELEMEN(NBPLI)%GETLV(HE,SS,LL,BLADE%BTAP,BLADE%EL1,BLADE%THP,BLADE%THP1D,BLADE%THP2D,OMG)

    UDPCHO=1

    RB=BLADE%ELEMEN(NBPLI)%XX+BLADE%EL1+HE
    RB=RB/BLADE%R

    CALL THIS%GETLINDV(INDVL,PSI,RB,BLADE%R,OMG,SOL)


    IF(BLADE%ELEMEN(NBPLI)%IS_TIP.EQ.1) THEN
        LDST=-BLADE%ELEMEN(NBPLI)%LDT(3)
        LDAT=-BLADE%ELEMEN(NBPLI)%LDT(2)
    ELSE 
        LDST=0.0D0
        LDAT=0.0D0
    END IF

    TX1=COS(LDST)*COS(LDAT)
    TAFG=(BLADE%THP+BLADE%BTAJ)*TX1+(BLADE%ELEMEN(NBPLI)%BTAT-BLADE%BTAJ)
    TAFGD1T=BLADE%THP1D*TX1
    TAFGD2T=BLADE%THP2D*TX1

    TAF=TAFG+BLADE%ELEMEN(NBPLI)%SPH
    TAFD1T=TAFGD1T+BLADE%ELEMEN(NBPLI)%SPH1T
    TAFD2T=TAFGD2T+BLADE%ELEMEN(NBPLI)%SPH2T

    GET_UDP:SELECT CASE(UDPCHO) 
    CASE(0)

!               UF=0.0D0
!           UF=INDVL
!           UF(1)=MU*OMG*R+UF(1)
!           UF(3)=-MU*TAN(ALP)*OMG*R+UF(3)
!
!
!! chopra 1981   
!        UDP(2)=OMG*(ss*ll+el1+he)*((1-SV1X**2/2.0)*COS(TAF)-SV1X*SW1X*SIN(TAF))+ &
!               COS(TAF)*(SV1T-SPH*(-1.0/2*BSCT)*SIN(TAF)+ &
!               OMG*(SU-SV1X*(-1.0/2*BSCT)*COS(TAF)-SW1X*(-1.0/2*BSCT)*SIN(TAF))-OMG*BTAP*(SW+(-1.0/2*BSCT)*SIN(TAF)))+ &
!               SIN(TAF)*(SW1T+SPH*(-1.0/2*BSCT)*COS(TAF)+OMG*BTAP*(SV+(-1.0/2*BSCT)*COS(TAF))+LAMIL*OMG*R)+ &
!               (SV1X*COS(TAF)+SW1X*SIN(TAF))*OMG*(SV+(-1.0/2*BSCT)*COS(TAF))
!               
!        UDP(3)=-OMG*(ss*ll+el1+he)*((1-SV1X**2/2.0)*SIN(TAF)-SV1X*SW1X*COS(TAF))- &
!               SIN(TAF)*(SV1T-SPH*(-1.0/2*BSCT)*SIN(TAF)+ &
!               OMG*(SU-SV1X*(-1.0/2*BSCT)*COS(TAF)-SW1X*(-1.0/2*BSCT)*SIN(TAF))-OMG*BTAP*(SW+(-1.0/2*BSCT)*SIN(TAF)))+ &
!               COS(TAF)*(SW1T+SPH*(-1.0/2*BSCT)*COS(TAF)+OMG*BTAP*(SV+(-1.0/2*BSCT)*COS(TAF))+LAMIL*OMG*R)- &
!               (SV1X*SIN(TAF)-SW1X*COS(TAF))*OMG*(SV+(-1.0/2*BSCT)*COS(TAF))      
    CASE(1)
        
        REAE(1)=BLADE%ELEMEN(NBPLI)%XX+BLADE%ELEMEN(NBPLI)%SU
        REAE(2)=BLADE%ELEMEN(NBPLI)%SV
        REAE(3)=BLADE%ELEMEN(NBPLI)%SW
        REAE=REAE+BLADE%ELEMEN(NBPLI)%HEE
         
        CALL VCP(TVCP1,BLADE%ELEMEN(NBPLI)%OMGE,REAE)

        VEAE(1)=BLADE%ELEMEN(NBPLI)%SU1T
        VEAE(2)=BLADE%ELEMEN(NBPLI)%SV1T
        VEAE(3)=BLADE%ELEMEN(NBPLI)%SW1T

        CALL VCP(TVCP2,BLADE%ELEMEN(NBPLI)%OMGE,VEAE)
!        TVCP2=2*TVCP2
        CALL VCP(TVCP3,BLADE%ELEMEN(NBPLI)%OMGE1D,REAE)
        
        VEAE=VEAE+TVCP1+BLADE%ELEMEN(NBPLI)%VBE

        VEAED1T(1)=BLADE%ELEMEN(NBPLI)%SU2T
        VEAED1T(2)=BLADE%ELEMEN(NBPLI)%SV2T
        VEAED1T(3)=BLADE%ELEMEN(NBPLI)%SW2T
        
!       CALL VCP(TVCP4,OMGE,TVCP1)
        
        VEAED1T=VEAED1T+TVCP2+TVCP3+BLADE%ELEMEN(NBPLI)%VBE1D
!       VEAED1T=VEAED1T+TVCP2+TVCP3+TVCP4+VBE1D

        UF=0.0D0
        UF=INDVL
        UF(1)=THIS%MU*OMG*BLADE%R+UF(1)
        UF(3)=-THIS%MU*TAN(THIS%ALP)*OMG*BLADE%R+UF(3)

        CALL TXYZ(TRN,3,PSI)
        CALL TXYZ(TBRBP,2,-BLADE%BTAP)
        CALL TXYZ(TBRTP,1,BLADE%THP)

        CALL MVMU2(TMV1,3,3,TRN,UF)
        CALL MVMU2(TMV1,3,3,TBRBP,TMV1)
        CALL MVMU2(TMV1,3,3,TBRTP,TMV1)
        CALL MVMU2(VAE,3,3,BLADE%ELEMEN(NBPLI)%TEB,TMV1)

        CALL TXYZD1T(TRND1T,3,PSI,OMG)
        CALL TXYZD1T(TBRTPD1T,1,BLADE%THP,BLADE%THP1D)

        CALL MVMU2(TMV1,3,3,TRND1T,UF)
        CALL MVMU2(TMV1,3,3,TBRBP,TMV1)
        CALL MVMU2(TMV1,3,3,TBRTP,TMV1)

        CALL MVMU2(TMV2,3,3,TRN,UF)
        CALL MVMU2(TMV2,3,3,TBRBP,TMV2)
        CALL MVMU2(TMV2,3,3,TBRTPD1T,TMV2)

        TMV1=TMV1+TMV2

        CALL MVMU2(VAED1T,3,3,BLADE%ELEMEN(NBPLI)%TEB,TMV1)

        CALL GETTDE(TDE,TDED1T,BLADE%ELEMEN(NBPLI)%BTAT,BLADE%ELEMEN(NBPLI)%SPH,&
                                                   BLADE%ELEMEN(NBPLI)%SPH1T,BLADE%ELEMEN(NBPLI)%SV1X,&
                                                   BLADE%ELEMEN(NBPLI)%SW1X,BLADE%ELEMEN(NBPLI)%SV1X1T,&
                                                   BLADE%ELEMEN(NBPLI)%SW1X1T)

        TMV1=VEAE-VAE
        CALL MVMU2(UDP,3,3,TDE,TMV1)
        CALL MVMU2(TMV2,3,3,TDED1T,TMV1)
        CALL MVMU2(UDPD1T,3,3,TDE,VEAED1T-VAED1T)

        UDPD1T=UDPD1T+TMV2
    END SELECT GET_UDP  

    CALL THIS%AEROMODEL(UDP,UDPD1T,NBPLIa,TAF,TAFD1T,TAFD2T,&
                    RB,PSI,SOL,TRN,TBRBP,TBRTP,BLADE%ELEMEN(NBPLI)%TEB,TDE)

    TDET=TRANSPOSE(TDE)
        
    CALL MVMU2(MLCST_ST(:,NBPLIA),3,3,TDET,MLCST_ST(:,NBPLIA))
    CALL MVMU2(FLCST_ST(:,NBPLIA),3,3,TDET,FLCST_ST(:,NBPLIA))
    
    
    CALL MAMU(TER,3,3,3,BLADE%ELEMEN(NBPLI)%TEB,BLADE%ELEMEN(NBPLI)%TBR)
     
    TM1=TRANSPOSE(TER)
      
    CALL MVMU2(MLCST_FR(:,NBPLIA),3,3,TM1,MLCST_ST(:,NBPLIA))
    CALL MVMU2(FLCST_FR(:,NBPLIA),3,3,TM1,FLCST_ST(:,NBPLIA))

END SUBROUTINE FMEL


subroutine GETLINDV(THIS,INDVL,PSI,RB,R,OMG,SOL)
implicit none
CLASS(AERO) :: THIS
INTEGER,INTENT(in) :: SOL
real(rdt),intent(in) :: psi,rb,R,OMG
real(rdt),intent(out) :: indvl(3)

real(rdt) :: psit
integer :: i

INTEGER :: NI,NJ
real(rdt) :: SSI,SSJ,ZERO=0.0

    IF(SOL.EQ.1) THEN
    
        INDVL=(/ZERO,ZERO,-THIS%LAMI*OMG*R/)

    ELSE IF(SOL.EQ.2) THEN
      
        PSIT=MOD(PSI,PI*2)

        call rxlct2(NI, SSI, THIS%nstp, PSIt, THIS%psiv) 
        call rxlct2(NJ, SSJ, NPTS, RB, RLUDAUD) 

        DO I=1,3
            CALL INTP2D(INDVL(I),SSI,SSJ,THIS%bvqa(I,NJ,NI,1),THIS%bvqa(I,NJ,NI+1,1),THIS%bvqa(I,NJ+1,NI,1),THIS%BVQA(I,NJ+1,NI+1,1))
            !CALL INTP2D(INDVL(I),SSI,SSJ,THIS%bvqa(I,NJ,NI,2),THIS%bvqa(I,NJ,NI+1,2),THIS%bvqa(I,NJ+1,NI,2),THIS%BVQA(I,NJ+1,NI+1,2))!更改共轴上下旋翼诱导速度输入
        END DO
    
    END IF

END SUBROUTINE

SUBROUTINE AEROMODEL(THIS,UDP,UDPD1T,LBI,TAF,TAFD1T,TAFD2T,RB,PSI,SOL,TRN,TBRBP,TBRTP,TEB,TDE)
IMPLICIT NONE
CLASS(AERO) :: THIS
REAL(RDT),INTENT(IN) :: UDP(3),UDPD1T(3)
REAL(RDT),INTENT(IN) :: RB,PSI,TAF,TAFD1T,TAFD2T
INTEGER,INTENT(IN) :: SOL 
REAL(RDT),INTENT(IN) :: TRN(3,3),TBRBP(3,3),TBRTP(3,3),TEB(3,3),TDE(3,3)
INTEGER,INTENT(IN) :: LBI 
REAL(RDT) UDR,UDR2,SNAA,CSAA
REAL(RDT) LAF,MAF,DAF,LAFX
REAL(RDT) TX1,TX2,TX3,TX4
REAL(RDT) TMV1(3)
REAL(RDT) :: LBN,LBM,LBC,LBAFA,LBQ,LBMA
REAL(RDT) :: FBG
REAL(RDT) :: FLCSTN(3),FLCSTTP(3),FLCSTNM1(3),MLCSTN(3),MLCSTNM1(3)
INTEGER I,J
REAL(RDT) ::  CFDF(3),CFDM(3),DELTAF(3)
REAL(RDT) :: SS,PSIT,AFA,ZERO=0.0
INTEGER :: NI
REAL(RDT),EXTERNAL:: FWM

    UDR2=UDP(2)**2+UDP(3)**2!参考速度的平方                                 
    UDR=SQRT(UDR2)
    LBMA=UDR/340.0
    AFA=-ATAN(UDP(3)/UDP(2))

    CALL THIS%SEF%TIPLOSS(fbg,THIS%nrbd, rb, THIS%lami)  
    CALL THIS%SEF%GETCL(THIS%ALC,LBMA,AFA)
    CALL THIS%SEF%GETCD(THIS%CD0,LBMA,AFA)

    GET_PDP:SELECT CASE(THIS%MKFM) 
    CASE(-3)
    
        CALL THIS%GREENBERGST(FLCSTN,MLCSTN,LAF,DAF,MAF,UDP,UDPD1T,TAF,TAFD1T,TAFD2T,FBG)
        
    CASE(-2)
        
        CALL THIS%AEROST(FLCSTN,MLCSTN,LAF,DAF,MAF,UDP,FBG,0)
            
    CASE(-1)
        
        CALL THIS%AEROST(FLCSTN,MLCSTN,LAF,DAF,MAF,UDP,FBG,1)
       
!********************************定常********************************************
    CASE(0)
        
        CALL THIS%AEROST(FLCSTN,MLCSTN,LAF,DAF,MAF,UDP,FBG,2)
    
        
!********************************准定常********************************************
    CASE(1)
        
        CALL THIS%GREENBERG(FLCSTN,MLCSTN,LAF,DAF,MAF,UDP,UDPD1T,TAF,TAFD1T,TAFD2T,FBG)

!********************************准定常加静态数据********************************************
    CASE(2)
        TX1=THIS%RHA*THIS%ALC*THIS%BSC
        TX2=UDP(2)-UDP(3)*TAF
        TX3=THIS%XAD-0.5*THIS%BSC
        TX4=TX2*(-UDP(3)+(THIS%BSC-THIS%XAD)*TAFD1T)

        LAF=0.5*THIS%BSC*(-UDPD1T(3)-TX3*TAFD2T)
        LAF=LAF+TX4
        LAF=TX1*LAF

        MAF=-0.5*THIS%BSC*(TX3*UDPD1T(3)+0.5*THIS%BSC*TX2*TAFD1T+(0.125*THIS%BSC**2+TX3**2)*TAFD2T)
        MAF=MAF+THIS%XAD*TX4
        MAF=TX1*MAF
        MAF=0.0

        DAF=THIS%RHA*THIS%CD0*THIS%BSC*UDR2
        
        IF( ABS(UDR) .LT. 1.0D-8 ) THEN 
            LAF=0.0                      !!!  UDR=0 ?
            MAF=0.0
            DAF=0.0
            SNAA=0.0
            CSAA=0.0
            write(*,*) udr
            pause 111
        ELSE
            SNAA=-UDP(3)/UDR              
            CSAA=UDP(2)/UDR
        END IF 


        IF( UDP(2) .LE. 0.0 ) THEN
      
!!!!!!      CALL AEROST(LAF,MAF,DAF,UDP,0)
            MAF=MAF*THIS%RHA*THIS%BSC**2*UDR2*2
            DAF=DAF*THIS%RHA*THIS%BSC*UDR2
            LAF=LAF*THIS%RHA*THIS%BSC*UDR2

        END IF
    
        LAF=LAF*FBG
    
        FLCSTN=(/ZERO,LAF*SNAA-DAF*CSAA,LAF*CSAA+DAF*SNAA/)
        MLCSTN=(/MAF,ZERO,ZERO/)
    

!*****************************非定常***********************************************
    CASE(3)
        IF( ABS(UDR) .LT. 1.0D-8 ) THEN 
            LAF=0.0                      !!!  UDR=0 ?
            MAF=0.0
            DAF=0.0
            SNAA=0.0
            CSAA=0.0
            write(*,*) udr
            pause 111
        ELSE
            SNAA=-UDP(3)/UDR              
            CSAA=UDP(2)/UDR
            LBAFA=-ATAN(UDP(3)/UDP(2))
        END IF 

        IF( UDP(2) .LE. 0.0 ) THEN

!           CALL AEROST(LAF,MAF,DAF,UDP,0)
!            MAF=MAF*RHA*BSCT**2*UDR2*2
!            DAF=DAF*RHA*BSCT*UDR2
!            LAF=LAF*RHA*BSCT*UDR2
!
            DAF=THIS%RHA*THIS%CD0*THIS%BSC*UDR2
            LAF=0.0
            MAF=0.0
            DAF=DAF
            

        ELSE IF(RB.LT.THIS%mu) THEN

            TX1=THIS%RHA*THIS%ALC*THIS%BSC
            TX2=UDP(2)-UDP(3)*TAF
            TX3=THIS%XAD-0.5*THIS%BSC
            TX4=TX2*(-UDP(3)+(THIS%BSC-THIS%XAD)*TAFD1T)

            LAF=0.5*THIS%BSC*(-UDPD1T(3)-TX3*TAFD2T)
            LAF=LAF+TX4
            LAF=TX1*LAF
            

            
            MAF=-0.5*THIS%BSC*(TX3*UDPD1T(3)+0.5*THIS%BSC*TX2*TAFD1T+(0.125*THIS%BSC**2+TX3**2)*TAFD2T)
            MAF=MAF+THIS%XAD*TX4
            MAF=TX1*MAF

            DAF=THIS%RHA*THIS%CD0*THIS%BSC*UDR2


        ELSE        
 
            LBQ=-(UDPD1T(3)*UDP(2)-UDP(3)*UDPD1T(2))/UDR**2

            
            CALL THIS%LB%LBPRE(LBI,0)
            CALL THIS%LB%LBmodel(LBN,LBM,LBC,LBAFA,LBQ,THIS%cchordrf,LBMA,0,0)
            CALL THIS%LB%LBPRE(LBI,1)

            MAF=LBM*THIS%RHA*THIS%BSC**2*UDR**2*2
            DAF=LBC*THIS%RHA*THIS%BSC*UDR**2
            LAF=LBN*THIS%RHA*THIS%BSC*UDR**2
            

        END IF    
        
        LAF=LAF*FBG    
        
        FLCSTN=(/ZERO,LAF*SNAA-DAF*CSAA,LAF*CSAA+DAF*SNAA/)
        MLCSTN=(/MAF,ZERO,ZERO/)
    
        
!*****************************CFD数据***********************************************
    CASE(4)

        IF( ABS(UDR) .LT. 1.0D-8 ) THEN 
            LAF=0.0                      !!!  UDR=0 ?
            MAF=0.0
            DAF=0.0
            SNAA=0.0
            CSAA=0.0
        ELSE
            SNAA=-UDP(3)/UDR              
            CSAA=UDP(2)/UDR
        END IF 
        
    IF(SOL.EQ.1) THEN
    
        CFDF=THIS%CFDFNR(:,LBI)
        CFDM=THIS%CFDMNR(:,LBI)
        CFDM=0.0D0
    ELSE IF(SOL.EQ.2) THEN
 
        PSIT=MOD(PSI,PI*2)

        call rxlct2(NI, SS, THIS%NSTP, PSIt, THIS%PSIV) 
 
        DO I=1,3
            CFDF(I)=FWM(THIS%CFDFPSI(I,LBI,NI),THIS%CFDFPSI(I,LBI,NI+1),SS)
            CFDM(I)=FWM(THIS%CFDMPSI(I,LBI,NI),THIS%CFDMPSI(I,LBI,NI+1),SS)
        END DO
        CFDM=0.0D0
    END IF

            !CALL MVMU2(TMV1,3,3,TRN,CFDF)
            CALL MVMU2(TMV1,3,3,TBRBP,CFDF)
            CALL MVMU2(TMV1,3,3,TBRTP,TMV1)
            CALL MVMU2(TMV1,3,3,TEB,TMV1)
            CALL MVMU2(CFDF,3,3,TDE,TMV1)  
                 
            !CALL MVMU2(TMV1,3,3,TRN,CFDM)
            CALL MVMU2(TMV1,3,3,TBRBP,CFDM)
            CALL MVMU2(TMV1,3,3,TBRTP,TMV1)
            CALL MVMU2(TMV1,3,3,TEB,TMV1)
            CALL MVMU2(CFDM,3,3,TDE,TMV1)  


        MAF=CFDM(1)*0.0
        DAF=THIS%RHA*THIS%CD0*THIS%BSC*UDR2
        LAF=CFDF(3)*(THIS%RHA*THIS%BSC*UDR2)
        
        FLCSTN=(/ZERO,LAF*SNAA-DAF*CSAA,LAF*CSAA+DAF*SNAA/)
        MLCSTN=(/MAF,ZERO,ZERO/)
  CASE(5)
        
        UDR2=UDP(2)**2+UDP(3)**2
        UDR=SQRT(UDR2)
    
        TX1=THIS%RHA*THIS%ALC*THIS%BSC
        TX2=UDP(2)-UDP(3)*TAF
        TX3=THIS%XAD-0.5*THIS%BSC
        TX4=TX2*(-UDP(3)+(THIS%BSC-THIS%XAD)*TAFD1T)

        LAF=0.5*THIS%BSC*(-UDPD1T(3)-TX3*TAFD2T)
        LAF=LAF+TX4
        LAF=TX1*LAF
        
        MAF=-0.5*THIS%BSC*(TX3*UDPD1T(3)+0.5*THIS%BSC*TX2*TAFD1T+(0.125*THIS%BSC**2+TX3**2)*TAFD2T)
        MAF=MAF+THIS%XAD*TX4
        MAF=TX1*MAF

        DAF=THIS%RHA*THIS%CD0*THIS%BSC*UDR2
   
        LAF=LAF*FBG
        
        PSIT=MOD(PSI,PI*2)
        call rxlct2(NI, SS, THIS%NSTP, PSIt, THIS%PSIV) 

        DO I=1,3
            CFDF(I)=FWM(THIS%CFDFPSI(I,LBI,NI),THIS%CFDFPSI(I,LBI,NI+1),SS)*THIS%RHA*THIS%BSC*UDR**2
            CFDM(I)=FWM(THIS%CFDMPSI(I,LBI,NI),THIS%CFDMPSI(I,LBI,NI+1),SS)*THIS%RHA*THIS%BSC**2*UDR**2*2
        END DO
        
        !CFDM=0.0
        !DELTAF=(/ZERO,-DAF,LAF/)+CFDF-THIS%FLCSTCOE(1:3,LBI)*THIS%RHA*THIS%BSC*UDR**2
        
        !FLCSTN=(/ZERO,-DELTAF(2),DELTAF(3)/)
        !MLCSTN=(/CFDM(1),ZERO,ZERO/)
		
        FLCSTN=(/ZERO,CFDF(2),CFDF(3)/)
        MLCSTN=(/CFDM(1),ZERO,ZERO/)
  CASE(7)
        
        UDR2=UDP(2)**2+UDP(3)**2
        UDR=SQRT(UDR2)
    
        TX1=THIS%RHA*THIS%ALC*THIS%BSC
        TX2=UDP(2)-UDP(3)*TAF
        TX3=THIS%XAD-0.5*THIS%BSC
        TX4=TX2*(-UDP(3)+(THIS%BSC-THIS%XAD)*TAFD1T)

        LAF=0.5*THIS%BSC*(-UDPD1T(3)-TX3*TAFD2T)
        LAF=LAF+TX4
        LAF=TX1*LAF
        
        MAF=-0.5*THIS%BSC*(TX3*UDPD1T(3)+0.5*THIS%BSC*TX2*TAFD1T+(0.125*THIS%BSC**2+TX3**2)*TAFD2T)
        MAF=MAF+THIS%XAD*TX4
        MAF=TX1*MAF

        DAF=THIS%RHA*THIS%CD0*THIS%BSC*UDR2
   
        LAF=LAF*FBG
        
        PSIT=MOD(PSI,PI*2)
        call rxlct2(NI, SS, THIS%NSTP, PSIt, THIS%PSIV) 

        DO I=1,3
            CFDF(I)=FWM(THIS%CFDFPSI(I,LBI,NI),THIS%CFDFPSI(I,LBI,NI+1),SS)*THIS%RHA*THIS%BSC*UDR**2
            CFDM(I)=FWM(THIS%CFDMPSI(I,LBI,NI),THIS%CFDMPSI(I,LBI,NI+1),SS)*THIS%RHA*THIS%BSC**2*UDR**2*2
        END DO
        
		CALL MVMU2(TMV1,3,3,TBRBP,CFDF)
		CALL MVMU2(TMV1,3,3,TBRTP,TMV1)
		CALL MVMU2(TMV1,3,3,TEB,TMV1)
		CALL MVMU2(CFDF,3,3,TDE,TMV1)  
		 
		!CALL MVMU2(TMV1,3,3,TBRBP,CFDM)
		!CALL MVMU2(TMV1,3,3,TBRTP,TMV1)
		!CALL MVMU2(TMV1,3,3,TEB,TMV1)
		!CALL MVMU2(CFDM,3,3,TDE,TMV1)  
		
        !CFDM=0.0
        !DELTAF=(/ZERO,-DAF,LAF/)+CFDF-THIS%FLCSTCOE(1:3,LBI)*THIS%RHA*THIS%BSC*UDR**2
        
        !FLCSTN=(/ZERO,-DELTAF(2),DELTAF(3)/)
        !MLCSTN=(/CFDM(1),ZERO,ZERO/)
		
        FLCSTN=(/ZERO,CFDF(2),CFDF(3)/)
        MLCSTN=(/MAF,ZERO,ZERO/)	
 CASE(6)
        !IF( ABS(UDR) .LT. 1.0D-8 ) THEN !增加SNAA和CSAA的计算
        !    LAF=0.0                      !!!  UDR=0 ?
        !    MAF=0.0
        !    DAF=0.0
        !    SNAA=0.0
        !    CSAA=0.0
        !ELSE
        !    SNAA=-UDP(3)/UDR              
        !    CSAA=UDP(2)/UDR
        !END IF
 
		!IF(RB.LT.1.0) THEN
		!	CALL THIS%GREENBERG(FLCSTN,MLCSTN,LAF,DAF,MAF,UDP,UDPD1T,TAF,TAFD1T,TAFD2T,FBG)
		!else
        
        PSIT=MOD(PSI,PI*2)
        call rxlct2(NI, SS, THIS%NSTP, PSIt, THIS%PSIV) 

        DO I=1,3
            !CFDF(I)=FWM(THIS%CFDFPSI(I,LBI,NI),THIS%CFDFPSI(I,LBI,NI+1),SS)*THIS%RHA*THIS%BSC*UDR**2
            !CFDM(I)=FWM(THIS%CFDMPSI(I,LBI,NI),THIS%CFDMPSI(I,LBI,NI+1),SS)*THIS%RHA*THIS%BSC**2*UDR**2*2
            CFDF(I)=FWM(THIS%CFDFPSI(I,LBI,NI),THIS%CFDFPSI(I,LBI,NI+1),SS)
            CFDM(I)=FWM(THIS%CFDMPSI(I,LBI,NI),THIS%CFDMPSI(I,LBI,NI+1),SS)
        END DO
        
        
        
        !CFDM=0.0
        !DELTAF=(/ZERO,-DAF,LAF/)+CFDF-THIS%FLCSTCOE(1:3,LBI)*THIS%RHA*THIS%BSC*UDR**2
        
        !FLCSTN=(/ZERO,-DELTAF(2),DELTAF(3)/)
        !MLCSTN=(/CFDM(1),ZERO,ZERO/)
		
        !LAF*SNAA-DAF*CSAA,LAF*CSAA+DAF*SNAA
        FLCSTN=(/ZERO,-CFDF(3),CFDF(2)/)
        MLCSTN=(/CFDM(1),ZERO,ZERO/)
 
		!end if
 
 
 
        !IF( ABS(UDR) .LT. 1.0D-8 ) THEN 
        !    LAF=0.0                      !!!  UDR=0 ?
        !    MAF=0.0
        !    DAF=0.0
        !    SNAA=0.0
        !    CSAA=0.0
        !    write(*,*) udr
        !    pause 111
        !ELSE
        !    SNAA=-UDP(3)/UDR              
        !    CSAA=UDP(2)/UDR
        !    LBAFA=-ATAN(UDP(3)/UDP(2))
        !END IF 

        !IF( UDP(2) .LE. 0.0 ) THEN

        !    DAF=THIS%RHA*THIS%CD0*THIS%BSC*UDR2
        !    LAF=0.0
        !    MAF=0.0
        !    DAF=DAF

        !ELSE IF(RB.LT.THIS%mu) THEN

        !    TX1=THIS%RHA*THIS%ALC*THIS%BSC
        !    TX2=UDP(2)-UDP(3)*TAF
        !    TX3=THIS%XAD-0.5*THIS%BSC
        !    TX4=TX2*(-UDP(3)+(THIS%BSC-THIS%XAD)*TAFD1T)

        !    LAF=0.5*THIS%BSC*(-UDPD1T(3)-TX3*TAFD2T)
        !    LAF=LAF+TX4
        !    LAF=TX1*LAF
        !    
        !    MAF=-0.5*THIS%BSC*(TX3*UDPD1T(3)+0.5*THIS%BSC*TX2*TAFD1T+(0.125*THIS%BSC**2+TX3**2)*TAFD2T)
        !    MAF=MAF+THIS%XAD*TX4
        !    MAF=TX1*MAF

        !    DAF=THIS%RHA*THIS%CD0*THIS%BSC*UDR2

        !ELSE        
       
		!	PSIT=MOD(PSI,PI*2)
		!	call rxlct2(NI, SS, THIS%NSTP, PSIt, THIS%PSIV) 
		!	DO I=1,3
		!		CFDF(I)=FWM(THIS%CFDFPSI(I,LBI,NI),THIS%CFDFPSI(I,LBI,NI+1),SS)*THIS%RHA*THIS%BSC*UDR**2
		!		CFDM(I)=FWM(THIS%CFDMPSI(I,LBI,NI),THIS%CFDMPSI(I,LBI,NI+1),SS)*THIS%RHA*THIS%BSC**2*UDR**2*2
		!	END DO
		!	
        !    TX1=THIS%RHA*THIS%ALC*THIS%BSC
        !    TX2=UDP(2)-UDP(3)*TAF
        !    TX3=THIS%XAD-0.5*THIS%BSC
        !    TX4=TX2*(-UDP(3)+(THIS%BSC-THIS%XAD)*TAFD1T)

        !    LAF=0.5*THIS%BSC*(-UDPD1T(3)-TX3*TAFD2T)
        !    LAF=LAF+TX4
        !    LAF=TX1*LAF
        !    
        !    MAF=-0.5*THIS%BSC*(TX3*UDPD1T(3)+0.5*THIS%BSC*TX2*TAFD1T+(0.125*THIS%BSC**2+TX3**2)*TAFD2T)
        !    MAF=MAF+THIS%XAD*TX4
        !    MAF=TX1*MAF
		!	
        !    DAF=THIS%RHA*THIS%CD0*THIS%BSC*UDR2
		!	
        !    LAF=CFDF(3)
        !    DAF=CFDF(2)
        !    MAF=CFDM(1)
		!END IF

        !FLCSTN=(/ZERO,LAF*SNAA-DAF*CSAA,LAF*CSAA+DAF*SNAA/)
        !MLCSTN=(/MAF,ZERO,ZERO/)
		
!*********************************************************************************
    END SELECT GET_PDP  
    
100 CONTINUE        

    FLCST_ST(:,LBI)=FLCSTN
    MLCST_ST(:,LBI)=MLCSTN
    
    !THIS%FLCSTCOE(:,LBI)=(/ZERO,-DAF,LAF/)
    !THIS%MLCSTCOE(:,LBI)=(/MAF,ZERO,ZERO/)
    !THIS%FLCSTCOE(:,LBI)=FLCSTN
    !THIS%MLCSTCOE(:,LBI)=MLCSTN
    !THIS%FLCSTCOE(:,LBI)=(/THIS%CFDFPSI(1:3,LBI,NI+1),(/ZERO,-DAF,LAF/)/(THIS%RHA*THIS%BSC*UDR2)*LBMA**2/)
    !THIS%FLCSTCOE(:,LBI)=(/(/ZERO,-DAF,LAF/)/(THIS%RHA*THIS%BSC*UDR2),(/ZERO,-DELTAF(2),DELTAF(3)/)/(THIS%RHA*THIS%BSC*UDR2)/)

    !THIS%FLCSTCOE(:,LBI)=(/(/ZERO,-DAF,LAF/)/(THIS%RHA*THIS%BSC*UDR2),(/ZERO,-DAF,LAF/)/(THIS%RHA*THIS%BSC*UDR2)*LBMA**2/)
    !THIS%FLCSTCOE(:,LBI)=(/FLCSTN(1:3)/(THIS%RHA*THIS%BSC*UDR2),(/ZERO,-DELTAF(2),DELTAF(3)/)/(THIS%RHA*THIS%BSC*UDR2)/)
    THIS%FLCSTCOE(:,LBI)=(/FLCSTN(1:3)/(THIS%RHA*THIS%BSC*UDR2),       (/FLCSTN(1:3)/)/(THIS%RHA*THIS%BSC*UDR2)*LBMA**2/)!二维力无量纲化
    THIS%MLCSTCOE(:,LBI)=(/MLCSTN(1:3)/(THIS%RHA*THIS%BSC**2*UDR**2*2),(/MLCSTN(1:3)/)/(THIS%RHA*THIS%BSC**2*UDR**2*2)*LBMA**2/)

!    FLCSTCOE(:,LBI)=FLCSTN/(RHA*BSCT*UDR2)
!    MLCSTCOE(:,LBI)=MLCSTN/(RHA*BSCT**2*UDR**2*2)
  
!    FLCST(2:3,LBI)=(/-DAF,LAF/)/(RHA*BSCT*UDR2)*(udr/340)**2

!!    FLCST(:,LBI)=pdp/(RHA*BSCT*UDR2)
!!    MLCST(:,LBI)=(/MAF,0.0,0.0/)/(RHA*BSCT**2*UDR**2*2)*(udr/340)**2
!   
!
!    AOALCST(LBI)=-ATAN(UDP(3)/UDP(2))*180/PI
        993 format(1000f26.12)
END SUBROUTINE

!*********************************************************************************
SUBROUTINE GREENBERG(THIS,FLCSTN,MLCSTN,LAF,DAF,MAF,UDP,UDPD1T,TAF,TAFD1T,TAFD2T,FBG)
IMPLICIT NONE
CLASS(AERO) :: THIS
REAL(RDT),INTENT(IN) :: UDP(3),UDPD1T(3)
REAL(RDT),INTENT(IN) :: TAF,TAFD1T,TAFD2T
REAL(RDT),INTENT(OUT) ::  FLCSTN(3),MLCSTN(3),LAF,DAF,MAF

REAL(RDT) UDR,UDR2,SNAA,CSAA
REAL(RDT) TX1,TX2,TX3,TX4
REAL(RDT) :: FBG,ZERO=0.0,xx

        UDR2=UDP(2)**2+UDP(3)**2
        UDR=SQRT(UDR2)
    
        TX1=THIS%RHA*THIS%ALC*THIS%BSC!Greenberg模型中a·ρA·b项
        TX2=UDP(2)-UDP(3)*TAF!Greenberg模型中来流速度V
        TX3=THIS%XAD-0.5*THIS%BSC!greenberg模型中升力计算公式xA-b/2项。
        TX4=TX2*(-UDP(3)+(THIS%BSC-THIS%XAD)*TAFD1T)!其中TAFD1T：来流迎角对时间的一阶导数；

        LAF=0.5*THIS%BSC*(-UDPD1T(3)-TX3*TAFD2T)
        LAF=LAF+TX4
        LAF=TX1*LAF
        
        MAF=-0.5*THIS%BSC*(TX3*UDPD1T(3)+0.5*THIS%BSC*TX2*TAFD1T+(0.125*THIS%BSC**2+TX3**2)*TAFD2T)
        MAF=MAF+THIS%XAD*TX4
        MAF=TX1*MAF

        DAF=THIS%RHA*THIS%CD0*THIS%BSC*UDR2!阻力用阻力系数计算，精度低
        
        IF( ABS(UDR) .LT. 1.0D-8 ) THEN 
            LAF=0.0                      !!!  UDR=0 ?
            MAF=0.0
            DAF=0.0
            SNAA=0.0
            CSAA=0.0
        ELSE
            SNAA=-UDP(3)/UDR              
            CSAA=UDP(2)/UDR
        END IF 

        IF( UDP(2) .LE. 0.0 ) THEN
            LAF=0.0
            MAF=0.0
            DAF=DAF
        END IF
    
        LAF=LAF*FBG!FBG：桨尖损失相关

        FLCSTN=(/ZERO,LAF*SNAA-DAF*CSAA,LAF*CSAA+DAF*SNAA/)
        MLCSTN=(/MAF,ZERO,ZERO/)

END SUBROUTINE
!*********************************************************************************
SUBROUTINE GREENBERGST(THIS,FLCSTN,MLCSTN,LAF,DAF,MAF,UDP,UDPD1T,TAF,TAFD1T,TAFD2T,FBG)
IMPLICIT NONE
CLASS(AERO) :: THIS
REAL(RDT),INTENT(IN) :: UDP(3),UDPD1T(3)
REAL(RDT),INTENT(IN) :: TAF,TAFD1T,TAFD2T
REAL(RDT),INTENT(OUT) ::  FLCSTN(3),MLCSTN(3),LAF,DAF,MAF

REAL(RDT) UDR,UDR2,SNAA,CSAA
REAL(RDT) TX1,TX2,TX3,TX4
REAL(RDT) AFA,LBMA,CLF,CMF,CDF
REAL(RDT) :: FBG,ZERO=0.0

        UDR2=UDP(2)**2+UDP(3)**2
        UDR=SQRT(UDR2)
        
        LBMA=UDR/340.0
        CALL GETAFA(AFA,UDP)
        CALL THIS%SEF%AEROST2(CLF,CMF,CDF,THIS%ALC,THIS%CD0,AFA,LBMA) 
        
        THIS%ALC=CLF
        
        TX1=THIS%RHA*THIS%ALC*THIS%BSC
        TX2=UDP(2)-UDP(3)*TAF
        TX3=THIS%XAD-0.5*THIS%BSC
        TX4=TX2*(-UDP(3)+(THIS%BSC-THIS%XAD)*TAFD1T)

        LAF=0.5*THIS%BSC*(-UDPD1T(3)-TX3*TAFD2T)
        LAF=LAF+TX4
        LAF=TX1*LAF
        
        MAF=-0.5*THIS%BSC*(TX3*UDPD1T(3)+0.5*THIS%BSC*TX2*TAFD1T+(0.125*THIS%BSC**2+TX3**2)*TAFD2T)
        MAF=MAF+THIS%XAD*TX4
        MAF=TX1*MAF

        DAF=THIS%RHA*THIS%CD0*THIS%BSC*UDR2
        
        IF( ABS(UDR) .LT. 1.0D-8 ) THEN 
            LAF=0.0                      !!!  UDR=0 ?
            MAF=0.0
            DAF=0.0
            SNAA=0.0
            CSAA=0.0
            write(*,*) udr
            pause 111
        ELSE
            SNAA=-UDP(3)/UDR              
            CSAA=UDP(2)/UDR
        END IF 

        IF( UDP(2) .LE. 0.0 ) THEN
            LAF=0.0
            MAF=0.0
            DAF=DAF
        END IF
    
        LAF=LAF*FBG

        FLCSTN=(/ZERO,LAF*SNAA-DAF*CSAA,LAF*CSAA+DAF*SNAA/)
        MLCSTN=(/MAF,ZERO,ZERO/)

END SUBROUTINE
!SUBROUTINE STAERO(FLCSTN,MLCSTN,RHA,ALCT,CD0T,BSCT,UDP,UDPD1T,TAF,TAFD1T,TAFD2T,XADT,FBG)
!USE GlobalDataFun
!IMPLICIT NONE
!REAL(RDT),INTENT(IN) :: UDP(3),UDPD1T(3)
!REAL(RDT),INTENT(IN) :: RHA, ALCT,CD0T
!REAL(RDT),INTENT(IN) :: BSCT,XADT,TAF,TAFD1T,TAFD2T
!REAL(RDT),INTENT(OUT) ::  FLCSTN(3),MLCSTN(3)
!
!REAL(RDT) UDR,UDR2,SNAA,CSAA
!REAL(RDT) LAF,MAF,DAF
!REAL(RDT) TX1,AFA
!REAL(RDT) :: FBG
!
!        AFA=-ATAN(UDP(3)/UDP(2))
!        UDR2=UDP(2)**2+UDP(3)**2
!       UDR=SQRT(UDR2)
!   
!       TX1=RHA*ALCT*BSCT*UDR2
!       LAF=TX1*AFA
!        
!       MAF=0.0
!
!       DAF=RHA*CD0T*BSCT*UDR2
!       
!        IF( ABS(UDR) .LT. 1.0D-8 ) THEN 
!            LAF=0.0                         !!!  UDR=0 ?
!            MAF=0.0
!            DAF=0.0
!            SNAA=0.0
!            CSAA=0.0
!            write(*,*) udr
!            pause 111
!        ELSE
!            SNAA=-UDP(3)/UDR             
!            CSAA=UDP(2)/UDR
!        END IF 
!
!        IF( UDP(2) .LE. 0.0 ) THEN
!            LAF=0.0
!            MAF=0.0
!            DAF=DAF
!        END IF
!   
!        LAF=LAF*FBG
!
!        FLCSTN=(/0.0,LAF*SNAA-DAF*CSAA,LAF*CSAA+DAF*SNAA/)
!        MLCSTN=(/MAF,0.0,0.0/)
!
!END
!
!SUBROUTINE STST(FLCSTN,MLCSTN,RHA,ALCT,CD0T,BSCT,UDP,UDPD1T,TAF,TAFD1T,TAFD2T,XADT,FBG)
!USE GlobalDataFun
!USE SECTIONF
!IMPLICIT NONE
!REAL(RDT),INTENT(IN) :: UDP(3),UDPD1T(3)
!REAL(RDT),INTENT(INOUT) :: RHA, ALCT,CD0T
!REAL(RDT),INTENT(IN) :: BSCT,XADT,TAF,TAFD1T,TAFD2T
!REAL(RDT),INTENT(OUT) ::  FLCSTN(3),MLCSTN(3)
!
!REAL(RDT) UDR,UDR2,SNAA,CSAA
!REAL(RDT) LAF,MAF,DAF,CM,CL,CD
!REAL(RDT) TX1,AFA,LBMA
!REAL(RDT) :: FBG
!
!        CALL GETAFA(AFA,UDP)
!        
!        UDR2=UDP(2)**2+UDP(3)**2
!       UDR=SQRT(UDR2)
!        LBMA=UDR/340.0
!        
!
!       IF( ABS(UDR) .LT. 1.0D-8 ) THEN 
!            LAF=0.0                         !!!  UDR=0 ?
!            MAF=0.0
!            DAF=0.0
!            SNAA=0.0
!            CSAA=0.0
!            write(*,*) udr
!            pause 111
!        ELSE
!            SNAA=-UDP(3)/UDR             
!            CSAA=UDP(2)/UDR
!        END IF 
!
!!        IF( UDP(2) .LE. 0.0 ) THEN
!!            LAF=0.0
!!            MAF=0.0
!!            DAF=DAF
!!        END IF
!   
!        LAF=LAF*FBG
!
!        FLCSTN=(/0.0,LAF*SNAA-DAF*CSAA,LAF*CSAA+DAF*SNAA/)
!        MLCSTN=(/MAF,0.0,0.0/)
!
!END
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
SUBROUTINE   AEROST(THIS,FLCSTN,MLCSTN,LAF,DAF,MAF,UDP,FBG,CHO)
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
implicit none
CLASS(AERO) :: THIS
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
INTEGER,INTENT(IN) :: CHO
REAL(RDT),INTENT(IN) :: UDP(3)
REAL(RDT),INTENT(IN) :: FBG
REAL(RDT),INTENT(OUT) ::  FLCSTN(3),MLCSTN(3),LAF,DAF,MAF
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->       
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->       READDATA1, READDATA2
!->       EIGANS, HOVER
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
REAL(RDT) UDR,UDR2,SNAA,CSAA
real(rdt) :: CLF,CMF,CDF
real(rdt) :: PI,AFA,LBMA,ZERO=0.0


        UDR2=UDP(1)**2+UDP(2)**2+UDP(3)**2
        UDR=SQRT(UDR2)
        LBMA=UDR/340.0
        PI=ATAN(1.0D0)*4

        IF( ABS(UDR) .LT. 1.0D-8 ) THEN 
            write(*,*) udr
            pause 111
        END IF 
        
    CALL GETAFA(AFA,UDP)
        
      
    IF(CHO.EQ.0) THEN
         AFA=AFA/PI*180
        CALL AEROST0(CLF,CMF,CDF,AFA,LBMA)
        AFA=AFA*PI/180
    ELSE IF(CHO.EQ.1) THEN
         AFA=AFA/PI*180
        CALL AEROST1(CLF,CMF,CDF,AFA,LBMA)
        AFA=AFA*PI/180
    ELSE IF(CHO.EQ.2) THEN
        CALL THIS%SEF%AEROST2(CLF,CMF,CDF,THIS%ALC,THIS%CD0,AFA,LBMA) 
        CLF=CLF*AFA
        cmf=0.0
!        write(*,997) clf,cmf,cdf,afa,lbma
    END IF    
    
    MAF=CMF*THIS%RHA*THIS%BSC**2*UDR2*2
    DAF=CDF*THIS%RHA*THIS%BSC*UDR2
    LAF=CLF*THIS%RHA*THIS%BSC*UDR2
    
        IF( ABS(UDR) .LT. 1.0D-8 ) THEN 
            LAF=0.0D0                        !!!  UDR=0 ?
            MAF=0.0D0
            DAF=0.0D0
            SNAA=0.0D0
            CSAA=0.0D0
            write(*,*) udr
            pause 111
        ELSE
            SNAA=-UDP(3)/UDR              
            CSAA=UDP(2)/UDR
        END IF 

!        IF( UDP(2) .LE. 0.0 ) THEN
!            LAF=0.0D0
!            MAF=0.0D0
!            DAF=DAF
!        END IF
    
        LAF=LAF*FBG

        FLCSTN=(/ZERO,LAF*SNAA-DAF*CSAA,LAF*CSAA+DAF*SNAA/)
        MLCSTN=(/MAF,ZERO,ZERO/)

997 format(20f7.3)
END SUBROUTINE




!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
SUBROUTINE AEROST0(CLF,CMF,CDF,AFA,LBMA)
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
implicit none
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
real(rdt),intent(in) :: AFA,LBMA
real(rdt),intent(out) :: CLF,CMF,CDF
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->       
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->       READDATA1, READDATA2
!->       EIGANS, HOVER
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES
real(rdt) :: SSI,SSJ,VIM1JM1,VIJM1,VIM1J,VIJ
character(len=100) :: copyfilenam

!    write(copyfilenam,*) trim(adjustl(inpre)), 'cl.txt' 
!    OPEN(10,FILE=copyfilenam)
!    write(copyfilenam,*) trim(adjustl(inpre)), 'cd.txt' 
!    OPEN(11,FILE=copyfilenam)
!    write(copyfilenam,*) trim(adjustl(inpre)), 'cm.txt' 
!    OPEN(12,FILE=copyfilenam)
    pause

    CALL GETSS(SSI,SSJ,VIM1JM1,VIJM1,VIJM1,VIJ,AFA,LBMA,10)
    CALL INTP2D(CLF,SSI,SSJ,VIM1JM1,VIJM1,VIM1J,VIJ)
    
    
    CALL GETSS(SSI,SSJ,VIM1JM1,VIJM1,VIJM1,VIJ,AFA,LBMA,11)
    CALL INTP2D(CDF,SSI,SSJ,VIM1JM1,VIJM1,VIM1J,VIJ)
    
    
    CALL GETSS(SSI,SSJ,VIM1JM1,VIJM1,VIJM1,VIJ,AFA,LBMA,12)
    CALL INTP2D(CMF,SSI,SSJ,VIM1JM1,VIJM1,VIM1J,VIJ)


    
    CLOSE(10)
    CLOSE(11)
    CLOSE(12) 
         
END SUBROUTINE

SUBROUTINE GETSS(SSI,SSJ,VIM1JM1,VIJM1,VIM1J,VIJ,AFA,LBMA,FNO)
IMPLICIT NONE
INTEGER,INTENT(IN) :: FNO
real(rdt),INTENT(IN) :: AFA,LBMA
real(rdt),INTENT(OUT) :: SSI,SSJ
real(rdt),INTENT(OUT) :: VIM1JM1,VIJM1,VIM1J,VIJ

real(rdt),ALLOCATABLE :: PSI(:),MA(:),CXTAB(:,:)
INTEGER :: I,IM1,J,JM1,M,N

    READ(FNO,*) M,N
    ALLOCATE(PSI(M),MA(N),CXTAB(M,N))
    
    
    READ(FNO,*) (MA(I),I=1,N)

    DO I=1,M
        READ(FNO,*) PSI(I),(CXTAB(I,J),J=1,N)
    END DO
    
    IF((LBMA.LT.0.0D0).OR.(LBMA.GT.1.0D0)) THEN
        WRITE(*,*) "ERROR IN GETSS"
        STOP
    END IF
    
    J=1
    DO WHILE(MA(J).LT.LBMA)
        J=J+1
    END DO
    JM1=J-1

    
    IF((AFA.LT.-180.0D0).OR.(AFA.GT.180.0D0)) THEN
        WRITE(*,*) "ERROR IN GETSS"
        STOP
    END IF
    
    IF(ABS(ABS(AFA)-180.0) .LT.1E-6) THEN
        I=1
        IM1=M
    ELSE 
        I=1
        DO WHILE(PSI(I).LT.AFA)
            I=I+1
        END DO
        IM1=I-1
    END IF
    
    VIM1JM1=CXTAB(IM1,JM1)
    VIM1J=CXTAB(IM1,J)
    VIJM1=CXTAB(I,JM1)
    VIJ=CXTAB(I,J)
     
    SSJ=(LBMA-MA(JM1))/(MA(J)-MA(JM1))
    SSI=(AFA-PSI(IM1))/(PSI(I)-PSI(IM1))
    
    DEALLOCATE(PSI,MA,CXTAB)
    
END SUBROUTINE


!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
SUBROUTINE AEROST1(CLF,CMF,CDF,AFA,LBMA)
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
implicit none
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
real(rdt),intent(in) :: AFA,LBMA
real(rdt),intent(out) :: CLF,CMF,CDF
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->       
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->       READDATA1, READDATA2
!->       EIGANS, HOVER
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:

real(rdt) :: SIG,AFAT
real(rdt) :: CLAFA,AFAL,K1,K2
real(rdt) :: CDI

        IF(AFA.LT.0) THEN
            SIG=-1
        ELSE
            SIG=1
        END IF
        
        AFAT=ABS(AFA)

        IF((AFAT.GE.0).AND.(AFAT.LE.15)) THEN
        
            IF(LBMA.LE.0.725) THEN
                CLAFA=0.1/SQRT(1-LBMA**2)-0.01*LBMA
                AFAL=15-16*LBMA
                IF(AFAT.LT.AFAL) THEN
                    CLF=CLAFA*AFAT
                ELSE 
                    K1=0.0233+0.342*LBMA**7.15
                    K2=2.05-0.95*LBMA
                    CLF=CLAFA*AFAT-K1*(AFAT-AFAL)**K2
                END IF
            ELSE IF(LBMA.GT.0.725) THEN
                CLAFA=0.677-0.744*LBMA
                AFAL=3.4
                IF(AFAT.LT.AFAL) THEN
                    CLF=CLAFA*AFAT
                ELSE 
                    K1=0.0575-0.144*(LBMA-0.725)**0.44
                    K2=2.05-0.95*LBMA
                    CLF=CLAFA*AFAT-K1*(AFAT-AFAL)**K2
                END IF                
            END IF
            
        ELSE IF((AFAT.GT.15).AND.(AFAT.LE.161)) THEN
        
            CLF=1.15*SIND(2*AFAT) 
            
        ELSE IF((AFAT.GT.161).AND.(AFAT.LE.173)) THEN
        
            CLF=-0.7
            
        ELSE IF((AFAT.GT.173).AND.(AFAT.LE.180)) THEN
        
            CLF=0.1*(AFAT-180)
            
        END IF
        
        CLF=SIG*CLF
        

        IF(AFAT.LE.10) THEN
            CDI=0.0081+(65.8*AFAT**2-0.226*AFAT**4+0.0046*AFAT**6)*10E-6
            IF(LBMA.LE.0.725) THEN
                AFAL=17-23.4*LBMA
                IF(AFAT.LT.AFAL) THEN
                    CDF=CDI
                ELSE 
                    K1=0.00066
                    K2=2.54
                    CDF=CDI+K1*(AFAT-AFAL)**K2
                END IF
            ELSE IF(LBMA.GT.0.725) THEN
                AFAL=0
                K1=0.00035
                K2=2.54
                CDF=CDI+K1*(AFAT-AFAL)**K2
                K1=21
                K2=3.2
                CDF=CDF+K1*(LBMA-0.725)**K2  
            END IF
        ELSE IF((AFAT.GT.10).AND.(AFAT.LE.180)) THEN
            CDF=1.03-1.02*COSD(2*AFAT)   
        END IF

        CMF=0.0D0
        
         
END SUBROUTINE

SUBROUTINE GETAFA(AFA,UDP)
IMPLICIT NONE
REAL(RDT),INTENT(IN) :: UDP(3)
REAL(RDT),INTENT(OUT) :: AFA

        IF( UDP(2) .GT. 0.0 ) THEN
            AFA=-ATAN(UDP(3)/UDP(2))
        ELSE IF( UDP(2) .LT. 0.0 ) THEN
            IF( UDP(3) .LT. 0.0) THEN
                AFA=PI-ATAN(UDP(3)/UDP(2))
            ELSE
                AFA=-PI-ATAN(UDP(3)/UDP(2))
            END IF
        END IF  
        
END SUBROUTINE


SUBROUTINE TECPLOT_BVQA(THIS)
IMPLICIT NONE
CLASS(AERO) :: THIS
INTEGER :: I,J,K

   OPEN(444,FILE='induced_vel.dat')
    WRITE(444,*) 'ZONE I=',NPTS,',J=',THIS%NSTP
    DO I=1,THIS%NSTP
        DO J=1,NPTS
            WRITE(444,"(1000F)") RLUDAUD(J)*SIN(THIS%PSIV(I)),-RLUDAUD(J)*COS(THIS%PSIV(I)),&
                                 (THIS%BVQA(K,J,I,1),K=1,3),(THIS%BVQA(K,J,I,2),K=1,3)!共轴
        END DO
    END DO
    close(444)
	
END SUBROUTINE
  
END MODULE
