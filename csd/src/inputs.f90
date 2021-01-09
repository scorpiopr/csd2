    MODULE INPUTS_CLASS
    USE GlobalDataFun
    IMPLICIT NONE
    !->+++++++++++++++++++++++++++++++++++++++++++++
    TYPE,PUBLIC :: IPTIN

        INTEGER :: NumofRotors,RotorID
        REAL(8) :: THETA(6)!共轴上下旋翼总距和周期变距
        INTEGER :: NLG,NPTS
        REAL(RDT) :: BTAP
        REAL(RDT) :: LDTIP(3)
        REAL(RDT) :: MU
        REAL(RDT) :: thetapa, delta
        INTEGER :: ISR
        REAL(RDT) :: CTC
        REAL(RDT) :: RHA
        REAL(RDT) :: TH0,TH1C,TH1S
        INTEGER :: NOD,RHS,NMN,NMNV,NMNW,NMNT,NMNA,NMNF
        INTEGER :: CNDN,HINGES(6),NFPCT,NHFLP,NPTCT,NHPIT,NLLCT,NHLAG
        INTEGER :: SOL
        INTEGER :: nrbd, afsen
        REAL(RDT) :: lockn, cchordrf, alcsrf
        REAL(RDT) :: kfp, kptbar, klgdpr,kfpspr
        REAL(RDT) :: cptbar, clgdpr,cfpspr
        REAL(RDT) :: ggt
        REAL(RDT) :: rdcm, rdck
        INTEGER :: mkim, mdyim
        INTEGER :: mktrim, mkfm


    CONTAINS

    PROCEDURE,PUBLIC :: READIPT

    END TYPE IPTIN
    !->+++++++++++++++++++++++++++++++++++++++++++++
    TYPE,PUBLIC :: egnin

        INTEGER :: rotan,fch,soc
        REAL(RDT) :: TFLCTRL

    CONTAINS

    PROCEDURE,PUBLIC :: READEGN

    END TYPE egnin
    !->+++++++++++++++++++++++++++++++++++++++++++++
    TYPE,PUBLIC :: fwdin

        INTEGER natmod, mrksl
        real(rdt) beta, gam
        integer sel, nit
        integer mktfe, npl, ntfe, nidtfe
        REAL(RDT) :: hdg
        INTEGER :: ndgop,nrd

    CONTAINS

    PROCEDURE,PUBLIC :: READfwd

    END TYPE fwdin
    !->+++++++++++++++++++++++++++++++++++++++++++++
    TYPE,PUBLIC :: secin

        REAL(RDT) :: MDIV,CLMAX,FSTALL,MSTALL,GSTALL,HSTALL
        REAL(RDT) :: CLF,CMAC,CMS,DEL0,DEL1,DEL2,DCDDM,MCRIT,ACRIT,ALFD,CDF
        integer  :: OPTIP,MACRT
        REAL(RDT)  :: BTIP

    CONTAINS

    PROCEDURE,PUBLIC :: READsec

    END TYPE secin
    !->+++++++++++++++++++++++++++++++++++++++++++++
    TYPE,PUBLIC :: trmin

        integer :: isinput
        REAL(RDT) :: cw
        real(rdt) :: rlxr, rlxv
        real(rdt) :: fcdfpdka
        real(rdt) :: rcgb(3), racb(3)
        real(rdt) :: fchbtg(3), mthbtg(3)

    CONTAINS

    PROCEDURE,PUBLIC :: READtrm

    END TYPE trmin
    !->+++++++++++++++++++++++++++++++++++++++++++++
    TYPE,PUBLIC :: ifwin

        INTEGER  UNITYPE,LINTYPE
        REAL(RDT) KHLMDA,KFLMDA
        real(rdt) FXLMDA,FYLMDA

    CONTAINS

    PROCEDURE,PUBLIC :: READifw

    END TYPE ifwin
    !->+++++++++++++++++++++++++++++++++++++++++++++
    TYPE,PUBLIC :: vtxin

        real(rdt) :: ri0
        integer ::nf,nw,nseg,nx
        real(rdt) :: vor_rootr,vor_tipr,rcct
        INTEGER :: IS_ELASTIC
        REAL(RDT) :: RELAX

    CONTAINS

    PROCEDURE,PUBLIC :: READvtx

    END TYPE vtxin
    !->+++++++++++++++++++++++++++++++++++++++++++++
    TYPE,PUBLIC :: cfdin

        real(rdt)::A0_B,A1_B,B1_B,A2_B,B2_B
        real(rdt)::AA0,AA1,BB1,AA2,BB2
        real(rdt)::EE0,EE1,FF1,EE2,FF2
        real(rdt)::R,OMIGA,ALF_TPP
        INTEGER::N_REV,N_B,N_F
        real(rdt)::K2,K4,CFL
        INTEGER::ITMAX
        real(rdt)::MU,MACH_TIP,GAMA,DELT_PSI

    CONTAINS

    PROCEDURE,PUBLIC :: READcfd

    END TYPE cfdin
    !->+++++++++++++++++++++++++++++++++++++++++++++
    TYPE,PUBLIC :: matin
        
        INTEGER :: NBPL
        REAL(RDT) :: R,EL1
        REAL(RDT) :: OMG
        INTEGER :: CHO
        REAL(RDT),ALLOCATABLE :: STIFF(:,:,:),MASS(:,:,:)
        REAL(RDT),ALLOCATABLE  :: BTA(:),TAO(:),RL(:)

    CONTAINS

    PROCEDURE,PUBLIC :: READmat

    END TYPE matin
    !->+++++++++++++++++++++++++++++++++++++++++++++
    TYPE,PUBLIC :: DFMin

        INTEGER :: TPE,VDOF,WDOF,FDOF,UDOF,ITN
        REAL(RDT) :: TOL

    CONTAINS

    PROCEDURE,PUBLIC :: READdfm

    END TYPE DFMin
    !->+++++++++++++++++++++++++++++++++++++++++++++
    CONTAINS
    !->+++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE READIPT(THIS,IIN)
    IMPLICIT NONE
    CLASS(IPTIN) :: THIS
    INTEGER,INTENT(IN) :: IIN

    INTEGER :: NumofRotors,RotorID,N
    REAL(8) :: THETA(6)!共轴上下旋翼总距和周期变距
    INTEGER :: NLG,NPTS
    REAL(RDT) :: BTAP
    REAL(RDT) :: LDTIP(3)
    REAL(RDT) :: MU
    REAL(RDT) :: thetapa, delta
    INTEGER :: ISR
    REAL(RDT) :: CTC
    REAL(RDT) :: RHA
    !REAL(RDT) :: TH0,TH1C,TH1S
    INTEGER :: NOD,RHS,NMN,NMNV,NMNW,NMNT,NMNA,NMNF
    INTEGER :: CNDN,HINGES(6),NFPCT,NHFLP,NPTCT,NHPIT,NLLCT,NHLAG
    INTEGER :: SOL
    INTEGER :: nrbd, afsen
    REAL(RDT) :: lockn, cchordrf, alcsrf
    REAL(RDT) :: kfp, kptbar, klgdpr,kfpspr
    REAL(RDT) :: cptbar, clgdpr,cfpspr
    REAL(RDT) :: ggt
    REAL(RDT) :: rdcm, rdck
    INTEGER :: mkim, mdyim
    INTEGER :: mktrim, mkfm

    NAMELIST /INPUT/ NumofRotors,RotorID,THETA,&
        NLG,&
        NPTS,&
        BTAP,&
        LDTIP,&
        MU,&
        thetapa, delta, isr,&
        CTC, &
        RHA,&
        !TH0,TH1C,TH1S,&
        NOD,RHS,NMN,NMNV,NMNW,NMNT,NMNA,NMNF, &
        CNDN,HINGES,NFPCT,NHFLP,NPTCT,NHPIT,NLLCT,NHLAG,&
        SOL,  &
        nrbd, afsen, &
        lockn, cchordrf, alcsrf, &
        kfp, kptbar, klgdpr,kfpspr,&
        cptbar, clgdpr,cfpspr,&
        ggt, &
        rdcm, rdck, &
        mkim, mdyim, &
        mktrim, mkfm

    REWIND(IIN)
    READ(IIN,NML=INPUT)
    THIS%NumofRotors=NumofRotors
    THIS%RotorID=RotorID
    THIS%THETA=THETA
    
    THIS%NLG=NLG
    THIS%NPTS=NPTS
    THIS%BTAP=BTAP
    THIS%LDTIP=LDTIP
    THIS%MU=MU
    THIS%thetapa=thetapa
    THIS%delta=delta
    THIS%isr=isr
    THIS%CTC =CTC
    THIS%RHA=RHA
    
    DO N=1,NumofRotors
        IF(RotorID.EQ.1)THEN!上旋翼或单旋翼
            THIS%TH0=THETA(1)
            THIS%TH1C=THETA(2)
            THIS%TH1S=THETA(3)
        ELSE IF(RotorID.EQ.2)THEN!下旋翼
            THIS%TH0=THETA(4)
            THIS%TH1C=THETA(5)
            THIS%TH1S=THETA(6)
        END IF
    END DO

    !THIS%TH0=TH0
    !THIS%TH1C=TH1C
    !THIS%TH1S=TH1S
    
    THIS%NOD=NOD
    THIS%RHS=RHS
    THIS%NMN=NMN
    THIS%NMNV=NMNV
    THIS%NMNW=NMNW
    THIS%NMNT=NMNT
    THIS%NMNA=NMNA
    THIS%NMNF=NMNF
    THIS%CNDN=CNDN
    THIS%HINGES=HINGES
    THIS%NFPCT=NFPCT
    THIS%NHFLP=NHFLP
    THIS%NPTCT=NPTCT
    THIS%NHPIT=NHPIT
    THIS%NLLCT=NLLCT
    THIS%NHLAG=NHLAG
    THIS%SOL=SOL
    THIS%nrbd=nrbd
    THIS%afsen=afsen
    THIS%lockn=lockn
    THIS%cchordrf=cchordrf
    THIS%alcsrf=alcsrf
    THIS%kfp=kfp
    THIS%kptbar=kptbar
    THIS%klgdpr=klgdpr
    THIS%kfpspr=kfpspr
    THIS%cptbar=cptbar
    THIS%clgdpr=clgdpr
    THIS%cfpspr=cfpspr
    THIS%ggt=GGT
    THIS%rdcm=RDCM
    THIS%rdck=RDCK
    THIS%mkim=mkim
    THIS%mdyim=mdyim
    THIS%mktrim=mktrim
    THIS%mkfm=mkfm

    IF(ISR.EQ.0) THEN!换算为弧度
        THIS%thetapa=THIS%thetapa*PI/180.0
        THIS%delta=THIS%delta*PI/180.0
        THIS%TH0=THIS%TH0*PI/180.0
        THIS%TH1C=THIS%TH1C*PI/180.0
        THIS%TH1S=THIS%TH1S*PI/180.0
    ELSE IF(ISR.NE.1) THEN
        WRITE(*,*) 'ERROR'
    END IF
    !	CALL MPIC%RECVI(THIS%NPTS,1,0)

    END SUBROUTINE
    !->+++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE READEGN(THIS,IIN)
    IMPLICIT NONE
    CLASS(egnin) :: THIS
    INTEGER,INTENT(IN) :: IIN

    INTEGER :: rotan
    INTEGER :: fch
    INTEGER :: soc
    REAL(RDT) :: TFLCTRL=5.0
    NAMELIST /eigen/  rotan,fch,soc,TFLCTRL

    REWIND(IIN)
    READ(IIN,NML=eigen)
    THIS%rotan=rotan
    THIS%fch=fch
    THIS%soc=soc
    THIS%TFLCTRL=TFLCTRL

    END SUBROUTINE
    !->+++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE READfwd(THIS,IIN)
    IMPLICIT NONE
    CLASS(fwdin) :: THIS
    INTEGER,INTENT(IN) :: IIN

    INTEGER natmod, mrksl
    real(rdt) beta, gam
    integer sel, nit
    integer mktfe, npl, ntfe, nidtfe
    REAL(RDT) :: hdg
    INTEGER :: ndgop,nrd

    NAMELIST /forward/ natmod,mrksl,&
        beta, gam, &
        sel, nit, &
        mktfe, npl, ntfe, nidtfe, &
        hdg, &
        ndgop, &
        nrd

    REWIND(IIN)
    READ(IIN,NML=forward)
    THIS%natmod=natmod
    THIS%mrksl=mrksl
    THIS%beta=beta
    THIS%gam=gam
    THIS%sel=sel
    THIS%nit=nit
    THIS%mktfe=mktfe
    THIS%npl=npl
    THIS%ntfe=ntfe
    THIS%nidtfe=nidtfe
    THIS%hdg=hdg
    THIS%ndgop=ndgop
    THIS%nrd=nrd

    END SUBROUTINE
    !->+++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE READifw(THIS,IIN)
    IMPLICIT NONE
    CLASS(ifwin) :: THIS
    INTEGER,INTENT(IN) :: IIN

    INTEGER  UNITYPE,LINTYPE
    REAL(RDT) KHLMDA,KFLMDA
    real(rdt) FXLMDA,FYLMDA

    NAMELIST /INFLOWPARA/     UNITYPE,KHLMDA,KFLMDA,&
        LINTYPE,FXLMDA,FYLMDA

    REWIND(IIN)
    READ(IIN,NML=INFLOWPARA)
    THIS%UNITYPE=UNITYPE
    THIS%KHLMDA=KHLMDA
    THIS%KFLMDA=KFLMDA
    THIS%LINTYPE=LINTYPE
    THIS%FXLMDA=FXLMDA
    THIS%FYLMDA=FYLMDA

    END SUBROUTINE
    !->+++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE READtrm(THIS,IIN)
    IMPLICIT NONE
    CLASS(trmin) :: THIS
    INTEGER,INTENT(IN) :: IIN

    integer :: isinput
    REAL(RDT) :: cw
    real(rdt) :: rlxr, rlxv
    real(rdt) :: fcdfpdka
    real(rdt) :: rcgb(3), racb(3)
    real(rdt) :: fchbtg(3), mthbtg(3)

    NAMELIST /trimipara/ isinput, &
        cw, &
        rlxr,rlxv,&
        fcdfpdka, &
        rcgb, &
        racb, &
        fchbtg, mthbtg

    REWIND(IIN)
    READ(IIN,NML=trimipara)
    THIS%isinput=isinput
    THIS%cw=cw
    THIS%rlxr=rlxr
    THIS%rlxv=rlxv
    THIS%fcdfpdka=fcdfpdka
    THIS%rcgb=rcgb
    THIS%racb=racb
    THIS%fchbtg=fchbtg
    THIS%mthbtg=mthbtg

    END SUBROUTINE
    !->+++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE READSEC(THIS,IIN)
    IMPLICIT NONE
    CLASS(secin) :: THIS
    INTEGER,INTENT(IN) :: IIN

    REAL(RDT) :: MDIV,CLMAX,FSTALL,MSTALL,GSTALL,HSTALL
    REAL(RDT) :: CLF,CMAC,CMS,DEL0,DEL1,DEL2,DCDDM,MCRIT,ACRIT,ALFD,CDF
    integer  :: OPTIP,MACRT
    REAL(RDT)  :: BTIP

    NAMELIST /SECTIONFPARA/    OPTIP,BTIP,MACRT,&
        MDIV,CLMAX,FSTALL,MSTALL,GSTALL,HSTALL,&
        CLF,CMAC,CMS,DEL0,DEL1,DEL2,&
        DCDDM,MCRIT,ACRIT,ALFD,CDF

    REWIND(IIN)
    READ(IIN,NML=SECTIONFPARA)
    THIS%OPTIP=OPTIP
    THIS%BTIP=BTIP
    THIS%MACRT=MACRT
    THIS%MDIV=MDIV
    THIS%CLMAX=CLMAX
    THIS%FSTALL=FSTALL
    THIS%MSTALL=MSTALL
    THIS%GSTALL=GSTALL
    THIS%HSTALL=HSTALL
    THIS%CLF=CLF
    THIS%CMAC=CMAC
    THIS%CMS=CMS
    THIS%DEL0=DEL0
    THIS%DEL1=DEL1
    THIS%DEL2=DEL2
    THIS%DCDDM=DCDDM
    THIS%MCRIT=MCRIT
    THIS%ACRIT=ACRIT
    THIS%ALFD=ALFD
    THIS%CDF=CDF

    END SUBROUTINE
    !->+++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE READvtx(THIS,IIN)
    IMPLICIT NONE
    CLASS(vtxin) :: THIS
    INTEGER,INTENT(IN) :: IIN

    real(rdt) :: ri0
    integer ::nf,nw,nseg,nx
    real(rdt) :: vor_rootr,vor_tipr,rcct
    integer :: IS_ELASTIC
    REAL(RDT) :: RELAX

    NAMELIST /vortex/ 	ri0,nf,nw,nseg,nx,vor_rootr,vor_tipr,rcct,IS_ELASTIC,RELAX

    REWIND(IIN)
    READ(IIN,NML=vortex)
    THIS%ri0=ri0
    THIS%nf=nf
    THIS%nw=nw
    THIS%nseg=nseg
    THIS%nx=nx
    THIS%vor_rootr=vor_rootr
    THIS%vor_tipr=vor_tipr
    THIS%rcct=rcct
    THIS%IS_ELASTIC=IS_ELASTIC
    THIS%RELAX=RELAX

    END SUBROUTINE
    !->+++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE READcfd(THIS,IIN)
    IMPLICIT NONE
    CLASS(cfdin) :: THIS
    INTEGER,INTENT(IN) :: IIN

    real(rdt)::A0_B,A1_B,B1_B,A2_B,B2_B
    real(rdt)::AA0,AA1,BB1,AA2,BB2
    real(rdt)::EE0,EE1,FF1,EE2,FF2
    real(rdt)::R,OMIGA,ALF_TPP
    INTEGER::N_REV,N_B,N_F
    real(rdt)::K2,K4,CFL
    INTEGER::ITMAX
    real(rdt)::MU,MACH_TIP,GAMA,DELT_PSI


    NAMELIST /CFD_PARA/ A0_B,A1_B,B1_B,A2_B,B2_B,AA0,AA1,BB1,AA2,BB2,&
        EE0,EE1,FF1,EE2,FF2,&
        R,OMIGA,ALF_TPP,&
        N_REV,N_B,N_F,&
        K2,K4,CFL,ITMAX,&
        MU,MACH_TIP,GAMA,DELT_PSI

    REWIND(IIN)
    READ(IIN,NML=CFD_PARA)
    THIS%A0_B=A0_B
    THIS%A1_B=A1_B
    THIS%B1_B=B1_B
    THIS%A2_B=A2_B
    THIS%B2_B=B2_B
    THIS%AA0=AA0
    THIS%AA1=AA1
    THIS%BB1=BB1
    THIS%AA2=AA2
    THIS%BB2=BB2
    THIS%EE0=EE0
    THIS%EE1=EE1
    THIS%FF1=FF1
    THIS%EE2=EE2
    THIS%FF2=FF2
    THIS%R=R
    THIS%OMIGA=OMIGA
    THIS%ALF_TPP=ALF_TPP
    THIS%N_REV=N_REV
    THIS%N_B=N_B
    THIS%N_F=N_F
    THIS%K2=K2
    THIS%K4=K4
    THIS%CFL=CFL
    THIS%ITMAX=ITMAX
    THIS%MU=MU
    THIS%MACH_TIP=MACH_TIP
    THIS%GAMA=GAMA
    THIS%DELT_PSI=DELT_PSI

    END SUBROUTINE
    !->+++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE READmat(THIS,IIN)
    IMPLICIT NONE
    CLASS(matin) :: THIS
    INTEGER,INTENT(IN) :: IIN

    INTEGER :: I,J,K,TP,N

    READ(IIN,*) THIS%OMG
    READ(IIN,*) THIS%R
    READ(IIN,*) THIS%EL1
    READ(IIN,*) THIS%NBPL
    
    N=THIS%NBPL
    ALLOCATE(THIS%STIFF(6,6,N),THIS%MASS(6,6,N))
    ALLOCATE(THIS%BTA(N),THIS%TAO(N),THIS%RL(N))
    DO K=1,N
        
        READ(IIN,*) TP
        READ(IIN,*) THIS%CHO
        READ(IIN,*) THIS%RL(K)
        READ(IIN,*) THIS%BTA(K)
        READ(IIN,*) THIS%TAO(K)

        DO I=1,6
            READ(IIN,*) (THIS%MASS(I,J,K),J=1,6)
        END DO

        THIS%STIFF(:,:,K)=0.0D0

        IF(THIS%CHO.EQ.0) THEN
            DO I=1,4
                READ(IIN,*) (THIS%STIFF(I,J,K),J=1,4)
            END DO
        ELSE
            DO I=1,6
                READ(IIN,*) (THIS%STIFF(I,J,K),J=1,6)
            END DO
        END IF
    END DO

    END SUBROUTINE

    !->+++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE READdfm(THIS,IIN)
    IMPLICIT NONE
    CLASS(dfmin) :: THIS
    INTEGER,INTENT(IN) :: IIN


    INTEGER :: TPE,VDOF,WDOF,FDOF,UDOF,ITN
    REAL(RDT) :: TOL

    NAMELIST /DFM_PARA/ TPE,VDOF,WDOF,FDOF,UDOF,ITN,TOL

    REWIND(IIN)
    READ(IIN,NML=DFM_PARA)
    THIS%TPE=TPE
    THIS%VDOF=VDOF
    THIS%WDOF=WDOF
    THIS%FDOF=FDOF
    THIS%UDOF=UDOF
    THIS%ITN=ITN
    THIS%TOL=TOL

    END SUBROUTINE
    !->+++++++++++++++++++++++++++++++++++++++++++++
    END MODULE
