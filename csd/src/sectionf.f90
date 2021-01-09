MODULE SECTIONF_CLASS
USE INPUT_CLASS
IMPLICIT NONE
!->+++++++++++++++++++++++++++++++++++++++++++++
TYPE,PUBLIC :: SECTIONF

    REAL(RDT) :: MDIV,CLMAX,FSTALL,MSTALL,GSTALL,HSTALL
    REAL(RDT) :: CLF,CMAC,CMS,DEL0,DEL1,DEL2,DCDDM,MCRIT,ACRIT,ALFD,CDF
    integer OPTIP,MACRT
    REAL(RDT) BTIP
    
    real(rdt) :: afas,cls
    
    CONTAINS

    PROCEDURE,PUBLIC ::  CONSTRUCT_SECTIONF
    PROCEDURE,PUBLIC :: GETCL
    PROCEDURE,PUBLIC ::  GETCD
    PROCEDURE,PUBLIC :: TIPLOSS
    PROCEDURE,PUBLIC :: getstangle
    PROCEDURE,PUBLIC :: getstcl
    PROCEDURE,PUBLIC :: getstcm
    PROCEDURE,PUBLIC :: getstcd
    PROCEDURE,PUBLIC :: REVERSECORRECT
    PROCEDURE,PUBLIC :: testsc
    PROCEDURE,PUBLIC :: AEROST2


    PROCEDURE,PRIVATE :: prandtiplossf1
    PROCEDURE,PRIVATE :: prandtiplossf

!    PROCEDURE,PUBLIC ::linearinflow
!
!    PROCEDURE,PRIVATE :: VEQ32
!    PROCEDURE,PRIVATE :: jonsoninflow


END TYPE SECTIONF
!->+++++++++++++++++++++++++++++++++++++++++++++
CONTAINS
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++PUBLIC++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE CONSTRUCT_SECTIONF(THIS,INPU)
IMPLICIT NONE
CLASS(SECTIONF) :: THIS
CLASS(INPUT) :: INPU

    THIS%OPTIP=INPU%SEC%OPTIP
    THIS%BTIP=INPU%SEC%BTIP
    THIS%MACRT=INPU%SEC%MACRT
    THIS%MDIV=INPU%SEC%MDIV
    THIS%CLMAX=INPU%SEC%CLMAX
    THIS%FSTALL=INPU%SEC%FSTALL
    THIS%MSTALL=INPU%SEC%MSTALL
    THIS%GSTALL=INPU%SEC%GSTALL
    THIS%HSTALL=INPU%SEC%HSTALL
    THIS%CLF=INPU%SEC%CLF
    THIS%CMAC=INPU%SEC%CMAC
    THIS%CMS=INPU%SEC%CMS
    THIS%DEL0=INPU%SEC%DEL0
    THIS%DEL1=INPU%SEC%DEL1
    THIS%DEL2=INPU%SEC%DEL2
    THIS%DCDDM=INPU%SEC%DCDDM
    THIS%MCRIT=INPU%SEC%MCRIT
    THIS%ACRIT=INPU%SEC%ACRIT
    THIS%ALFD=INPU%SEC%ALFD
    THIS%CDF=INPU%SEC%CDF
    THIS%ALFD=THIS%ALFD/57.3

END SUBROUTINE
    
SUBROUTINE GETCL(THIS,CLA,MA,AFA)
IMPLICIT NONE
CLASS(SECTIONF) :: THIS

REAL(RDT),INTENT(OUT) :: CLA
REAL(RDT),INTENT(IN) :: MA,AFA
    
REAL(RDT) :: COE1
    
    IF(THIS%MACRT.EQ.0) 	RETURN

    IF(MA.LE.THIS%MDIV) THEN
        COE1=1.0/SQRT(1-MA**2)
    ELSE IF(MA.GT.THIS%MDIV .AND. MA.LE.(THIS%MDIV+0.1)) THEN
        COE1=(1-MA)/((1-THIS%MDIV)*SQRT(1-THIS%MDIV**2))
    ELSE
        COE1=(1-MA)/((1-THIS%MDIV)*SQRT(1-THIS%MDIV**2))+(MA-THIS%MDIV-0.1)/(1-THIS%MDIV-0.1)
    END IF
    CLA=CLA*COE1
        
END SUBROUTINE
        
SUBROUTINE GETCD(THIS,CD,MA,AFA)
IMPLICIT NONE
CLASS(SECTIONF) :: THIS

REAL(RDT),INTENT(OUT) :: CD
REAL(RDT),INTENT(IN) :: MA,AFA
    
REAL(RDT) :: MC,DELTACD
    
    IF(THIS%MACRT.EQ.0) 	RETURN
    MC=MAX(0.0D0,THIS%MCRIT*(1-ABS(AFA)/THIS%ACRIT))
    DELTACD=MAX(0.0D0,THIS%DCDDM*(MA-MC))
    CD=THIS%DEL0+THIS%DEL1*AFA+THIS%DEL2*AFA**2+DELTACD
        
END SUBROUTINE
    
    
SUBROUTINE  TIPLOSS(THIS,fbg, nb, rb, ladi) 
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
IMPLICIT NONE
CLASS(SECTIONF) :: THIS

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
integer,intent(in) :: nb
real(rdt),intent(in) :: rb, ladi
real(rdt),intent(out) :: fbg

real(rdt) :: b
   IF(THIS%OPTIP.EQ.0) THEN
        FBG=1.0
    ELSE 
        if(ABS(THIS%BTIP-0.0).gt.1e-3) then
            b=THIS%BTIP
        else
            b=1-ladi/dble(nb)*2*log(2.0)
        end if
        if(THIS%OPTIP.eq.1) then
            if(rb.lt.b) then
                fbg=0
            else
                fbg=1
            end if
        else if(THIS%OPTIP.eq.2) then
	        call THIS%prandtiplossf(fbg,nb, rb, ladi,b)  
	    end if
    END IF
    
end subroutine

  


!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> Prandtl's tip-loss function.
subroutine prandtiplossf1(THIS,fbg, nb, rb, ladi) 
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
IMPLICIT NONE
CLASS(SECTIONF) :: THIS

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
integer nb
real(rdt) rb, ladi
real(rdt) fbg

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
real(rdt), external:: paai
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->        
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
real(rdt) f

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

      f=0.5d0*dble(1)*(1.0D0-rb)/ladi
      fbg=2.0D0*ACOS(EXP(-f))/paai() 
end subroutine prandtiplossf1
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> Prandtl's tip-loss function.
subroutine prandtiplossf(THIS,fbg, nb, rb, ladi,b) 
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
IMPLICIT NONE
CLASS(SECTIONF) :: THIS

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
integer nb
real(rdt),intent(in) :: rb, ladi,b
real(rdt),intent(out) :: fbg

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
real(rdt), external:: paai
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->        
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
real(rdt) f

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

      f=0.5d0**((1.0D0-rb)/(1.0D0-b))
      fbg=2.0D0*ACOS(f)/paai() 
      
end subroutine prandtiplossf
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


subroutine getstangle(THIS,ma,cl)
IMPLICIT NONE
CLASS(SECTIONF) :: THIS
real(rdt),intent(in) :: ma,cl

    THIS%cls=THIS%CLMAX*min(1.0,((1-ma)+THIS%FSTALL*(ma-THIS%MSTALL))/(1-THIS%MSTALL))
    THIS%afas=THIS%cls/cl
    
end subroutine


subroutine getstcl(THIS,cl,afa,ma)
IMPLICIT NONE
CLASS(SECTIONF) :: THIS
real(rdt),intent(in) :: ma,afa
real(rdt),intent(inout) :: cl


real(rdt) :: tx1,tx2,tx3,tx4

    
    if(ABS(afa).gt.45/57.3) then
        cl=THIS%CLF*SIN(2*AFA)
    ELSE
        TX1=THIS%HSTALL*THIS%cls
        TX2=THIS%CLF*SIN(2*ABS(AFA))
        TX3=MAX(TX1,TX2)
        
        TX1=(THIS%GSTALL*THIS%afas-ABS(AFA))*THIS%cls+(ABS(AFA)-THIS%afas)*THIS%HSTALL*THIS%cls
        TX2=TX1/(THIS%GSTALL*THIS%afas-THIS%afas)
        TX4=MAX(TX3,TX2)
        CL=SIGN(TX4,AFA)
    END IF
    
end subroutine

subroutine getstcm(THIS,cm,afa,ma)
IMPLICIT NONE
CLASS(SECTIONF) :: THIS
real(rdt),intent(in) :: ma,afa
real(rdt),intent(inout) :: cm


real(rdt) :: CDD,cmf,tx1,tx2,tx3,tx4

    CALL THIS%GETCD(CDD,MA,THIS%ALFD)
    CMF=-0.25*(CDD+THIS%CDF)
    if(ABS(afa).LT.60/57.3) then
        
        TX1=(60.0/57.3-ABS(afa))*THIS%CMS
        TX2=(ABS(AFA)-THIS%afas)*0.75*CMF
        TX3=(TX1+TX2)/(60/57.3-THIS%afas)
        
        CM=SIGN(TX3,AFA)
    ELSE
    
        TX1=(90.0/57.3-ABS(afa))*0.75*CMF
        TX2=(ABS(AFA)-60/57.3)*CMF
        TX3=(TX1+TX2)/(30/57.3)
        
        CM=SIGN(TX3,AFA)
    END IF
    
end subroutine

subroutine getstcd(THIS,cd,afa,ma)
IMPLICIT NONE
CLASS(SECTIONF) :: THIS
real(rdt),intent(in) :: ma,afa
real(rdt),intent(inout) :: cd


real(rdt) :: CDD,tx1

    CALL THIS%GETCD(CDD,MA,THIS%ALFD)
   
   tx1=(ABS(afa)-THIS%ALFD)/(90/57.3-THIS%ALFD)*90/57.3
    cd=CDD+THIS%CDF*SIN(tx1)
    
end subroutine


SUBROUTINE  REVERSECORRECT(THIS,CL,CD,CM,AFA)
IMPLICIT NONE
CLASS(SECTIONF) :: THIS
REAL(RDT),INTENT(IN) :: CL,CD,AFA
REAL(RDT),INTENT(OUT) :: CM

REAL(RDT) AFAE

    IF(ABS(AFA).GT.3.1415926/2) THEN
        AFAE=AFA-SIGN(3.14159,AFA)
        CM=CM+0.5*COS(AFAE)*CL+0.5*SIN(AFAE)*CD
    END IF
    
END SUBROUTINE





subroutine testsc(THIS)
IMPLICIT NONE
CLASS(SECTIONF) :: THIS

real(rdt) afa,cl,cm,cd
real(rdt) lbma,alct,CD0T
integer i

!    call THIS%CONSTRUCT_SECTIONF()
    
    OPEN(11,FILE='../test/section.dat')
    lbma=0.4
    alct=5.7
    cd0t=0.01
    
    do i=1,1000
    
        afa=-3.14+6.28/1000*i
    
     CALL THIS%GETCL(ALCT,LBMA,AFA)
     CALL THIS%GETCD(CD0T,LBMA,AFA)
     
    CALL THIS%getstangle(LBMA,ALCT)
	    
	    IF(ABS(AFA).LE.THIS%afas) THEN
	        CL=ALCT*AFA
	        CM=0.0
	    ELSE
	        CALL THIS%getstcl(CL,afa,LBMA)
	        CALL THIS%getstcm(CM,afa,LBMA)
		END IF
		
	    IF(ABS(AFA).LE.THIS%ALFD) THEN
	        CD=CD0T
	    ELSE
	        CALL THIS%getstcd(CD,afa,LBMA)   
	    END IF
	    
	    CALL THIS%REVERSECORRECT(CL,CD,CM,AFA)
	    
        write(11,997) afa*57.3,cl,cd,cm
    end do
    close(11)
 997 FORMAT(10F26.12)
end subroutine

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
SUBROUTINE AEROST2(THIS,CLFX,CMFX,CDFX,ALCT,CD0T,AFA,LBMA)
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
implicit none
CLASS(SECTIONF) :: THIS
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
real(rdt),intent(in) :: AFA,LBMA,ALCT,CD0T
real(rdt),intent(out) :: CLFX,CMFX,CDFX

REAL(RDT) CM,CL,CD

        CALL THIS%getstangle(LBMA,ALCT)
        
	    IF(ABS(AFA).LE.THIS%ALFD) THEN
	        CDFX=CD0T
	    ELSE
	        CALL THIS%getstcd(CD,afa,LBMA)   
	        CDFX=CD
	    END IF
	    
	    IF(ABS(AFA).LE. THIS%AFAS) THEN
	        CLFX=ALCT
	        CMFX= THIS%CMAC
	    ELSE
	        CALL THIS%getstcl(CL,afa,LBMA)
	        CALL THIS%getstcm(CM,afa,LBMA)
!	        CALL REVERSECORRECT(CL,CD,CM,AFA)
	        CLFX=CL/AFA
	        CMFX=CM
	    END IF
		
END SUBROUTINE

END MODULE