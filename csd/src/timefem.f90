module TFEM_CLASS
USE EIGENBEAM_CLASS
USE AERO_CLASS
IMPLICIT NONE
!->+++++++++++++++++++++++++++++++++++++++++++++
TYPE,PUBLIC:: TFEM

    integer npl, ntfe, nidtfe
    integer ngtfe, nltfe
    REAL(RDT) qnu
    real(rdt) sdt, dltpsid, dltm
    REAL(RDT), ALLOCATABLE:: xtfe(:), xeltfe(:), xtfep(:)
    REAL(RDT), ALLOCATABLE:: xtfex(:), xtfexx(:), xtfepx(:), xtfepxx(:)
    REAL(RDT), ALLOCATABLE:: Kmtg(:, :), Pvtg(:)
    REAL(RDT), ALLOCATABLE:: Kmtl(:, :), Pvtl(:)
    REAL(RDT), ALLOCATABLE:: Hmt(:, :, :), Hdt(:, :, :), Hd2t(:, :, :)
    REAL(RDT), ALLOCATABLE:: Hmtgs(:, :), Hdtgs(:, :), Hd2tgs(:, :), sstfm(:)
    
    CONTAINS
	
    PROCEDURE,PUBLIC :: CONSTRUCT_TFEM
    PROCEDURE,PUBLIC :: tfemini
    PROCEDURE,PUBLIC :: timefmwp
    PROCEDURE,PUBLIC :: timefm
    PROCEDURE,PUBLIC :: kpeletfm
    PROCEDURE,PUBLIC :: MTXTFMWP
    PROCEDURE,PUBLIC :: mtxtfm
    PROCEDURE,PUBLIC :: procuqstfm

END TYPE TFEM
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++PUBLIC++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
CONTAINS
SUBROUTINE CONSTRUCT_TFEM(THIS,INPU)    
implicit none
CLASS(TFEM) :: THIS
TYPE(INPUT) :: INPU

	THIS%npl=INPU%FWD%npl
	THIS%ntfe=INPU%FWD%ntfe
	THIS%nidtfe=INPU%FWD%nidtfe

END SUBROUTINE
!->+++++++++++++++++++++++++++++++++++++++++++++
subroutine tfemini(this,nhf,nstp,omg,q,yst)
USE Mmyfunifs
implicit none
class(tfem) :: this
integer,intent(in) :: nhf,nstp
real(rdt),intent(in) :: omg
real(rdt),intent(inout) :: q(nhf),yst(2*nhf)
integer i, ii, igs
integer :: j,jj,err
real(rdt) idmt(nhf, nhf)
real(rdt), parameter:: oned=1.0D0

!real(rdt),external :: lagbaswp,lagbasdtwp2,lagbasd2twp2


    this%ngtfe=nstp*nhf
    this%nltfe=(this%npl+1)*nhf


     ALLOCATE(this%xtfe(this%ngtfe), this%xeltfe(this%nltfe), this%xtfep(this%ngtfe), &
              this%xtfex(this%ngtfe), this%xtfexx(this%ngtfe), this%xtfepx(this%ngtfe),&
              this%xtfepxx(this%ngtfe), this%Kmtg(this%ngtfe, this%ngtfe), &
              this%Pvtg(this%ngtfe), this%Kmtl(this%nltfe, this%nltfe), this%Pvtl(this%nltfe), &
              this%Hmt(nhf, this%nltfe, nlg), this%Hdt(nhf, this%nltfe, nlg), &
              this%Hd2t(nhf, this%nltfe, nlg), this%Hmtgs(this%npl+1, nlg), &
              this%Hdtgs(this%npl+1, nlg), this%Hd2tgs(this%npl+1, nlg), this%sstfm(nlg), &
              STAT=ERR)  
	IF(ERR.NE.0) THEN
		WRITE(*, 998)
		STOP
	END IF       

        q=yst(1:nhf)  

        jj=0
        do 550 j=1, this%npl*this%ntfe+1
            this%xtfep(jj+1:jj+nhf)=q
            jj=jj+nhf
550     continue
        !>   xtfep=0.0D0
        this%xtfe=this%xtfep
                  
    this%dltpsid=360.0D0/this%ntfe 
    this%dltm=this%dltpsid*ACOS(-oned)/180.0D0/omg
    this%sdt=oned/this%dltm
	
      do 10 i=1, nlg
            this%sstfm(i)=0.5d0*(GT(i, nlg)+1.0D0)
     
            this%Hmtgs(:, i)=lagbaswp(this%npl, this%sstfm(i))
            this%Hdtgs(:, i)=lagbasdtwp2(this%npl, this%sstfm(i), this%sdt) 
            !-> Hdtgs(:, i)=lagbasdtwp(npl, ss(i), sdt)  
            this%Hd2tgs(:, i)=lagbasd2twp2(this%npl, this%sstfm(i), this%sdt) 
10  continue

      call eyes(idmt, nhf)
      
      do 20 igs=1, nlg
            do 30 i=1, this%npl+1
                  ii=(i-1)*nhf 
                  this%Hmt(1:nhf,  ii+1:ii+nhf, igs)=this%Hmtgs(i, igs)*idmt
                  this%Hdt(1:nhf,  ii+1:ii+nhf, igs)=this%Hdtgs(i, igs)*idmt
                  this%Hd2t(1:nhf,  ii+1:ii+nhf, igs)=this%Hd2tgs(i, igs)*idmt
30        continue   
20  continue   
	
	
998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
end subroutine 
!->+++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE timefmwp(THIS,BLADE,AEROF,TH0,TH1C,TH1S,SOL,omg)
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
implicit none
CLASS(TFEM) :: THIS
CLASS(EIGENBEAM) :: BLADE
CLASS(AERO) :: AEROF

INTEGER,INTENT(IN) :: SOL
REAL(RDT),INTENT(IN) :: TH0,TH1C,TH1S

REAL(RDT),INTENT(IN) :: OMG
REAL(RDT), EXTERNAL:: DNRM2
integer ii
REAL(RDT) :: ONE=1.0
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> main body:	

    do 210 ii=1, THIS%nidtfe
        call THIS%timefm(BLADE,AEROF,TH0,TH1C,TH1S,OMG,SOL,omg)
        THIS%qnu=DNRM2(THIS%ngtfe, THIS%xtfe-THIS%xtfep, 1)/SQRT(THIS%ngtfe*ONE)            
        THIS%xtfep=THIS%xtfe
        write(*,*) ii,THIS%qnu
210 continue

998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
end subroutine timefmwp
!->+++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE timefm(THIS,BLADE,AEROF,TH0,TH1C,TH1S,OMG,SOL,omga)
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
implicit none
CLASS(TFEM) :: THIS
CLASS(EIGENBEAM) :: BLADE
CLASS(AERO) :: AEROF

INTEGER,INTENT(IN) :: SOL
REAL(RDT),INTENT(IN) :: TH0,TH1C,TH1S,OMG

real(RDT), intent(in):: omga

character(80) filenam
integer i, ii
integer, parameter:: ione=1
real(rdt), parameter:: oned=1.0D0

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> main body:	

      THIS%Kmtg=0.0D0
      THIS%Pvtg=0.0D0 
      ii=0 
      do 50 i=1, THIS%ntfe
            THIS%xeltfe(:)=THIS%xtfe(ii+1:ii+THIS%nltfe)
            call THIS%kpeletfm(BLADE,AEROF,TH0,TH1C,TH1S,OMG,SOL,I)
            THIS%Kmtg(ii+1:ii+THIS%nltfe, ii+1:ii+THIS%nltfe)=THIS%Kmtg(ii+1:ii+THIS%nltfe, ii+1:ii+THIS%nltfe)+THIS%Kmtl(:, :)
            THIS%Pvtg(ii+1:ii+THIS%nltfe)=THIS%Pvtg(ii+1:ii+THIS%nltfe)+THIS%Pvtl(:)
            ii=ii+THIS%npl*BLADE%nhf 
50  continue


      call lincon2(THIS%xtfe, THIS%ngtfe, BLADE%NHF, THIS%Kmtg, THIS%Pvtg)
!      call linsolver2(xtfe, 1, ngtfe, Kmtg, Pvtg)
!      xtfe(1:nhf)=xtfe(ngtfe-nhf+1:ngtfe)
      
      !filenam='xxreds0'
      !call outmatx(filenam, ngtfe, 1, Pvtg)   
      
      !filenam='kkk'
      !call outmatx(filenam, ngtfe, ngtfe, Kmtg)   
      
!      filenam='xxtfe'
!      call outmatx(filenam, ngtfe, 1, xtfe)   
!            
!      call DGEMV('N', ngtfe, ngtfe, -oned, Kmtg, ngtfe, xtfe, ione, oned, Pvtg, ione)  
!      filenam='xxreds'
!      call outmatx(filenam, ngtfe, 1, Pvtg)      
       
end subroutine timefm
!->+++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE kpeletfm(THIS,BLADE,AEROF,TH0,TH1C,TH1S,OMG,SOL,itfe)
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
implicit none
CLASS(TFEM) :: THIS
CLASS(EIGENBEAM) :: BLADE
CLASS(AERO) :: AEROF

INTEGER,INTENT(IN) :: SOL
REAL(RDT),INTENT(IN) :: TH0,TH1C,TH1S,OMG
integer, intent(in):: itfe

integer i
integer err
real(rdt), parameter:: oned=1.0D0
real(rdt), allocatable:: Kmtlt(:, :), Pvtlt(:)

	ALLOCATE(Kmtlt(THIS%nltfe, THIS%nltfe), Pvtlt(THIS%nltfe), &
            STAT=ERR)

	IF(ERR.NE.0) THEN
		WRITE(*, 998)
		STOP
	END IF
	
      THIS%Kmtl=0.0D0
      THIS%Pvtl=0.0D0
      
      DO 10 i=1, NLG
            call THIS%mtxtfmwp(Kmtlt, Pvtlt, BLADE,AEROF,TH0,TH1C,TH1S,OMG,SOL,&
                            itfe, this%sstfm(i), THIS%Hmt(:, :, i), THIS%Hdt(:, :, i), THIS%Hd2t(:, :, i))
	THIS%Kmtl=THIS%Kmtl+Kmtlt*GC(I,NLG)
	THIS%Pvtl=THIS%Pvtl+Pvtlt*GC(I,NLG)
		
10	CONTINUE

	THIS%Kmtl=THIS%Kmtl*THIS%dltm/2.0D0
	THIS%Pvtl=THIS%Pvtl*THIS%dltm/2.0D0


	DEALLOCATE(Kmtlt, Pvtlt, &
            STAT=ERR)

	IF(ERR.NE.0) THEN
		WRITE(*, 999)
		STOP
	END IF     
	 
998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')	
end subroutine kpeletfm
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> the driver / wrapper.
!-> get the matrices used by time finite method.
!->       H[T]*K*H+H*C*Hdt-Hdt[T]*M*Hdt, H[T]*F.
SUBROUTINE mtxtfmwp(THIS,Kmtl, Pvtl, BLADE,AEROF,TH0,TH1C,TH1S,OMG,SOL,itfe, ss, Hmt, Hdt, Hd2t)
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
implicit none
CLASS(TFEM) :: THIS
CLASS(EIGENBEAM) :: BLADE
CLASS(AERO) :: AEROF

INTEGER,INTENT(IN) :: SOL
REAL(RDT),INTENT(IN) :: TH0,TH1C,TH1S,OMG

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
INTEGER, INTENT(IN)::  itfe
real(RDT), intent(in):: ss
REAL(RDT), INTENT(in):: Hmt(BLADE%nhf, THIS%nltfe), Hdt(BLADE%nhf, THIS%nltfe), Hd2t(BLADE%nhf, THIS%nltfe)
REAL(RDT), INTENT(OUT):: Kmtl(THIS%nltfe, THIS%nltfe), Pvtl(THIS%nltfe)
REAL(RDT) :: PPSI
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->       mtxtfm
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
integer i, ii
real(rdt), parameter:: oned=1.0D0
real(rdt) psid

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

      psid=(itfe-1)*THIS%dltpsid+ss*THIS%dltpsid
      ppsi=psid*ACOS(-ONED)/180.0D0 		
      
      CALL MVMU2(BLADE%q, BLADE%nhf, THIS%nltfe, Hmt, THIS%xeltfe)
      CALL MVMU2(BLADE%q1t, BLADE%nhf, THIS%nltfe, Hdt, THIS%xeltfe)
      CALL MVMU2(BLADE%q2t, BLADE%nhf, THIS%nltfe, Hd2t, THIS%xeltfe)      
      
      if( BLADE%natmod .eq. 1 ) then 
            goto 131
      end if 
   
    BLADE%UQS=BLADE%q
    BLADE%UQS1T=BLADE%q1t
    BLADE%UQS2T=BLADE%q2t


    CALL BLADE%UPDATE()
    CALL AEROF%MLEGGAUSS_AERO(OMG,BLADE,TH0,TH1C,TH1S,PPSI,SOL)
    CALL BLADE%MLEGGAUSS_ST_BEAM(TH0,TH1C,TH1S,PPSI,OMG,SOL)
    CALL BLADE%FORMULATE_ST()

      BLADE%ucmc=BLADE%ucmc+&
                                        BLADE%rdcm*BLADE%ummc+&
                                        BLADE%rdck*BLADE%ukmc
    
      BLADE%mmt=BLADE%UMMC
      BLADE%cmt=BLADE%UCMC
      BLADE%kmt=BLADE%UKMC
      BLADE%fv=-BLADE%UFVC
      goto 200
      
131 continue

      CALL MVMU2(BLADE%UQS, BLADE%NDOFC, BLADE%NMN, BLADE%VTEG, BLADE%q)
      CALL MVMU2(BLADE%UQS1T, BLADE%NDOFC, BLADE%NMN, BLADE%VTEG, BLADE%q1t)
      CALL MVMU2(BLADE%UQS2T, BLADE%NDOFC, BLADE%NMN, BLADE%VTEG, BLADE%q2t)

    CALL BLADE%UPDATE()
    CALL AEROF%MLEGGAUSS_AERO(OMG,BLADE,TH0,TH1C,TH1S,PPSI,SOL)
    CALL BLADE%MLEGGAUSS_ST_BEAM(TH0,TH1C,TH1S,PPSI,OMG,SOL)
    CALL BLADE%FORMULATE_ST()
      
      BLADE%ucmc=BLADE%ucmc+&
                                        BLADE%rdcm*BLADE%ummc+&
                                        BLADE%rdck*BLADE%ukmc
            
       CALL BLADE%GETMKCFT()
      goto 200

      
200 continue
    
      call this%mtxtfm(Kmtl, Pvtl, BLADE,Hmt, Hdt)

end subroutine mtxtfmwp
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:jyuhjh
!-> get the matrices used by time finite method.
!->       H[T]*K*H+H*C*Hdt-Hdt[T]*M*Hdt, H[T]*F.
SUBROUTINE mtxtfm(THIS,Kmtl, Pvtl, BLADE, Hmt, Hdt)
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
implicit none
CLASS(TFEM) :: THIS
CLASS(EIGENBEAM) :: BLADE
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
REAL(RDT), INTENT(in):: Hmt(BLADE%NHF, this%nltfe), Hdt(BLADE%NHF, this%nltfe)
REAL(RDT), INTENT(OUT):: Kmtl(this%nltfe, this%nltfe), Pvtl(this%nltfe)
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->       MATMLYTP, MVMU2  
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
integer err
real(rdt), allocatable:: mtt(:, :), HmtT(:, :), HdtT(:, :)

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

    ALLOCATE(mtt(this%nltfe, this%nltfe), HmtT(this%nltfe, BLADE%NHF), HdtT(this%nltfe, BLADE%NHF), &
        STAT=ERR)

    IF(ERR.NE.0) THEN
        WRITE(*, 998)
        STOP
    END IF


    HmtT=transpose(Hmt)
    HdtT=transpose(Hdt)

    call MATMLYTP(mtt, this%nltfe, BLADE%NHF, BLADE%NHF, this%nltfe, HmtT, BLADE%Kmt, Hmt) 
    Kmtl=mtt
    call MATMLYTP(mtt, this%nltfe, BLADE%NHF, BLADE%NHF, this%nltfe, HmtT, BLADE%Cmt, Hdt) 
    Kmtl=Kmtl+mtt
    call MATMLYTP(mtt, this%nltfe, BLADE%NHF, BLADE%NHF, this%nltfe, HdtT, BLADE%Mmt, Hdt) 
    Kmtl=Kmtl-mtt

    call MVMU2(Pvtl, this%nltfe, BLADE%NHF, HmtT, BLADE%Fv)     


    DEALLOCATE(mtt, HmtT, HdtT, &
        STAT=ERR)

    IF(ERR.NE.0) THEN
        WRITE(*, 999)
        STOP
    END IF     
	 
998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
end subroutine mtxtfm
!->+++++++++++++++++++++++++++++++++++++++++++++

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> to output the response.
subroutine procuqstfm(THIS,BLADE,UQSA, UQS1TA,UQS2TA,nstp) 
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
USE Mmyfunifs
implicit none
CLASS(TFEM) :: THIS
CLASS(EIGENBEAM) :: BLADE
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
integer, intent(in) :: NSTP
REAL(RDT),INTENT(OUT) :: UQSA(BLADE%NDOFC,NSTP),UQS1TA(BLADE%NDOFC,NSTP),UQS2TA(BLADE%NDOFC,NSTP)
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->       lagbasdtwp2, lagbasdtwp, lagbasd2twp2
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->       MVMU2, outmatx
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
integer i, ii, j, jj, kk, err
real(rdt),ALLOCATABLE :: hhx(:,:),hhtx(:, :),hh2tx(:, :)
real(rdt) ss
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> main body:



    ALLOCATE(hhx(this%npl+1, this%npl),hhtx(this%npl+1, this%npl),hh2tx(this%npl+1, this%npl),&
        stat=err)
	    IF(ERR.NE.0) THEN
		    WRITE(*, 998)
		    STOP
	    END IF
        do 10 i=1, this%npl
            ss=1.0D0/this%npl*(i-1)
            hhx(:, i)=lagbaswp(this%npl, ss)
            hhtx(:, i)=lagbasdtwp2(this%npl, ss, THIS%sdt) 
            hh2tx(:, i)=lagbasd2twp2(this%npl, ss, THIS%sdt) 
10  continue


      ii=0
      do 550 i=1, THIS%ntfe
            THIS%xeltfe(:)=THIS%xtfe(ii+1:ii+THIS%nltfe)
            ii=ii+THIS%npl*BLADE%nhf 

            do 540 j=1, THIS%npl            
                  jj=0
                  BLADE%q=0.0D0
                  BLADE%q1t=0.0D0
                  BLADE%q2t=0.0D0
                  do 530 kk=1, THIS%npl+1
                        BLADE%q=BLADE%q+THIS%xeltfe(jj+1:jj+BLADE%nhf)*hhx(kk, j)
                        BLADE%q1t=BLADE%q1t+THIS%xeltfe(jj+1:jj+BLADE%nhf)*hhtx(kk, j)
                        BLADE%q2t=BLADE%q2t+THIS%xeltfe(jj+1:jj+BLADE%nhf)*hh2tx(kk, j)
                        jj=jj+BLADE%nhf
530             continue

                  if( BLADE%natmod .eq. 1 ) then
                          CALL MVMU2(BLADE%UQS, BLADE%NDOFC, BLADE%NMN, BLADE%VTEG, BLADE%q)
                          CALL MVMU2(BLADE%UQS1T, BLADE%NDOFC, BLADE%NMN, BLADE%VTEG, BLADE%q1t)
                          CALL MVMU2(BLADE%UQS2T, BLADE%NDOFC, BLADE%NMN, BLADE%VTEG, BLADE%q2t)
                  else if(BLADE%natmod .eq. 0 ) then
                        BLADE%UQS=BLADE%q
                        BLADE%UQS1T=BLADE%q1t
                        BLADE%UQS2T=BLADE%q2t
                  end if
                 
                  UQSA(:, (i-1)*THIS%NPL+j)=BLADE%uqs
                  UQS1TA(:, (i-1)*THIS%NPL+j)=BLADE%uqs1T
                  UQS2TA(:, (i-1)*THIS%NPL+j)=BLADE%uqs2T   
540       continue
550 continue

                  UQSA(:, nstp)=UQSA(:, 1)
                  UQS1TA(:, nstp)=UQS1TA(:, 1)
                  UQS2TA(:, nstp)=UQS2TA(:, 1)  

      
    deallocate(hhx,hhtx,hh2tx,&
    stat=err)
	IF(ERR.NE.0) THEN
		WRITE(*, 999)
		STOP
	END IF


998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
end subroutine procuqstfm
!->+++++++++++++++++++++++++++++++++++++++++++++


END MODULE 
