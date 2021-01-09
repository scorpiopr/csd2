
!* THE ROOT-MEAN-SQUARE OF A VECTOR.
FUNCTION rms(DM,V)
use GlobalDataFun
implicit none

integer, intent(in)::dm
real(RDT), intent(in)::v(dm)
real(RDT) rms

real(RDT),external::vip

	rms=dsqrt(vip(dm,v,v)/dble(dm))

END 



FUNCTION PHES75(NBPL,UQDC,SMD,UVC2,NAD,NBD,RR,RL,BTA,UQS,MHIV2)
use GlobalDataFun
implicit none

INTEGER, INTENT(IN):: NBPL, UQDC, SMD, UVC2, NAD, NBD 
REAL(RDT), INTENT(IN):: RR
REAL(RDT), INTENT(IN):: RL(NBPL), BTA(NBPL), UQS(UQDC)
INTEGER, INTENT(IN):: MHIV2(SMD)
REAL(RDT) PHES75

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->       
real(RDT),external:: FWM, VIP
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->       VV
!->+++++++++++++++++++++++++++++++++++++++++++++

INTEGER UVC, UVCC
INTEGER NBPLI, NBPLIP1
INTEGER NBOOD
REAL(RDT) XQ(SMD), QS(SMD)
REAL(RDT) QSPH(NBD), TPQ(NBD)
REAL(RDT) SS, LL, BTAT
REAL(RDT) SPH

      UVC=0  
      DO 25 NBPLI=NBPL-1,1,-1
            NBPLIP1=NBPLI+1 
            IF( RL(NBPLI) .LE. 0.75*RR ) THEN
                XQ(1:SMD)=UQS(UVC+1:UVC+SMD)
                QS(MHIV2(1:SMD))=XQ(1:SMD)
                LL=RL(NBPLIP1)-RL(NBPLI)
                SS=(0.75*RR-RL(NBPLI))/LL
                EXIT 
            END IF                  
            UVC=UVC+UVC2
25  CONTINUE

    BTAT=FWM(BTA(NBPLI),BTA(NBPLIP1),SS)

    QSPH(:)=QS(2*NAD+1:2*NAD+NBD)      
    NBOOD=0 
    CALL VV(TPQ,NBOOD,NBD,SS,LL)    
    SPH=VIP(NBD,TPQ,QSPH)

    PHES75=BTAT+SPH
END FUNCTION PHES75



function triseries(a0, a1c, a1s, alp)
use GlobalDataFun
implicit none

REAL(RDT), INTENT(IN):: a0, a1c, a1s, alp
REAL(RDT) triseries

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->       
!->+++++++++++++++++++++++++++++++++++++++++++++

     triseries=a0+a1c*COS(alp)+a1s*SIN(alp)
end function triseries



!* THE INNER PRODUCT OF THE TWO VECTORS.
!* CAN ALSO USE THE INTRINSIC FUNCTION <DOT_PRODUCT> TO GET THE RESULT.
FUNCTION VIP(DM,V1,V2)
use GlobalDataFun
implicit none

INTEGER, INTENT(IN):: DM
REAL(RDT), INTENT(IN):: V1(DM),V2(DM)
REAL(RDT) VIP

INTEGER I

      VIP=0.0
    DO 5 I=1,DM
        VIP=VIP+V1(I)*V2(I)
5   CONTINUE

END FUNCTION VIP



!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> BELOW SUBROUTINES 



!->+++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> get the identity matrix.
SUBROUTINE eyes(mtx, dm)
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
INTEGER, INTENT(IN):: dm
REAL(RDT), INTENT(OUT):: mtx(dm, dm)

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->     
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
integer i

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

      mtx=0.0D0
      do i=1, dm
            mtx(i, i)=1.0D0
      end do

end subroutine eyes
!->+++++++++++++++++++++++++++++++++++++++++++++



! TO OBTAIN THE MATRIXES OF TED AND ITS TIME DERIVATIVE.
SUBROUTINE GETTDE(TDE,TDED1T,BTA,SPH,SPHD1T,SVD1X,SWD1X,SVD1XD1T,SWD1XD1T)
use GlobalDataFun
implicit none

INTEGER, PARAMETER:: DM=3
REAL(RDT),INTENT(IN)::BTA,SPH,SPHD1T,SVD1X,SWD1X,SVD1XD1T,SWD1XD1T
REAL(RDT),INTENT(OUT)::TDE(DM,DM),TDED1T(DM,DM)

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->       
REAL(RDT),EXTERNAL:: FABSC
!->+++++++++++++++++++++++++++++++++++++++++++++



REAL(RDT) TT,CST,SNT,CSB,SNB,TACP,TACPD1T
REAL(RDT) X1, X2, X1D1T, X2D1T
REAL(RDT) Y1, Y2, Y1D1T, Y2D1T

    TT=BTA+SPH
    CST=COS(TT)
    SNT=SIN(TT)
    CSB=COS(BTA)
    SNB=SIN(BTA) 
    X1=FABSC(-SWD1X,SVD1X,BTA)
    X2=FABSC(SVD1X,SWD1X,BTA)
    X1D1T=FABSC(-SWD1XD1T,SVD1XD1T,BTA)
    X2D1T=FABSC(SVD1XD1T,SWD1XD1T,BTA)
    TACP=X1*X2
    TACPD1T=X1D1T*X2+X1*X2D1T

      Y1=FABSC(-SVD1X,-SWD1X,TT)
      Y2=FABSC(-SWD1X,SVD1X,TT)
      Y1D1T=FABSC(-SVD1XD1T,-SWD1XD1T,TT)
      Y2D1T=FABSC(-SWD1XD1T,SVD1XD1T,TT)

    TDE(1,1)=1.0D0
    TDE(1,2)=SVD1X
    TDE(1,3)=SWD1X
    TDE(2,1)=Y1
    TDE(2,2)=CST
    TDE(2,3)=SNT
    TDE(3,1)=Y2
    TDE(3,2)=-SNT+TACP*CSB
    TDE(3,3)=CST+TACP*SNB

    TDED1T(1,1)=0.0D0
    TDED1T(1,2)=SVD1XD1T
    TDED1T(1,3)=SWD1XD1T
    TDED1T(2,1)=Y1D1T+SPHD1T*Y2
    TDED1T(2,2)=-SPHD1T*SNT
    TDED1T(2,3)=SPHD1T*CST
    TDED1T(3,1)=Y2D1T-SPHD1T*Y1
    TDED1T(3,2)=-SPHD1T*CST+TACPD1T*CSB
    TDED1T(3,3)=-SPHD1T*SNT+TACPD1T*SNB

END SUBROUTINE GETTDE
!->+++++++++++++++++++++++++++++++++++++++++++++



!->+++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> the modified Gram-Schmidt process.
!-> n <= m.
subroutine GramSchmidtPMd(m, n, A, AUM) 
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
integer, intent(in):: m, n
real(RDT), intent(in):: A(m, n)
real(RDT), intent(out):: AUM(m, n)

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
REAL(RDT), EXTERNAL:: DNRM2, DDOT
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->        
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
character(80) filenam
integer i, j
integer incx, incy
real(rdt) xnm, xdd

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

      incx=1
      incy=1
       
      AUM=A
      do 10 i=1, n
            do 20 j=1, i-1
                  xdd=DDOT(m, AUM(:, i), incx, AUM(:, j), incy)
                  AUM(:, i)=AUM(:, i)-xdd*AUM(:, j)       
20        continue
            xnm=DNRM2(m, AUM(:, i), 1)
            AUM(:, i)=AUM(:, i)/xnm            
10  continue

end subroutine GramSchmidtPMd
!->+++++++++++++++++++++++++++++++++++++++++++++



!->+++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> the [classical] Gram-Schmidt process.
!-> n <= m.
subroutine GramSchmidtP(m, n, A, AUM) 
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
integer, intent(in):: m, n
real(RDT), intent(in):: A(m, n)
real(RDT), intent(out):: AUM(m, n)

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
REAL(RDT), EXTERNAL:: DNRM2, DDOT
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->        
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
character(80) filenam
integer i, j
integer incx, incy
real(rdt) xnm, xdd

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

      incx=1
      incy=1
       
      AUM=A
      do 10 i=1, n
            do 20 j=1, i-1
                  xdd=DDOT(m, A(:, i), incx, AUM(:, j), incy)
                  AUM(:, i)=AUM(:, i)-xdd*AUM(:, j)       
20        continue
            xnm=DNRM2(m, AUM(:, i), 1)
            AUM(:, i)=AUM(:, i)/xnm            
10  continue

end subroutine GramSchmidtP
!->+++++++++++++++++++++++++++++++++++++++++++++



!->+++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> to solves overdetermined [or underdetermined] real linear systems.
subroutine LeastSquareM(XS, mm, nn, pp, AA, BB, lsmd)
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
integer, intent(in):: mm, nn, pp, lsmd
real(RDT), intent(in):: AA(mm, nn), BB(mm, pp)
real(RDT), intent(out):: XS(nn, pp)

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->       DGELSY, DGELS
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
CHARACTER TRANS

integer INFO, LDA, LDB, LWORK, M, N, NRHS, RANK
integer mn
integer err, nout

real(RDT) RCOND

integer, allocatable:: JPVT(:)
real(RDT), allocatable:: A(:, :), B(:, :)
real(RDT), allocatable:: WORK(:)

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

      M=mm 
      N=nn
      NRHS=pp
      LDA=m 
      LDB=m
      mn=min0(M,N)
     
      goto(96, 97) lsmd      
96  continue   

      LWORK=MAX( MN+3*N+1, 2*MN+NRHS )    
      
      allocate(A(mm, nn), B(mm, pp), &
            work(lwork), & 
            JPVT(nn), & 
            stat=err)
    if(err.ne.0) then
        write(*, 998)
          stop
      end if

      A=AA
      B=BB 
      
      RCOND=1.0D-10      
      JPVT=0  
      call DGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, &
            WORK, LWORK, INFO )
      goto 100
      
97  continue      
      TRANS = 'N' 
      LWORK=MN+ max0( MN, NRHS ) 
      
      allocate(A(mm, nn), B(mm, pp), &
            work(lwork), & 
            JPVT(nn), & 
            stat=err)
    if(err.ne.0) then
        write(*, 998)
          stop
      end if

      A=AA
      B=BB 
      
      call DGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, &
            INFO )    
      goto 100
            
100 continue
      XS(1:nn, 1:pp)=B(1:nn, 1:pp)
      
      deallocate(A, B, &
            work, & 
            JPVT, & 
            stat=err)
    if(err.ne.0) then
        write(*, 999)
          stop
      end if      
      return 
      
998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')   
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')     
end subroutine LeastSquareM
!->+++++++++++++++++++++++++++++++++++++++++++++



!->+++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> solve the linear systems (amt*x=p) with constraints: 
!->       x(1:m)=x(n-m+1:n).   [///p(1:m)=p(n-m+1:n)///]
!->       n > 2m.
subroutine lincon2(x, n, m, amt, p)
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
integer, intent(in):: n, m
real(RDT), intent(in):: amt(n, n)
real(RDT), intent(in):: p(n)
real(RDT), intent(out):: x(n)

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->       mtxinv2, MAMU, linsolver2
!->       DGEMM, DGEMV          
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
integer, parameter:: ione=1
integer np, nc
integer err
integer lsmd
real(RDT), allocatable:: ax(:, :)

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> main body:
      
      np=n-m

    ALLOCATE(ax(n, np), &
            STAT=ERR)

    IF(ERR.NE.0) THEN
        WRITE(*, 998)
        STOP
    END IF
    
    ax=amt(1:n, 1:np)
    ax(1:n, 1:m)=ax(1:n, 1:m)+amt(1:n, np+1:n)

      lsmd=1
      nc=1 
      call LeastSquareM(x, n, np, nc, ax, p, lsmd)
      x(np+1:n)=x(1:m) 
      
    DEALLOCATE(ax, &
            STAT=ERR)

    IF(ERR.NE.0) THEN
        WRITE(*, 999)
        STOP
    END IF     
     
998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
end subroutine lincon2
!->+++++++++++++++++++++++++++++++++++++++++++++



!->+++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> solve the linear systems (amt*x=p) with constraints: 
!->       x(1:m)=x(n-m+1:n).   [///p(1:m)=p(n-m+1:n)///]
!->       n > 2m.
!->   this routine is obsolete.
subroutine lincon(x, n, m, amt, p)
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
integer, intent(in):: n, m
real(RDT), intent(in):: amt(n, n)
real(RDT), intent(in):: p(n)
real(RDT), intent(out):: x(n)

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->       mtxinv2, MAMU, linsolver2
!->       DGEMM, DGEMV          
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
integer, parameter:: ione=1
integer np
integer err
real(RDT), parameter:: oned=1.0D0
real(RDT), allocatable:: a11(:, :), a12(:, :), a21(:, :), a22(:, :)
real(RDT), allocatable:: a12x(:, :)
real(RDT), allocatable:: p1(:), p2(:)

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

      
      np=n-m

    ALLOCATE(a11(np, np), a12(np, m), a21(m, np), a22(m, m), &
          a12x(np, m), &
          p1(np), p2(m), &
            STAT=ERR)

    IF(ERR.NE.0) THEN
        WRITE(*, 998)
        STOP
    END IF
    
    a11=amt(1:np, 1:np)
    a12=amt(1:np, np+1:n)
    a21=amt(np+1:n, 1:np)
    a22=amt(np+1:n, np+1:n)
    p1=p(1:np)
    p2=p(np+1:n)
    
      call mtxinv2(a22, 1, m, a22) 
      write(*, *) 989 
      call MAMU(a12x, np, m, m, a12, a22)
      call DGEMM('N', 'N', np, np, m, -oned, a12x, np, a21, m, oned, a11, np)
      call DGEMV('N', np, m, -oned, a12x, np, p2, ione, oned, p1, ione) 

      call linsolver2(x, 1, np, a11, p1)  
      x(np+1:n)=x(1:m) 
      write(*, *) 9892
      
    DEALLOCATE(a11, a12, a21, a22, &
          a12x, &
          p1, p2, & 
            STAT=ERR)

    IF(ERR.NE.0) THEN
        WRITE(*, 999)
        STOP
    END IF     
     
998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
end subroutine lincon
!->+++++++++++++++++++++++++++++++++++++++++++++



!->+++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> to solve the systems of linear equations.
subroutine linsolver2(xsol, ichs, mdm, mtrx, bvct)
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
!  include 'link_f90_static.h'
!USE LIN_SOL_GEN_INT
!USE LSARG_INT
!USE LSLRG_INT
!USE LSADS_INT
!USE LSLDS_INT
!USE WRRRN_INT
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
integer, intent(in):: ichs, mdm
real(rdt), intent(in):: mtrx(mdm, mdm)
real(rdt), intent(in):: bvct(mdm)
real(rdt), intent(out):: xsol(mdm)

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
DOUBLE PRECISION, EXTERNAL:: DLANGE
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->       DGETRF, DGECON, DGETRS, DGERFS
!->       DPOTRF, DPOCON, DPOTRS, DPORFS
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
CHARACTER UPLO, NORMS, TRANS
INTEGER INFO, LDA, N, M, LWORK, NRHS, LDB
INTEGER LDAF, LDX

integer err

integer, allocatable:: IPIV(:), IWORK(:)


DOUBLE PRECISION ANORM, RCOND
DOUBLE PRECISION FERR, BERR

DOUBLE PRECISION, allocatable:: WORK(:)
DOUBLE PRECISION, allocatable::  A(:, :), B(:, :), AF(:, :), X(:, :)



      ALLOCATE(A( mdm, mdm), B( mdm, 1), AF( mdm, mdm), X(mdm, 1), &
            WORK( 4*mdm ), & 
            IPIV(mdm), IWORK(mdm), & 
            STAT=ERR)

      IF(ERR.NE.0) THEN
            WRITE(*, 998)
            STOP
      END IF
    

   ! bmtx(:, 1)=bvct(:)
  !  CALL D_LIN_SOL_GEN(mtrx, bmtx, xmtx)
   ! xsol(:)=xmtx(:, 1)
    ! solves A(m, m) * X(m, n)=B(m, n)
    
     A=mtrx      
         
     M=mdm
     N=mdm
     LDA=mdm     
     
     NORMS='1'
     
        ANORM=DLANGE( NORMS, M, N, A, LDA, WORK ) 

     goto(5, 15) ichs
!*   ichs=1, real general system; 
!*   ichs=2, real symmetric positive definite system of linear equations.   
  
5   continue     
     call DGETRF( M, N, A, LDA, IPIV, INFO )

     call DGECON( NORMS, N, A, LDA, ANORM, RCOND, WORK, IWORK, &
                              INFO )    
          
     TRANS='N'
     NRHS=1
     B(:, 1)=bvct(:)
     LDB=mdm
     call DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )     
 
     LDAF=mdm
     LDX=mdm
     AF=A
     A=mtrx
     X=B
     B(:, 1)=bvct(:)
     call DGERFS( TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, &
                             X, LDX, FERR, BERR, WORK, IWORK, INFO )  
     xsol=x(:, 1)                                                     
    !-> Solves a real general system of linear equations with iterative refinement.
    return
   
   
    
    ! CALL D_LSLRG(mtrx, bvct, xsol)
    !-> Solves a real general system of linear equations without iterative refinement.



15  continue   
     UPLO='L'
     call DPOTRF( UPLO, N, A, LDA, INFO )   

     call DPOCON( UPLO, N, A, LDA, ANORM, RCOND, WORK, IWORK, &
                              INFO )            

     NRHS=1
     B(:, 1)=bvct(:)
     LDB=mdm
     call DPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO )           

     LDAF=mdm
     LDX=mdm
     AF=A
     A=mtrx
     X=B
     B(:, 1)=bvct(:)
     call DPORFS( UPLO, N, NRHS, A, LDA, AF, LDAF, B, LDB, X, &
                             LDX, FERR, BERR, WORK, IWORK, INFO ) 
     xsol=x(:, 1)                       
    !-> Solves a real symmetric positive definite system of linear equations with iterative refinement.
    return
   
    
    !CALL D_LSLDS(mtrx, bvct, xsol)
    !-> Solves a real symmetric positive definite system of linear equations without iterative refinement. 
      
   !  CALL WRRRN ('X', xsol, 1, mdm, 1)

    DEALLOCATE(A, B, AF, X, &
            WORK, & 
            IPIV, IWORK, & 
            STAT=ERR)

    IF(ERR.NE.0) THEN
        WRITE(*, 999)
        STOP
    END IF     
     
998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
end subroutine linsolver2
!->+++++++++++++++++++++++++++++++++++++++++++++



SUBROUTINE low2full(n, A)
implicit none

INTEGER, INTENT(IN):: N
DOUBLE PRECISION, INTENT(INOUT):: A(N, N)

integer i, j

     do 5 i=1, n-1
          do 10 j=i+1, n
               A(i, j)=A(j, i)
10      continue           
5   continue    

end SUBROUTINE low2full



!->+++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> THE MULTIPLICATION OF TWO MATRICES.
!-> CAN ALSO USE THE INTRINSIC FUNCTION <MATMUL> TO GET THE RESULT.
SUBROUTINE MAMU(MR,M,P,N,M1,M2)
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
INTEGER,INTENT(IN)::M,P,N
REAL(RDT),INTENT(IN)::M1(M,P),M2(P,N)
REAL(RDT),INTENT(OUT)::MR(M,N)

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->       
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
INTEGER I,J,K

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

!$OMP PARALLEL DO PRIVATE(I,J,K)
      DO 10 I=1,M
        DO 20 J=1,N
            MR(I,J)=0.0D0
            DO 30 K=1,P
                MR(I,J)=MR(I,J)+M1(I,K)*M2(K,J)
30            CONTINUE
20      CONTINUE
10  CONTINUE
!$OMP END PARALLEL DO

END SUBROUTINE MAMU
!->+++++++++++++++++++++++++++++++++++++++++++++



!->+++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> THE MULTIPLICATION OF TWO MATRICES.
!-> CAN ALSO USE THE INTRINSIC FUNCTION <MATMUL> TO GET THE RESULT.
SUBROUTINE MAMU2(MR,M,P,N,M1,M2)
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
INTEGER,INTENT(IN)::M,P,N
REAL(RDT),INTENT(IN)::M1(M,P),M2(P,N)
REAL(RDT),INTENT(OUT)::MR(M,N)

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->       
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
INTEGER ERR
INTEGER I,J,K
REAL(RDT), ALLOCATABLE:: mt1(:, :), mt2(:, :)

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

    ALLOCATE(MT1(M,P), MT2(P,N), &
            STAT=ERR)

    IF(ERR.NE.0) THEN
        WRITE(*,998)
        STOP
    END IF
    
    MT1=M1
    MT2=M2
      CALL MAMU(MR,M,P,N,MT1,MT2)


    DEALLOCATE(MT1, MT2, &
            STAT=ERR)

    IF(ERR.NE.0) THEN
        WRITE(*,999)
        STOP
    END IF     
     
998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
END SUBROUTINE MAMU2
!->+++++++++++++++++++++++++++++++++++++++++++++



!->+++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> THE consecutive [/successive] MULTIPLICATION OF THREE MATRICES.
subroutine MATMLYTP(mtt, nm, np, nq, nn, mta, mtb, mtc) 
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
integer, intent(in):: nm, np, nq, nn
real(RDT), intent(in):: mta(nm, np), mtb(np, nq), mtc(nq, nn)
real(RDT), intent(out):: mtt(nm, nn)

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->       MAMU
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
REAL(RDT), ALLOCATABLE:: mtbc(:,:)
integer err

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

    ALLOCATE(mtbc(np, nn), &
            STAT=ERR)

    IF(ERR.NE.0) THEN
        WRITE(*,998)
        STOP
    END IF
    
    
      call MAMU(mtbc,np,nq,nn,mtb,mtc)
      call MAMU(mtt,nm,np,nn,mta,mtbc)
     
     
    DEALLOCATE(mtbc, &
            STAT=ERR)

    IF(ERR.NE.0) THEN
        WRITE(*,999)
        STOP
    END IF     
     
998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
end subroutine MATMLYTP
!->+++++++++++++++++++++++++++++++++++++++++++++



!->+++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> the inverse of a matrix.
subroutine mtxinv2(invmtx, ichs, mdm, mtrx)
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
integer, intent(in):: ichs, mdm
real(rdt), intent(in):: mtrx(mdm, mdm)
real(rdt), intent(out):: invmtx(mdm, mdm)

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->       MAMU
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
CHARACTER UPLO
INTEGER INFO, LDA, N, M, LWORK
INTEGER IPIV( mdm)
integer err
DOUBLE PRECISION, ALLOCATABLE:: A(:, :), work(:)

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

    ALLOCATE(A(mdm, mdm), work(mdm), &
            STAT=ERR)

    IF(ERR.NE.0) THEN
        WRITE(*,998)
        STOP
    END IF
    
    
     A=mtrx 
     
     goto(5, 15) ichs
!*   ichs=1, real general system; 
!*   ichs=2, real symmetric positive definite system of linear equations.   
  
5   continue     
     M=mdm
     N=mdm
     LDA=mdm
     LWORK=mdm
     call DGETRF( M, N, A, LDA, IPIV, INFO )
     if(info .eq. 0) then
          call DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
          if(info .eq. 0) then     
               invmtx=A   
          else
               write(*, *) " ERROR!  566."
               pause        
          end if
     else
          write(*, *) " ERROR!  567."
          pause
     end if 
    !-> for the general matrix.
     return
   
  
15  continue   
     UPLO='L'
     N=mdm  
     LDA=mdm    
     call DPOTRF( UPLO, N, A, LDA, INFO )    
     CALL DPOTRI( UPLO, N, A, LDA, INFO )  
     call low2full(mdm, A)   
     invmtx=A  
     !-> for the real symmetric positive definite matrix.
     return      

    DEALLOCATE(A, work, &
            STAT=ERR)

    IF(ERR.NE.0) THEN
        WRITE(*,999)
        STOP
    END IF     
     
998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
end subroutine mtxinv2
!->+++++++++++++++++++++++++++++++++++++++++++++



!->+++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> THE MULTIPLICATION OF A MATRIX AND A VECTOR.
SUBROUTINE MVMU(MR,M,N,MT,VT)
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
INTEGER,INTENT(IN)::M,N
REAL(RDT),INTENT(IN)::MT(M,N),VT(N)
REAL(RDT),INTENT(OUT)::MR(M)

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->       
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
INTEGER I,K

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

      DO 10 I=1,M
        MR(I)=0.0D0
        DO 20 K=1,N
            MR(I)=MR(I)+MT(I,K)*VT(K)
20      CONTINUE
10  CONTINUE

END SUBROUTINE MVMU
!->+++++++++++++++++++++++++++++++++++++++++++++



!->+++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> THE MULTIPLICATION OF A MATRIX AND A VECTOR.
SUBROUTINE MVMU1(MR,M,N,MT,VT)
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
INTEGER,INTENT(IN)::M,N
REAL(RDT),INTENT(IN)::MT(M,N),VT(N)
REAL(RDT),INTENT(OUT)::MR(M)

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->       
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
INTEGER I,K

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

      MR=0.0D0
      DO 10 K=1,N
            DO 20 I=1,M
                MR(I)=MR(I)+MT(I,K)*VT(K)
20      CONTINUE
10  CONTINUE

END SUBROUTINE MVMU1
!->+++++++++++++++++++++++++++++++++++++++++++++



!->+++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> THE MULTIPLICATION OF A MATRIX AND A VECTOR.
SUBROUTINE MVMU2(MR,M,N,MT,VT)
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
INTEGER,INTENT(IN)::M,N
REAL(RDT),INTENT(IN)::MT(M,N),VT(N)
REAL(RDT),INTENT(OUT)::MR(M)

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->       MVMU1  
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
REAL(RDT), ALLOCATABLE:: VTT(:)
INTEGER I,K
integer err

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

    ALLOCATE(VTT(N), &
            STAT=ERR)

    IF(ERR.NE.0) THEN
        WRITE(*, 998)
        STOP
    END IF
    
      VTT=VT            
      CALL MVMU1(MR,M,N,MT,VTT)

    DEALLOCATE(VTT, &
            STAT=ERR)

    IF(ERR.NE.0) THEN
        WRITE(*, 999)
        STOP
    END IF     
     
998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
END SUBROUTINE MVMU2
!->+++++++++++++++++++++++++++++++++++++++++++++



!->+++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> output a (m x n) matrix.
subroutine outmatx(filenam, m, n, mtx)
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
character*(*), intent(in):: filenam
integer, intent(in):: m, n
real(RDT), intent(in):: mtx(m, n)

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->       
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
character(len=80) copyfilenam
integer i, j

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

      write(copyfilenam,*) '../', trim(filenam), '.dat'
      !>open(unit=97, file=copyfilenam)
      open(97, file=copyfilenam)


      WRITE(97, *)  m, n
      WRITE(97, *) 
      do 588 i=1, m
            WRITE(97, 993) i, (mtx(i, j), j=1, n)
588 continue

      close(97)
993 format(I12, 1000f36.12)      
end subroutine outmatx
!->+++++++++++++++++++++++++++++++++++++++++++++



!->+++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> output a (m x n) matrix.
subroutine outmatxint(filenam, m, n, mtx)
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
character*(*), intent(in):: filenam
integer, intent(in):: m, n
integer, intent(in):: mtx(m, n)

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->       
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
character(len=80) copyfilenam
integer i, j

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

      write(copyfilenam,*) '../', trim(filenam), '.dat'
      !>open(unit=97, file=copyfilenam)
      open(97, file=copyfilenam)


      WRITE(97, *)  m, n
      WRITE(97, *) 
      do 588 i=1, m
            WRITE(97, 992) i, (mtx(i, j), j=1, n)
588 continue

      close(97)
992 format(1000I12)      
end subroutine outmatxint
!->+++++++++++++++++++++++++++++++++++++++++++++



subroutine rotationrateT(omg, v)
use GlobalDataFun
implicit none

integer, parameter:: dm=3

real(rdt), intent(in):: v(dm)
real(rdt), intent(out):: omg(dm, dm)

     omg=0.0D0
     omg(1, 2)=-v(3)
     omg(1, 3)=v(2)
     omg(2, 1)=v(3)
     omg(2, 3)=-v(1)
     omg(3, 1)=-v(2)
     omg(3, 2)=v(1)
   
end subroutine rotationrateT



!* DEL:DELTA
!* TM:TRANSFORM MATRIX
SUBROUTINE TXYZ(TM,L,DEL)
use GlobalDataFun
implicit none

INTEGER,INTENT(IN)::L
REAL(RDT),INTENT(IN)::DEL
REAL(RDT),INTENT(OUT)::TM(3,3)

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++

REAL(RDT) A,B
        
    TM=0.0D0
    A=COS(DEL)
    B=SIN(DEL)

    GET_AXIS: SELECT CASE(L)
    CASE(1)
        TM(1,1)=1.0D0            
        TM(2,2)=A
        TM(2,3)=B
        TM(3,2)=-B
        TM(3,3)=A
    CASE(2)          
        TM(1,1)=A            
        TM(1,3)=-B
        TM(2,2)=1.0D0
        TM(3,1)=B
        TM(3,3)=A
    CASE(3)
        TM(1,1)=A            
        TM(1,2)=B
        TM(2,1)=-B
        TM(2,2)=A
        TM(3,3)=1.0D0
    END SELECT GET_AXIS
    
END SUBROUTINE TXYZ



!* DEL:DELTA
!* TM:TRANSFORM MATRIX
SUBROUTINE TXYZD1T(TM,L,DEL,DELD1T)
use GlobalDataFun
implicit none

INTEGER,INTENT(IN)::L
REAL(RDT),INTENT(IN)::DEL,DELD1T
REAL(RDT),INTENT(OUT)::TM(3,3)

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++

REAL(RDT) A,B
        
    TM=0.0D0
    A=-SIN(DEL)
    B=COS(DEL)

    GET_AXIS: SELECT CASE(L)
    CASE(1)          
        TM(2,2)=A
        TM(2,3)=B
        TM(3,2)=-B
        TM(3,3)=A
      CASE(2)          
        TM(1,1)=A            
        TM(1,3)=-B
        TM(3,1)=B
        TM(3,3)=A
      CASE(3)
        TM(1,1)=A            
        TM(1,2)=B
        TM(2,1)=-B
        TM(2,2)=A
      END SELECT GET_AXIS
    
    TM=DELD1T*TM

END SUBROUTINE TXYZD1T



!* THE CROSS PRODUCT OF THE 2 3-DIMENSIONAL VECTORS.
SUBROUTINE VCP(VR,V1,V2)
use GlobalDataFun
implicit none

INTEGER, PARAMETER:: DM=3
REAL(RDT),INTENT(IN)::V1(DM),V2(DM)
REAL(RDT),INTENT(OUT)::VR(DM)

    VR(1)=V1(2)*V2(3)-V1(3)*V2(2)
    VR(2)=V1(3)*V2(1)-V1(1)*V2(3)
    VR(3)=V1(1)*V2(2)-V1(2)*V2(1)

END SUBROUTINE VCP
 
      
      



!* THE OUTER PRODUCT OF THE TWO VECTORS.
SUBROUTINE VOP(MR,DM1,DM2,V1,V2)
use GlobalDataFun
implicit none

INTEGER,INTENT(IN)::DM1,DM2
REAL(RDT),INTENT(IN)::V1(DM1),V2(DM2)
REAL(RDT),INTENT(OUT)::MR(DM1,DM2)

INTEGER I,J

    DO 10 J=1,DM2
        DO 5 I=1,DM1
            MR(I,J)=V1(I)*V2(J)
5        CONTINUE
10  CONTINUE    
    
END 


SUBROUTINE VOP2(MR,DM1,DM2,V1,V2)
use GlobalDataFun
implicit none

INTEGER, INTENT(IN):: DM1,DM2
REAL(RDT), INTENT(IN):: V1(DM1),V2(DM2)
REAL(RDT),INTENT(OUT)::MR(DM1,DM2)

INTEGER I,J

    DO 5 J=1,DM2
        DO 10 I=1,DM1
            MR(I,J)=V1(I)*V2(J)
10      CONTINUE
5   CONTINUE    
    
END SUBROUTINE VOP2


!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
function paai() 
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
real(rdt) paai

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->        
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
real(rdt), parameter:: oned=1.0D0

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

      paai=ACOS(-oned)

end function paai
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> the average or mean value of array.
subroutine mean2(dmn, atm, m, n, p, mk) 
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
integer, intent(in):: m, n, p, mk
real(rdt), intent(in):: atm(m, n)
real(rdt), intent(out):: dmn(p)

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->       mean      
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
integer i
integer err
real(rdt), allocatable:: atsp(:, :)

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

      if(m .lt. 1 .or. n .lt. 1) then 
            goto 100 
      end if
       
      if(mk .eq. 1) then 
            if(p .ne. m) then 
                  goto 100 
            else
                  call  mean(dmn, atm, m, n) 
                  return
            end if  
      else if(mk .eq. 2 ) then 
            if(p .ne. n) then 
                  goto 100 
            else
                  ALLOCATE( atsp(n, m), &
                        STAT=ERR)
                              
                IF(ERR.NE.0) THEN
                    WRITE(*, 998)
                    STOP
                END IF 
                
                  atsp=transpose(atm) 
                  call  mean(dmn, atsp, n, m) 

                  DEALLOCATE( atsp, &
                        STAT=ERR)
                              
                IF(ERR.NE.0) THEN
                    WRITE(*, 999)
                    STOP
                END IF
    
                  return
            end if  
      else
            goto 100 
      end if   

100 continue
      write(*, *) 'error.', 900010
      return       

998 FORMAT(1X,'ALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
999 FORMAT(1X,'DEALLOCATING ARRAY ERROR,PROGRAM STOPPED!')
end subroutine mean2
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> the average or mean value of array.
subroutine mean(dmn, atm, m, n) 
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
integer, intent(in):: m, n
real(rdt), intent(in):: atm(m, n)
real(rdt), intent(out):: dmn(m)

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->        
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
integer i

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

      dmn=0.0D0
      do 5 i=1, n
            dmn(:)=dmn(:)+atm(:, i)     
5    continue      
      dmn=dmn/dble(n)
      
end subroutine mean
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> Prandtl's tip-loss function.
subroutine prandtiplossf1(fbg, nb, rb, ladi) 
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
implicit none

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
subroutine prandtiplossf(fbg, nb, rb, ladi) 
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
implicit none

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
real(rdt) f,b

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> main body:
      b=1-ladi/dble(nb)*2*log(2.0)
      f=0.5d0**((1.0D0-rb)/(1.0D0-b))
      fbg=2.0D0*ACOS(f)/paai() 
      
end subroutine prandtiplossf
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
FUNCTION FWM(KL,KR,WWT)
use GlobalDataFun
implicit none

REAL(RDT), INTENT(IN):: KL,KR,WWT
REAL(RDT) FWM

    FWM=KL*(1.0-WWT)+KR*WWT
    !->++++++++++++++
      !FWM=KL
END FUNCTION FWM

FUNCTION FABSC(KC,KS,BETA)
use GlobalDataFun
implicit none

REAL(RDT), INTENT(IN):: KC,KS,BETA
REAL(RDT) FABSC

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->       
!->+++++++++++++++++++++++++++++++++++++++++++++

    FABSC=KC*COS(BETA)+KS*SIN(BETA)

END FUNCTION FABSC



!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
subroutine rxlct(imk, ss, ll, n, rx, rl) 
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
integer, intent(in):: n
real(rdt), intent(in):: rx, rl(n)
integer, intent(out):: imk
real(rdt), intent(out):: ss, ll

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->        
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
integer i, inx

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

      imk=-1
      do 5 i=1, n-1
            inx=i+1 
            if( rl(i) .le. rx .and. rl(inx) .ge. rx ) then   
                  imk=i
                  ll=rl(inx)-rl(i)
                  ss=(rx-rl(i))/ll 
                  exit  
            end if    
5   continue

    if( imk .eq. -1.OR.imk.EQ.n ) then
            WRITE(*,*) RX
            WRITE(*,*) RL
          write(*, *) "Error!", 999963
          PAUSE
          stop
    end if
            
end subroutine rxlct
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
subroutine rxlct2(imk, ss, n, rx, rl) 
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
integer, intent(in):: n
real(rdt), intent(in):: rx, rl(n)
integer, intent(out):: imk
real(rdt), intent(out):: ss

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->        
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
integer i, inx
real(rdt) :: ll
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

      imk=-1
      do 5 i=1, n-1
            inx=i+1 
            if( rl(i) .le. rx .and. rl(inx) .ge. rx ) then   
                  imk=i
                  ll=rl(inx)-rl(i)
                  ss=(rx-rl(i))/ll 
                  exit  
            end if    
5   continue

    if( imk .eq. -1.and.rx.lt.rl(1)) then
        !if(ABS((rx-rl(1))/rx).gt.0.05) then
        !    write(*,*) 'double warning!'
        !    write(*,*) RX
        !    WRITE(*,*) RL(1)
        !    !PAUSE
        !end if
        !write(*,*) 'warning!'
        !WRITE(*,*) RX
        !WRITE(*,*) RL(1)
        imk=1
        ss=0.0
    else if(imk.EQ.-1.and.rx.gt.rl(n)) then
        !if(abs((rx-rl(n))/rx).gt.0.05) then
        !    write(*,*) 'double warning!'
        !    write(*,*) RX
        !    WRITE(*,*) RL(n)
        !    !PAUSE
        !end if
        !write(*,*) 'warning!'
        !WRITE(*,*) RX
        !WRITE(*,*) RL(n)
        imk=n-1
        ss=1.0
    end if
            
end subroutine rxlct2
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
FUNCTION INTOP(N,CL,AIM,LENTH)
use GlobalDataFun
IMPLICIT NONE
INTEGER :: N
REAL(RDT),INTENT(IN) :: AIM,LENTH(N),CL(N)
REAL(RDT)  :: INTOP
INTEGER :: LBI,LBIP1
REAL(RDT) :: SS
REAL(RDT),EXTERNAL :: FWM
    call rxlct2(LBI, ss, N, AIM, LENTH) 

    LBIP1=LBI+1
    INTOP=FWM(CL(LBI),CL(LBIP1),SS)
    
END FUNCTION

SUBROUTINE INTP2D(CXF,SSI,SSJ,VIM1JM1,VIJM1,VIM1J,VIJ)
use GlobalDataFun
IMPLICIT NONE
real(rdt),INTENT(IN) :: SSI,SSJ,VIM1JM1,VIJM1,VIM1J,VIJ
real(rdt),INTENT(OUT) :: CXF

real(rdt) :: COEIM1JM1,COEIM1J,COEIJM1,COEIJ


    COEIM1JM1=(1.0-SSJ)*(1.0-SSI)
    COEIM1J=SSJ*(1.0-SSI)
    COEIJM1=(1.0-SSJ)*SSI
    COEIJ=SSJ*SSI
    
    IF(ABS(COEIM1JM1+COEIM1J+COEIJM1+COEIJ).LT.1D-6) THEN
        WRITE(*,*) 'INTP2D EEROR'
    END IF
    
    CXF=COEIM1JM1*VIM1JM1+COEIM1J*VIM1J+COEIJM1*VIJM1+COEIJ*VIJ

END SUBROUTINE

SUBROUTINE INTP3D(CXF,SSI,SSJ,SSK,VIM1JM1KM1,VIJM1KM1,VIM1JKM1,VIJKM1,VIM1JM1K,VIJM1K,VIM1JK,VIJK)
use GlobalDataFun
IMPLICIT NONE
real(rdt),INTENT(IN) :: SSI,SSJ,SSK,VIM1JM1KM1,VIJM1KM1,VIM1JKM1,VIJKM1,VIM1JM1K,VIJM1K,VIM1JK,VIJK
real(rdt),INTENT(OUT) :: CXF
    
    CXF=(1.0-SSK)*(1.0-SSJ)*(1.0-SSI)*VIM1JM1KM1+&
		(1.0-SSK)*(1.0-SSJ)*SSI*VIJM1KM1+&
		(1.0-SSK)*SSJ*(1.0-SSI)*VIM1JKM1+&
		(1.0-SSK)*SSJ*SSI*VIJKM1+&
		SSK*(1.0-SSJ)*(1.0-SSI)*VIM1JM1K+&
		SSK*(1.0-SSJ)*SSI*VIJM1K+&
		SSK*SSJ*(1.0-SSI)*VIM1JK+&
		SSK*SSJ*SSI*VIJK

END SUBROUTINE

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> C = v + A*B.
subroutine abplsv(rc, dm, dn, npt, rea, ted, cdlct) 
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
!-> dm == 3.
integer, intent(in):: dm, dn, npt
real(rdt), intent(in):: rea(dm), ted(dm, dn), cdlct(dn, npt)
real(rdt), intent(out):: rc(dm, npt)

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->       MVMU2        
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
integer i
real(rdt) vtp1(dn), vtp2(dm)

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

      do 5 i=1, npt
            vtp1=cdlct(:, i)  
            CALL MVMU2(vtp2, dm, dn, ted, vtp1) 
            rc(:, i)=rea+vtp2       
5    continue      

end subroutine abplsv 
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
subroutine avadd(a, v, m, n, p, mk) 
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
integer, intent(in):: m, n, p, mk
real(rdt), intent(in):: v(p)
real(rdt), intent(inout):: a(m, n)

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->        
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
integer i

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

      if(m .lt. 1 .or. n .lt. 1) then 
            goto 100 
      end if
       
      if(mk .eq. 1) then 
            if(p .ne. m) then 
                  goto 100 
            else
                  do 5 i=1, n
                        a(:, i)=a(:, i)+v(:) 
5                continue      
                  return
            end if  
      else if(mk .eq. 2 ) then 
            if(p .ne. n) then 
                  goto 100 
            else
                  do 15 i=1, m
                        a(i, :)=a(i, :)+v(:) 
15              continue    
                  return
            end if  
      else
            goto 100 
      end if   

100 continue
      write(*, *) 'error.', 900011
      return       
      
end subroutine avadd
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> to solve nonlinear equations , f(x)=0.
!-> Newton-Raphson method, single step
subroutine newtonrap_sp(xv, n, fxv, jacb, kk)
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
integer, intent(in):: n
real(rdt), intent(inout):: xv(n)
real(rdt), intent(in):: fxv(n), jacb(n, n), kk

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->        
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
real(rdt) dlxv(n), jacbiv(n, n)

!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

   
      call linsolver2(dlxv, 1, n, jacb, -fxv) 
      
      xv=xv+kk*dlxv
      
      return
end subroutine newtonrap_sp
!->+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!************************************************************
FUNCTION RtoChar(n) RESULT(char)
use GlobalDataFun
    REAL(rdt),INTENT(IN):: n
    CHARACTER(100):: char

     WRITE(char, '(F6.3)') n

END FUNCTION RtoChar
!***********************************************************
 !->+++++++++++++++++++++++++++++++++++++++++++++
subroutine dichotomy(x,tree,m,n)
use GlobalDataFun
implicit none
! input
REAL(KIND=rdt),intent(in):: x ! 需要确定位置的数值
integer,intent(in):: m ! 数组序列的维数
REAL(KIND=rdt),intent(in):: tree(m) ! 数组序列
! output
integer,intent(out):: n ! x在数组序列中的位置

integer:: i,s_b,s_e,s_mid

n=-1

if(tree(1)<tree(m)) then ! 升序
s_b=1
s_e=m
do i=1,1000
  s_mid=(s_e+s_b)/2
  if(x>=tree(s_b).and.x<=tree(s_mid)) then
    s_e=s_mid
  else
    s_b=s_mid
  end if
  if(s_e-s_b==1)exit
end do
n=s_b
else  ! 降序
s_b=1
s_e=m
do i=1,1000
  s_mid=(s_e+s_b)/2
  if(x<=tree(s_b).and.x>=tree(s_mid)) then
    s_e=s_mid
  else
    s_b=s_mid
  end if
  if(s_e-s_b==1)exit
end do
n=s_b
end if

end subroutine dichotomy
FUNCTION FIND_MAX(N,MV,M)
USE GlobalDataFun
IMPLICIT NONE
REAL(RDT) FIND_MAX
INTEGER,INTENT(OUT) :: M
INTEGER,INTENT(IN) :: N
REAL(RDT),INTENT(IN) :: MV(N)

REAL(RDT) TEMP
INTEGER I
    TEMP=MV(1)
    DO I=2,N
        IF(TEMP.LT.MV(I)) THEN
            M=I
            TEMP=MV(I)
        END IF
    END DO
    
    FIND_MAX=TEMP
END FUNCTION


FUNCTION FIND_MIN(N,MV,M)
USE GlobalDataFun
IMPLICIT NONE
REAL(RDT) FIND_MIN
INTEGER,INTENT(OUT) :: M
INTEGER,INTENT(IN) :: N
REAL(RDT),INTENT(IN) :: MV(N)

REAL(RDT) TEMP
INTEGER I
    TEMP=MV(1)
    DO I=2,N
        IF(TEMP.GT.MV(I)) THEN
            M=I
            TEMP=MV(I)
        END IF
    END DO
    
    FIND_MIN=TEMP
END FUNCTION
