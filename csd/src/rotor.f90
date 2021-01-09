!MODULE ROTOR_CLASS
!USE EIGEN_CLASS
!IMPLICIT NONE
!!->+++++++++++++++++++++++++++++++++++++++++++++
!TYPE,PUBLIC :: ROTOR
!    
!    INTEGER ::SOL
!    REAL(RDT) :: R
!    REAL(RDT) :: TH0,TH1C,TH1S,PSI,OMG
!    INTEGER :: NRBD
!
!    TYPE(EIGEN),ALLOCATABLE :: BEAMS(:)
!
!    CONTAINS
!    PROCEDURE,PUBLIC :: SOLCASE
!    PROCEDURE,PUBLIC :: CONSTRUCT_ROTOR
!
!END TYPE ROTOR
!!->++++++++++++++++++++++++++++++++++++++++++++
!CONTAINS
!!->+++++++++++++++++++++++++++++++++++++++++++++
!!->+++++++++++++++++++++++++++++++++++++++++++++
!!->+++++++++++++++PUBLIC++++++++++++++++++++++++++
!!->+++++++++++++++++++++++++++++++++++++++++++++
!SUBROUTINE SOLCASE(THIS,INPU)
!USE INPUT_CLASS
!IMPLICIT NONE
!CLASS(ROTOR) :: THIS
!CLASS(INPUT) :: INPU
!INTEGER :: I,j,k,imk
!    WRITE(*,*) '+++++++++++++++++++++++++++++++++++++++++++++'
!    WRITE(*,*) '构建计算模型'
!    CALL THIS%CONSTRUCT_ROTOR(INPU)
!!->++++++++++++++++++++++++++++++++++++++++++++
!    WRITE(*,*) '形成单元刚度阵'
!    THIS%TH0=0.0D0
!    THIS%TH1C=0.0D0
!    THIS%TH1S=0.0D0
!    THIS%PSI=0.0D0
!
!    DO I=1,THIS%NRBD
!        CALL THIS%BEAMS(I)%MLEGGAUSS_BEAM(THIS%R,&
!                                                                            THIS%TH0,&
!                                                                            THIS%TH1C,&
!                                                                            THIS%TH1S,&
!                                                                            THIS%PSI,&
!                                                                            THIS%OMG,&
!                                                                            THIS%SOL)
!    END DO
!!->++++++++++++++++++++++++++++++++++++++++++++
!    WRITE(*,*) '形成总刚度阵'
!    DO I=1,THIS%NRBD
!        CALL THIS%BEAMS(I)%FORMULATE()
!    END DO
!
!    DO I=1,THIS%NRBD
!        CALL THIS%BEAMS(I)%GIVECONS()
!!        open(11,file='tk.dat')
!!        open(12,file='tm.dat') 
!!        open(13,file='tc.dat') 
!!        open(14,file='tf.dat') 
!!        do k=1,THIS%BEAMS(I)%NDOFC
!!            write(11,995) (THIS%BEAMS(I)%UKMC(k,j),j=1,THIS%BEAMS(I)%NDOFC)
!!            write(12,995) (THIS%BEAMS(I)%UMMC(k,j),j=1,THIS%BEAMS(I)%NDOFC)
!!            write(13,995) (THIS%BEAMS(I)%UCMC(k,j),j=1,THIS%BEAMS(I)%NDOFC)
!!            write(14,995) THIS%BEAMS(I)%UFVC(k)
!!        end do
!!
!!        close(11)  
!!        close(12) 
!!        close(13)
!!        close(14)
!
!    END DO
!
!    DO I=1,THIS%NRBD
!        CALL THIS%BEAMS(I)%CONSTRUCT_EIGEN(INPU)
!        CALL THIS%BEAMS(I)%eigensolv(imk)
!        WRITE(*,*) THIS%BEAMS(I)%VEGF/THIS%OMG
!    END DO
!
!
!995 format(400f26.12)
!!->++++++++++++++++++++++++++++++++++++++++++++
!END SUBROUTINE
!!->++++++++++++++++++++++++++++++++++++++++++++
!SUBROUTINE CONSTRUCT_ROTOR(THIS,INPU)
!USE INPUT_CLASS
!IMPLICIT NONE
!CLASS(ROTOR) :: THIS
!CLASS(INPUT) :: INPU
!
!INTEGER :: I
!
!    THIS%SOL=INPU%IPT%SOL
!    THIS%R=INPU%MAT%R
!    THIS%TH0=INPU%IPT%TH0
!    THIS%TH1C=INPU%IPT%TH1C
!    THIS%TH1S=INPU%IPT%TH1S
!    THIS%OMG=INPU%MAT%OMG
!
!    IF(THIS%SOL.EQ.0) THEN
!        THIS%NRBD=1
!    ELSE
!        THIS%NRBD=INPU%IPT%NRBD
!    END IF
!
!    ALLOCATE(THIS%BEAMS(THIS%NRBD))
!
!    DO I=1,THIS%NRBD
!        CALL THIS%BEAMS(I)%CONSTRUCT_BEAM(INPU)
!    END DO
!
!END SUBROUTINE
!
!END MODULE