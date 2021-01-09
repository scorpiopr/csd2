    SUBROUTINE UU(VTU,ID,DM,XX,LL)
    use GlobalDataFun
    implicit none

    INTEGER,INTENT(IN)::ID,DM
    REAL(RDT),INTENT(IN)::XX,LL
    REAL(RDT),INTENT(OUT)::VTU(DM)

    !->+++++++++++++++++++++++++++++++++++++++++++++
    !-> FUNCTIONS INVOKED:
    !->
    REAL(RDT),EXTERNAL::PHC1,PHC2,PHC3,PHC4
    REAL(RDT),EXTERNAL::PHC1D1,PHC1D2,PHC1D3,PHC1D4
    REAL(RDT),EXTERNAL::PHC2D1,PHC2D2,PHC2D3,PHC2D4
    !->+++++++++++++++++++++++++++++++++++++++++++++

    IF( DM .NE. 4 ) THEN
        WRITE(*, *) 'ERROR, INPUT FOR SHAPE FUNCTIONS.'
        PAUSE 90001
    END IF

    GET_OOD:SELECT CASE(ID)
    CASE(0)
        VTU(1)=PHC1(XX,LL) 		   ! Another form can be used in which the coeffs of
        VTU(2)=PHC2(XX,LL) 		   ! interpolating polynomials may be treated as a matrix.
        VTU(3)=PHC3(XX,LL)
        VTU(4)=PHC4(XX,LL)		   !-> NAD==4
    CASE(1)
        VTU(1)=PHC1D1(XX,LL)
        VTU(2)=PHC1D2(XX,LL)
        VTU(3)=PHC1D3(XX,LL)
        VTU(4)=PHC1D4(XX,LL)
    CASE(2)
        VTU(1)=PHC2D1(XX,LL)
        VTU(2)=PHC2D2(XX,LL)
        VTU(3)=PHC2D3(XX,LL)
        VTU(4)=PHC2D4(XX,LL)
    END SELECT GET_OOD

    END SUBROUTINE

    SUBROUTINE VV(VTV,ID,DM,XX,LL)
    use GlobalDataFun
    implicit none

    INTEGER,INTENT(IN)::ID,DM
    REAL(RDT),INTENT(IN)::XX,LL
    REAL(RDT),INTENT(OUT)::VTV(DM)

    !->+++++++++++++++++++++++++++++++++++++++++++++
    !-> FUNCTIONS INVOKED:
    !->
    REAL(RDT),EXTERNAL::PHQ1,PHQ2,PHQ3
    REAL(RDT),EXTERNAL::PHQ1D1,PHQ1D2,PHQ1D3
    REAL(RDT),EXTERNAL::PHQ2D1,PHQ2D2,PHQ2D3
    !->+++++++++++++++++++++++++++++++++++++++++++++


    IF( DM .NE. 3 ) THEN
        WRITE(*, *) 'ERROR, INPUT FOR SHAPE FUNCTIONS.'
        PAUSE 90001
    END IF

    GET_OOD:SELECT CASE(ID)
    CASE(0)
        VTV(1)=PHQ1(XX,LL)
        VTV(2)=PHQ2(XX,LL)
        VTV(3)=PHQ3(XX,LL) 			   !-> NBD==3
    CASE(1)
        VTV(1)=PHQ1D1(XX,LL)
        VTV(2)=PHQ1D2(XX,LL)
        VTV(3)=PHQ1D3(XX,LL)
    CASE(2)
        VTV(1)=PHQ2D1(XX,LL)
        VTV(2)=PHQ2D2(XX,LL)
        VTV(3)=PHQ2D3(XX,LL)
    END SELECT GET_OOD

    END SUBROUTINE



    FUNCTION PHC1(XX,LL)
    use GlobalDataFun
    implicit none

    REAL(RDT), INTENT(IN):: XX,LL
    REAL(RDT) PHC1
    PHC1=1-3*XX**2+2*XX**3
    END FUNCTION


    FUNCTION PHC2(XX,LL)
    use GlobalDataFun
    implicit none

    REAL(RDT), INTENT(IN):: XX,LL
    REAL(RDT) PHC2
    PHC2=(XX-2*XX**2+XX**3)*LL
    END FUNCTION


    FUNCTION PHC3(XX,LL)
    use GlobalDataFun
    implicit none

    REAL(RDT), INTENT(IN):: XX,LL
    REAL(RDT) PHC3
    PHC3=3*XX**2-2*XX**3
    END FUNCTION


    FUNCTION PHC4(XX,LL)
    use GlobalDataFun
    implicit none

    REAL(RDT), INTENT(IN):: XX,LL
    REAL(RDT) PHC4
    PHC4=(-XX**2+XX**3)*LL
    END FUNCTION


    FUNCTION PHC1D1(XX,LL)
    use GlobalDataFun
    implicit none

    REAL(RDT), INTENT(IN):: XX,LL
    REAL(RDT) PHC1D1
    PHC1D1=(-6*XX+6*XX**2)/LL
    END FUNCTION


    FUNCTION PHC1D2(XX,LL)
    use GlobalDataFun
    implicit none

    REAL(RDT), INTENT(IN):: XX,LL
    REAL(RDT) PHC1D2
    PHC1D2=1-4*XX+3*XX**2
    END FUNCTION


    FUNCTION PHC1D3(XX,LL)
    use GlobalDataFun
    implicit none

    REAL(RDT), INTENT(IN):: XX,LL
    REAL(RDT) PHC1D3
    PHC1D3=(6*XX-6*XX**2)/LL
    END FUNCTION


    FUNCTION PHC1D4(XX,LL)
    use GlobalDataFun
    implicit none

    REAL(RDT), INTENT(IN):: XX,LL
    REAL(RDT) PHC1D4
    PHC1D4=-2*XX+3*XX**2
    END FUNCTION


    FUNCTION PHC2D1(XX,LL)
    use GlobalDataFun
    implicit none

    REAL(RDT), INTENT(IN):: XX,LL
    REAL(RDT) PHC2D1
    PHC2D1=(-6+12*XX)/LL/LL
    END FUNCTION


    FUNCTION PHC2D2(XX,LL)
    use GlobalDataFun
    implicit none

    REAL(RDT), INTENT(IN):: XX,LL
    REAL(RDT) PHC2D2
    PHC2D2=(-4+6*XX)/LL
    END FUNCTION


    FUNCTION PHC2D3(XX,LL)
    use GlobalDataFun
    implicit none

    REAL(RDT), INTENT(IN):: XX,LL
    REAL(RDT) PHC2D3
    PHC2D3=(6-12*XX)/LL/LL
    END FUNCTION


    FUNCTION PHC2D4(XX,LL)
    use GlobalDataFun
    implicit none

    REAL(RDT), INTENT(IN):: XX,LL
    REAL(RDT) PHC2D4
    PHC2D4=(-2+6*XX)/LL
    END FUNCTION



    FUNCTION PHQ1(XX,LL)
    use GlobalDataFun
    implicit none

    REAL(RDT), INTENT(IN):: XX,LL
    REAL(RDT) PHQ1
    PHQ1=1-3*XX+2*XX**2
    END FUNCTION


    FUNCTION PHQ2(XX,LL)
    use GlobalDataFun
    implicit none

    REAL(RDT), INTENT(IN):: XX,LL
    REAL(RDT) PHQ2
    PHQ2=4*XX-4*XX**2
    END FUNCTION


    FUNCTION PHQ3(XX,LL)
    use GlobalDataFun
    implicit none

    REAL(RDT), INTENT(IN):: XX,LL
    REAL(RDT) PHQ3
    PHQ3=-XX+2*XX**2
    END FUNCTION


    FUNCTION PHQ1D1(XX,LL)
    use GlobalDataFun
    implicit none

    REAL(RDT), INTENT(IN):: XX,LL
    REAL(RDT) PHQ1D1
    PHQ1D1=(-3+4*XX)/LL
    END FUNCTION


    FUNCTION PHQ1D2(XX,LL)
    use GlobalDataFun
    implicit none

    REAL(RDT), INTENT(IN):: XX,LL
    REAL(RDT) PHQ1D2
    PHQ1D2=(4-8*XX)/LL
    END FUNCTION


    FUNCTION PHQ1D3(XX,LL)
    use GlobalDataFun
    implicit none

    REAL(RDT), INTENT(IN):: XX,LL
    REAL(RDT) PHQ1D3
    PHQ1D3=(-1+4*XX)/LL
    END FUNCTION


    FUNCTION PHQ2D1(XX,LL)
    use GlobalDataFun
    implicit none

    REAL(RDT), INTENT(IN):: XX,LL
    REAL(RDT) PHQ2D1
    PHQ2D1=4/LL/LL
    END FUNCTION


    FUNCTION PHQ2D2(XX,LL)
    use GlobalDataFun
    implicit none

    REAL(RDT), INTENT(IN):: XX,LL
    REAL(RDT) PHQ2D2
    PHQ2D2=-8/LL/LL
    END FUNCTION


    FUNCTION PHQ2D3(XX,LL)
    use GlobalDataFun
    implicit none

    REAL(RDT), INTENT(IN):: XX,LL
    REAL(RDT) PHQ2D3
    PHQ2D3=4/LL/LL
    END FUNCTION