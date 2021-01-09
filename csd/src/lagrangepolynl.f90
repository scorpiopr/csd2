!->+++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> the interface for functions.
module Mmyfunifs

	interface

            function lagbasd2twp2(n, x, xdt) result(xk)
                  use GlobalDataFun
                  integer, intent(in):: n
                  real(RDT), intent(in):: x, xdt
                  real(RDT) xk(n+1)
            end function
           
            
            function lagbasdtwp2(n, x, xdt) result(xk)
                  use GlobalDataFun
                  integer, intent(in):: n
                  real(RDT), intent(in):: x, xdt
                  real(RDT) xk(n+1)
            end function lagbasdtwp2 
           
            
            function lagbasdtwp(n, x, xdt) result(xk)
                  use GlobalDataFun
                  integer, intent(in):: n
                  real(RDT), intent(in):: x, xdt
                  real(RDT) xk(n+1)
            end function lagbasdtwp
            
            
            function lagbaswp(n, x) result(xk) 
                  use GlobalDataFun 
                  integer, intent(in):: n
                  real(RDT), intent(in):: x
                  real(RDT) xk(n+1)
            end function lagbaswp
            
            
            function lagbas(n, m, pp, x) result(xk) 
                  use GlobalDataFun 
                  integer, intent(in):: n, m
                  real(RDT), intent(in):: pp(n+1, m), x
                  real(RDT) xk(m)      
            end function lagbas


            function polydiff(n, p) result(xk) 
                  use GlobalDataFun
                  integer, intent(inout):: n
                  real(RDT), intent(in):: p(n+1)
                  real(RDT) xk(n)
            end function polydiff
            
      end interface
end module
!->+++++++++++++++++++++++++++++++++++++++++++++

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> the driver / wrapper.
!-> get the second derivative values of Lagrange basis polynomials 
!-> of degree n evaluated at x with xdt.
function lagbasd2twp2(n, x, xdt) result(xk)
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
use Mmyfunifs
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
integer, intent(in):: n
real(RDT), intent(in):: x, xdt
real(RDT) xk(n+1)

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->       lagbas, polydiff
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->       
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
integer i, j
integer m
real(RDT) pp3d2t(2, 4), pp4d2t(3, 5), pp5d2t(4, 6)
real(RDT) pp3dt(3, 4), pp4dt(4, 5), pp5dt(5, 6)
real(RDT) pp3(4, 4), pp4(5, 5), pp5(6, 6)
real(RDT) pp3x(4, 4), pp4x(5, 5), pp5x(6, 6)

data ((pp3(j, i), j=1, 4), i=1, 4) /&
      -9.0D0, 18.0D0, -11.0D0, 2.0D0, &
      27.0D0, -45.0D0, 18.0D0, 0.0D0, &
      -27.0D0, 36.0D0, -9.0D0, 0.0D0, &
      9.0D0, -9.0D0, 2.0D0, 0.0D0/
      
data ((pp4(j, i), j=1, 5), i=1, 5) /&
      32.0D0, -80.0D0, 70.0D0, -25.0D0, 3.0D0, &
      -128.0D0, 288.0D0, -208.0D0, 48.0D0, 0.0D0, &
      192.0D0, -384.0D0, 228.0D0, -36.0D0, 0.0D0, &
      -128.0D0, 224.0D0, -112.0D0, 16.0D0, 0.0D0, &
      32.0D0, -48.0D0, 22.0D0, -3.0D0, 0.0D0/

data ((pp5(j, i), j=1, 6), i=1, 6) /&
      -625.0D0, 1875.0D0, -2125.0D0, 1125.0D0, -274.0D0, 24.0D0, &
      3125.0D0, -8750.0D0, 8875.0D0, -3850.0D0, 600.0D0, 0.0D0, &
      -6250.0D0, 16250.0D0, -14750.0D0, 5350.0D0, -600.0D0, 0.0D0, &
      6250.0D0, -15000.0D0, 12250.0D0, -3900.0D0, 400.0D0, 0.0D0, &
      -3125.0D0, 6875.0D0, -5125.0D0, 1525.0D0, -150.0D0, 0.0D0, &
      625.0D0, -1250.0D0, 875.0D0, -250.0D0, 24.0D0, 0.0D0/

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

      pp3x=pp3/2.0D0
      pp4x=pp4/3.0D0
      pp5x=pp5/24.0D0

      lagord: select case(n)
      case(3)
            do 30 i=1, n+1
                  m=n 
                  pp3dt(:, i)=polydiff(m, pp3x(:, i))
                  m=n-1
                  pp3d2t(:, i)=polydiff(m, pp3dt(:, i))
30        continue

            xk=lagbas(n-2, n+1, pp3d2t, x)
      case(4)
            do 40 i=1, n+1
                  m=n 
                  pp4dt(:, i)=polydiff(m, pp4x(:, i))
                  m=n-1
                  pp4d2t(:, i)=polydiff(m, pp4dt(:, i))                  
40        continue      

            xk=lagbas(n-2, n+1, pp4d2t, x)
      case(5)
            do 50 i=1, n+1
                  m=n 
                  pp5dt(:, i)=polydiff(m, pp5x(:, i))
                  m=n-1
                  pp5d2t(:, i)=polydiff(m, pp5dt(:, i))         
50        continue     
 
            xk=lagbas(n-2, n+1, pp5d2t, x) 
      end select lagord
      
      !-> lagbasd2twp2=lagbasd2twp2*xdt*xdt+lagbasdtwp2*xd2t      
      !->///    xd2t=0.0D0 
      xk=xk*xdt*xdt
      
end function lagbasd2twp2
!->+++++++++++++++++++++++++++++++++++++++++++++



!->+++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> the driver / wrapper.
!-> get the [first] derivative values of Lagrange basis polynomials 
!-> of degree n evaluated at x with xdt.
function lagbasdtwp2(n, x, xdt) result(xk)
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
use Mmyfunifs
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
integer, intent(in):: n
real(RDT), intent(in):: x, xdt
real(RDT) xk(n+1)

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->       lagbas, polydiff
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->       
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
integer i, j
integer m
real(RDT) pp3dt(3, 4), pp4dt(4, 5), pp5dt(5, 6)
real(RDT) pp3(4, 4), pp4(5, 5), pp5(6, 6)
real(RDT) pp3x(4, 4), pp4x(5, 5), pp5x(6, 6)

data ((pp3(j, i), j=1, 4), i=1, 4) /&
      -9.0D0, 18.0D0, -11.0D0, 2.0D0, &
      27.0D0, -45.0D0, 18.0D0, 0.0D0, &
      -27.0D0, 36.0D0, -9.0D0, 0.0D0, &
      9.0D0, -9.0D0, 2.0D0, 0.0D0/
      
data ((pp4(j, i), j=1, 5), i=1, 5) /&
      32.0D0, -80.0D0, 70.0D0, -25.0D0, 3.0D0, &
      -128.0D0, 288.0D0, -208.0D0, 48.0D0, 0.0D0, &
      192.0D0, -384.0D0, 228.0D0, -36.0D0, 0.0D0, &
      -128.0D0, 224.0D0, -112.0D0, 16.0D0, 0.0D0, &
      32.0D0, -48.0D0, 22.0D0, -3.0D0, 0.0D0/

data ((pp5(j, i), j=1, 6), i=1, 6) /&
      -625.0D0, 1875.0D0, -2125.0D0, 1125.0D0, -274.0D0, 24.0D0, &
      3125.0D0, -8750.0D0, 8875.0D0, -3850.0D0, 600.0D0, 0.0D0, &
      -6250.0D0, 16250.0D0, -14750.0D0, 5350.0D0, -600.0D0, 0.0D0, &
      6250.0D0, -15000.0D0, 12250.0D0, -3900.0D0, 400.0D0, 0.0D0, &
      -3125.0D0, 6875.0D0, -5125.0D0, 1525.0D0, -150.0D0, 0.0D0, &
      625.0D0, -1250.0D0, 875.0D0, -250.0D0, 24.0D0, 0.0D0/

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

      pp3x=pp3/2.0D0
      pp4x=pp4/3.0D0
      pp5x=pp5/24.0D0

      lagord: select case(n)
      case(3)
            do 30 i=1, n+1
                  m=n 
                  pp3dt(:, i)=polydiff(m, pp3x(:, i))
30        continue

            xk=lagbas(n-1, n+1, pp3dt, x)
      case(4)
            do 40 i=1, n+1
                  m=n 
                  pp4dt(:, i)=polydiff(m, pp4x(:, i))
40        continue      

            xk=lagbas(n-1, n+1, pp4dt, x)
      case(5)
            do 50 i=1, n+1
                  m=n 
                  pp5dt(:, i)=polydiff(m, pp5x(:, i))
50        continue     
 
            xk=lagbas(n-1, n+1, pp5dt, x) 
      end select lagord
      
      xk=xk*xdt
      
end function lagbasdtwp2
!->+++++++++++++++++++++++++++++++++++++++++++++



!->+++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> the driver / wrapper.
!-> get the [first] derivative values of Lagrange basis polynomials 
!-> of degree n evaluated at x with xdt.
function lagbasdtwp(n, x, xdt) result(xk)
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
use Mmyfunifs
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
integer, intent(in):: n
real(RDT), intent(in):: x, xdt
real(RDT) xk(n+1)

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->       lagbas
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->       
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
integer i, j
real(RDT) pp3dt(4, 4), pp4dt(5, 5), pp5dt(6, 6)
real(RDT) pp3dtx(4, 4), pp4dtx(5, 5), pp5dtx(6, 6)

data ((pp3dt(j, i), j=1, 4), i=1, 4) /&
      0.0D0, -27.0D0, 36.0D0, -11.0D0, &
      0.0D0, 81.0D0, -90.0D0, 18.0D0, &
      0.0D0, -81.0D0, 72.0D0, -9.0D0, &
      0.0D0, 27.0D0, -18.0D0, 2.0D0/
      
data ((pp4dt(j, i), j=1, 5), i=1, 5) /&
      0.0D0, 128.0D0, -240.0D0, 140.0D0, -25.0D0, &
      0.0D0, -512.0D0, 864.0D0, -416.0D0, 48.0D0, &
      0.0D0, 768.0D0, -1152.0D0, 456.0D0, -36.0D0, &
      0.0D0, -512.0D0, 672.0D0, -224.0D0, 16.0D0, &
      0.0D0, 128.0D0, -144.0D0, 44.0D0, -3.0D0/

data ((pp5dt(j, i), j=1, 6), i=1, 6) /&
      0.0D0, -3125.0D0, 7500.0D0, -6375.0D0, 2250.0D0, -274.0D0, &
      0.0D0, 15625.0D0, -35000.0D0, 26625.0D0, -7700.0D0, 600.0D0, &
      0.0D0, -31250.0D0, 65000.0D0, -44250.0D0, 10700.0D0, -600.0D0, &
      0.0D0, 31250.0D0, -60000.0D0, 36750.0D0, -7800.0D0, 400.0D0, &
      0.0D0, -15625.0D0, 27500.0D0, -15375.0D0, 3050.0D0, -150.0D0, &
      0.d00, 3125.0D0, -5000.0D0, 2625.0D0, -500.0D0, 24.0D0/

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

      pp3dtx=pp3dt/2.0D0
      pp4dtx=pp4dt/3.0D0
      pp5dtx=pp5dt/24.0D0
      
      lagord: select case(n)
      case(3)
            xk=lagbas(n, n+1, pp3dtx, x)
      case(4)
            xk=lagbas(n, n+1, pp4dtx, x)
      case(5)
            xk=lagbas(n, n+1, pp5dtx, x) 
      end select lagord
      
      xk=xk*xdt
      
end function lagbasdtwp
!->+++++++++++++++++++++++++++++++++++++++++++++



!->+++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> the driver / wrapper.
!-> get the values of Lagrange basis polynomials of degree n evaluated at x.
function lagbaswp(n, x) result(xk) 
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
use Mmyfunifs
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
integer, intent(in):: n
real(RDT), intent(in):: x
real(RDT) xk(n+1)

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
!->       lagbas
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->       
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
integer i, j
real(RDT) pp3(4, 4), pp4(5, 5), pp5(6, 6)
real(RDT) pp3x(4, 4), pp4x(5, 5), pp5x(6, 6)

data ((pp3(j, i), j=1, 4), i=1, 4) /&
      -9.0D0, 18.0D0, -11.0D0, 2.0D0, &
      27.0D0, -45.0D0, 18.0D0, 0.0D0, &
      -27.0D0, 36.0D0, -9.0D0, 0.0D0, &
      9.0D0, -9.0D0, 2.0D0, 0.0D0/
      
data ((pp4(j, i), j=1, 5), i=1, 5) /&
      32.0D0, -80.0D0, 70.0D0, -25.0D0, 3.0D0, &
      -128.0D0, 288.0D0, -208.0D0, 48.0D0, 0.0D0, &
      192.0D0, -384.0D0, 228.0D0, -36.0D0, 0.0D0, &
      -128.0D0, 224.0D0, -112.0D0, 16.0D0, 0.0D0, &
      32.0D0, -48.0D0, 22.0D0, -3.0D0, 0.0D0/

data ((pp5(j, i), j=1, 6), i=1, 6) /&
      -625.0D0, 1875.0D0, -2125.0D0, 1125.0D0, -274.0D0, 24.0D0, &
      3125.0D0, -8750.0D0, 8875.0D0, -3850.0D0, 600.0D0, 0.0D0, &
      -6250.0D0, 16250.0D0, -14750.0D0, 5350.0D0, -600.0D0, 0.0D0, &
      6250.0D0, -15000.0D0, 12250.0D0, -3900.0D0, 400.0D0, 0.0D0, &
      -3125.0D0, 6875.0D0, -5125.0D0, 1525.0D0, -150.0D0, 0.0D0, &
      625.0D0, -1250.0D0, 875.0D0, -250.0D0, 24.0D0, 0.0D0/

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

      pp3x=pp3/2.0D0
      pp4x=pp4/3.0D0
      pp5x=pp5/24.0D0
      

      lagord: select case(n)
      case(3)
            xk=lagbas(n, n+1, pp3x, x)
      case(4)
            xk=lagbas(n, n+1, pp4x, x)
      case(5)
            xk=lagbas(n, n+1, pp5x, x) 
      end select lagord
      
end function lagbaswp
!->+++++++++++++++++++++++++++++++++++++++++++++



!->+++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> get the values of m Lagrange basis polynomials of degree n evaluated at x.
function lagbas(n, m, pp, x) result(xk) 
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
integer, intent(in):: n, m
real(RDT), intent(in):: pp(n+1, m), x
real(RDT) xk(m)

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> FUNCTIONS INVOKED:
!->
real(RDT), external:: polyval
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> SUBROUTINES INVOKED:
!->       
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> LOCAL VARIABLES:
integer i

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> main body:

      do 10 i=1, m
            xk(i)=polyval(n, pp(:, i), x)
10  continue      

end function lagbas
!->+++++++++++++++++++++++++++++++++++++++++++++



!->+++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> differentiate a polynomial of degree n.
function polydiff(n, p) result(xk) 
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
integer, intent(inout):: n
real(RDT), intent(in):: p(n+1)
real(RDT) xk(n)

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

      do 10 i=n, 1, -1
            xk(i)=(n+1-i)*p(i) 
10  continue      
      n=n-1 
      
end function polydiff
!->+++++++++++++++++++++++++++++++++++++++++++++



!->+++++++++++++++++++++++++++++++++++++++++++++
!-> the purpose:
!-> output the value of a polynomial of degree n evaluated at x.
function polyval(n, p, x)
!->+++++++++++++++++++++++++++++++++++++++++++++
!-> modules used:
use GlobalDataFun
implicit none

!->+++++++++++++++++++++++++++++++++++++++++++++
!-> dummy arguments / interface: 
integer, intent(in):: n
real(RDT), intent(in):: p(n+1), x
real(RDT) polyval

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

      polyval=p(1)
      do 10 i=1, n
            polyval=x*polyval+p(i+1) 
10  continue      

end function polyval
!->+++++++++++++++++++++++++++++++++++++++++++++