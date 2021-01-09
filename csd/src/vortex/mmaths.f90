module mmaths
USE GlobalDataFun
	implicit none
	contains
		!xl�ӳ�����֪��������������ab------------------------
		subroutine xl(a,b,c)
			implicit none
			real(rdt) a(3),b(3),c(3)
			c(:)=b(:)-a(:)
		end subroutine
		!cosx�ӳ�����֪����������н�����-----------------------
		subroutine cosx(a,b,c)
			implicit none 
			real(rdt) a(3),b(3),c,fen,fen1
			fen=SQRT((a(1)**2+a(2)**2+a(3)**2)*(b(1)**2+b(2)**2+b(3)**2))
			if(fen>0.0000000001)then		!��ֹĳ��������С���������
				c=(a(1)*b(1)+a(2)*b(2)+a(3)*b(3))/fen 
			else
                c=(a(1)*b(1)+a(2)*b(2)+a(3)*b(3))/(fen+1)
			end if
		end subroutine
		!wi�ӳ�����֪������a,b,������a*b-------------------------
		subroutine wi(a,b,c)
			implicit none
			real(rdt) a(3),b(3),c(3)
			c(1)=a(2)*b(3)-a(3)*b(2)
			c(2)=a(3)*b(1)-a(1)*b(3)
			c(3)=a(1)*b(2)-a(2)*b(1)
		end subroutine
		!fang�ӳ�����֪һ����,���䵥λ��������-----------------------
		subroutine fang(a,b)
			implicit none
			real(rdt) a(3),b(3),lab,lab1
			lab=SQRT(a(1)**2+a(2)**2+a(3)**2)
			if(lab>0.0000000001)then		!��ֹĳ��������С���������
				b=a/lab
			else
                b=a/(lab+1)
			end if
		end subroutine
		!length�ӳ�����֪���㣬����������-----------------------------
		subroutine length(a,b,l)
		    implicit none
		    real(rdt) a(3),b(3),l
		    l=SQRT((a(1)-b(1))**2+(a(2)-b(2))**2+(a(3)-b(3))**2)
		end subroutine

        !inverse�ӳ����������� 
        subroutine inverse(a,n,l)
        implicit none
        real(rdt)::a(n,n),t,d
        integer::n,l,i,j,k,is(n),js(n)
	    l=1
        do k=1,n
	        d=0.0
            do i=k,n
                do j=k,n
                    if (ABS(a(i,j))>d) then
                        d=ABS(a(i,j))
	                    is(k)=i
	                    js(k)=j
                    end if
                end do
            end do
            if(d+1.0==1.0) then
	            l=0
	            write(*,*) 'Error---The matrix is not inv'
	            return
            end if
            do j=1,n
	            t=a(k,j)
	            a(k,j)=a(is(k),j)
	            a(is(k),j)=t
            end do
            do i=1,n
	            t=a(i,k)
	            a(i,k)=a(i,js(k))
	            a(i,js(k))=t
            end do
	        a(k,k)=1/a(k,k)
            do j=1,n
                if (j/=k) then
	                a(k,j)=a(k,j)*a(k,k)
                end if
            end do
            do i=1,n
                if(i/=k) then
                    do j=1,n
                        if(j/=k) then
            	            a(i,j)=a(i,j)-a(i,k)*a(k,j)
                        end if
                    end do
                end if
            end do
            do i=1,n
                if (i/=k) then
                    a(i,k)=-a(i,k)*a(k,k)
                end if
            end do
        end do
        do k=n,1,-1
            do j=1,n
	            t=a(k,j)
	            a(k,j)=a(js(k),j)
	            a(js(k),j)=t
            end do
            do i=1,n
	            t=a(i,k)
	            a(i,k)=a(i,is(k))
	            a(i,is(k))=t
            end do
        end do
	    return
        end subroutine
        
end module
