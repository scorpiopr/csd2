module get_vi_y
    use LinearAlgebra
    !use global
    !use cl_cd
    !use vi_hover
    !use mkl95_precision
    !use mkl95_blas
    !use blas95
    implicit none
    
    REAL(8),parameter :: pi=3.141592653
    !REAL(8) :: dpsi,rr0,r_cut=0.1117,rr=5.4864,dr,omega=33.73,vx,vz,gamma=0.0,phi70=0.0,dphi0=0.0,rho=1.121,cc,aa,hh=0.762,mu=0.3963,alpha_s=7.5!共轴
    !INTEGER :: kk=3,np=180,nr=50!共轴
    REAL(8) :: dpsi,rr0,r_cut,rr,dr,omega,vx,vz,gamma,phi70,dphi0,rho,cc,aa,hh,mu,alpha_s!gamma：提前操纵角；phi70：0.7R处安装角；dphi0：桨叶扭度
    INTEGER :: kk,np,nr
    
    PRIVATE pi,dpsi,rr0,r_cut,rr,dr,omega,vx,vz,gamma,phi70,dphi0,rho,cc,aa,hh,mu,alpha_s,kk,np,nr
    
    CONTAINS
    
    SUBROUTINE COAXIALINFLOWMAIN(viu,vil,wi_bar)
        IMPLICIT NONE
        REAL(8) :: THETA(6)!上下旋翼总距和周期变距
        REAL(8) :: viu(:,:),vil(:,:),wi_bar
        
        CALL READCOAXIALTRIM(THETA)
        CALL get_vi_forward(THETA,viu,vil,wi_bar)
    
    END SUBROUTINE
    
    SUBROUTINE READCOAXIALTRIM(THETA)
        IMPLICIT NONE
        REAL(8) :: THETA(:)

        OPEN(151,FILE=trim('COAXIALTRIMINPUT.cfg') )
        NAMELIST /COAXIALTRIM/ THETA,r_cut,rr,omega,gamma,phi70,dphi0,rho,hh,mu,alpha_s,kk,np,nr

        REWIND(151)
        READ(151,NML=COAXIALTRIM)
        !WRITE(*,NML=COAXIALTRIM)
        
    END SUBROUTINE
   
    function get_ad(wd)
        implicit none
        real(8)::get_ad,wd(nr,np)
        integer::i,j
        get_ad=0.0
        do i=1,nr,1
            do j=1,np,1
                if(wd(i,j)/=0.0) then
                    get_ad=get_ad+pi*((rr0+(rr-rr0)*i/nr)**2-(rr0+(rr-rr0)*(i-1)/nr)**2)/np
                end if
            end do
        end do
        return
    end function
    subroutine flag_xy(xs,ys,ii,jj,flag_rs)
        implicit none
        integer::flag_rs,ii,jj
        real(8)::xs,ys
        real(8)::rs,psi
        flag_rs=0
        rs=sqrt(xs**2+ys**2)-dr/2.0            
        if((rs>=rr0).and.(rs<rr)) then
            ii=int((rs-rr0)/dr)+1
        else if (rs==rr) then           !用了等号来判断
            ii=nr
        else 
            return
        end if
        !
        if(xs==0.0) then
            if(ys>0) then
                psi=0.5*180.0
            else
                psi=1.5*180.0
            end if
        else    
            if((xs>0.0).and.(ys>0.0)) then
                psi=atand(ys/xs)
            else if(xs<0.0)then
                psi=180.0+atand(ys/xs)
            else if((xs>0.0).and.(ys<=0.0)) then
                psi=2*180.0+atand(ys/xs)
            end if
        end if
        jj=int(psi/dpsi)+1
        if(jj>=np)then
            jj=np
        end if
        flag_rs=1
        !write(*,*)'flag_rs',flag_rs,ii,jj
        return
    end subroutine flag_xy
    subroutine get_epsilon(ad,bw,ep,chi)    !ad,bw->ep,chi
        implicit none
        real(8)::ad(2),bw(2),ep(2),chi(2)
        real(8)::chi2(2),dchi(2),chi_f(2),y(2),y2(2),y_f(2)
        real(8)::vn(2)
        real(8)::res_chi(2),residual_chi,res_dchi
        integer::i_epsilon
        !赋初值
        y_f(:)=(/0.0,0.0/)
        y(:)=(/0.0,0.0/);y2(:)=(/0.0,0.0/)
        chi(:)=(/0.0,0.0/);chi2(:)=(/0.0,0.0/)
        vn(:)=(/0.0,0.0/)
        !**********************************************
        dchi=(/10.0,10.0/)
        res_dchi=10.0
        residual_chi=1.0
        i_epsilon=0
        do while((residual_chi>0.1).and.(abs(res_dchi)>0.1))
            i_epsilon=i_epsilon+1
            ep(1)=(1+(-hh)/sqrt(rr**2*cosd(chi(1))**2+hh**2))
            ep(2)=(1+(hh)/sqrt(rr**2*cosd(chi(2))**2+hh**2))
            vn(1)=vz+bw(1)+(ad(1)/aa)*ep(1)*bw(2)
            vn(2)=vz+bw(2)+(ad(2)/aa)*ep(2)*bw(1)
            y(1)=chi(1)-atand(vx/vn(1))
            y(2)=chi(2)-atand(vx/vn(2))
            !
            chi2=chi+dchi
            !
            ep(1)=(1+(-hh)/sqrt(rr**2*cosd(chi2(1))**2+hh**2))
            ep(2)=(1+(hh)/sqrt(rr**2*cosd(chi2(2))**2+hh**2))
            vn(1)=vz+bw(1)+(ad(1)/aa)*ep(1)*bw(2)
            vn(2)=vz+bw(2)+(ad(2)/aa)*ep(2)*bw(1)
            y2(1)=chi2(1)-atand(vx/vn(1))
            y2(2)=chi2(2)-atand(vx/vn(2))
            !
            res_chi=abs(y_f-y2)
            residual_chi=max(res_chi(1),res_chi(2))
            !
            dchi(1)=(y_f(1)-y2(1))*(chi2(1)-chi(1))/(y2(1)-y(1))
            dchi(2)=(y_f(2)-y2(2))*(chi2(2)-chi(2))/(y2(2)-y(2))
            res_dchi=max(abs(dchi(1)),abs(dchi(2)))
            !write(*,*)'dchi',dchi
            chi=chi2
            !write(*,*)'res get_epsilon',i_epsilon,residual_chi
            !write(43,*)'res get_epsilon',i_epsilon,residual_chi
            if((i_epsilon>5).and.(residual_chi<0.5))then
                residual_chi=0.0
            else if((i_epsilon>8).and.(residual_chi<1.0))then
                residual_chi=0.0
            else if(i_epsilon>10)then
                residual_chi=0.0
            end if
        end do
        chi=chi2
        ep(1)=(1+(-hh)/sqrt(rr**2*cosd(chi2(1))**2+hh**2))
        ep(2)=(1+(hh)/sqrt(rr**2*cosd(chi2(2))**2+hh**2))
        return
    end subroutine get_epsilon
    subroutine bw_w(bw,wu,wl,wud,wld,chi)    !wud,wld,bw->epsilon,chi->wu,wl
        implicit none
        real(8)::bw(2),wu(nr,np),wl(nr,np),wud(nr,np),wld(nr,np)
        real(8)::ad(2),chi(2),kapx(2),kapy(2)
        integer::i,j
        !ad(1)=get_ad(wld)
        !ad(2)=get_ad(wud)
        if(mu>0.00001)then
            kapx(1)=(4.0/3.0)*(1-cosd(chi(1))-1.8*mu**2)/(sind(chi(1)))
            kapy(1)=-2.0*mu
            kapx(2)=(4.0/3.0)*(1-cosd(chi(2))-1.8*mu**2)/(sind(chi(2)))
            kapy(2)=2.0*mu
        else
            write(*,*)'前进比约为0，请选用悬停模块计算诱导速度！'
            write(*,*)'mu:',mu
            pause
            stop
        end if
        do i=1,nr
            do j=1,np
                wu(i,j)=bw(1)*(1.0+kapx(1)*(rr0+i*dr)/rr*cosd(j*dpsi)+kapy(1)*(rr0+i*dr)/rr*sind(j*dpsi))
                wl(i,j)=bw(2)*(1.0+kapx(2)*(rr0+i*dr)/rr*cosd(j*dpsi)+kapy(2)*(rr0+i*dr)/rr*sind(j*dpsi))
            end do
        end do
    end subroutine bw_w
    subroutine bw_wd(bw,wu,wl,wud,wld,wud2,wld2,ep)
        implicit none
        real(8)::bw(2),wu(nr,np),wl(nr,np),wud(nr,np),wld(nr,np),wud2(nr,np),wld2(nr,np),ep(2)
        real(8)::ad(2),r,psi,x,y,xs,ys,vn(2)
        real(8)::gam(2)
        integer::flag_rs,ii,jj
        integer::i,j
        !
        do i=1,nr
            do j=1,np
                wud2(i,j)=0.0;wld2(i,j)=0.0
            end do
        end do
        !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        if(mu>4.0)then
            do i=1,nr
                do j=1,np
                    wud2(i,j)=0.0;wld2(i,j)=0.0
                end do
            end do
            write(*,*)'由于前进比过大，直接认为上下旋翼间无干扰'
            return
        end if
        !
        ad(1)=get_ad(wld)
        ad(2)=get_ad(wud)
        !write(40,*)'ad',ad
        vn(1)=vz+bw(1)+(ad(1)/aa)*ep(1)*bw(2)
        vn(2)=vz+bw(2)+(ad(2)/aa)*ep(2)*bw(1)
        gam(1)=sqrt(vn(1)/(vz+ep(2)*bw(1)+(ad(1)/aa)*bw(2)))    
        gam(2)=sqrt(vn(2)/(vz+ep(1)*bw(2)+(ad(2)/aa)*bw(1)))
        do i=1,nr
            do j=1,np
                r=rr0+i*dr
                psi=j*dpsi
                x=r*cosd(psi)
                y=r*sind(psi)
                if(mu>5.0)then
                    wld2(i,j)=0.0;wud2(i,j)=0.0
                else
                    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                    
                    xs=gam(1)*x+(hh)*(vx/vn(1))
                    ys=gam(1)*y
                    call flag_xy(xs,ys,ii,jj,flag_rs)                   
                    if(flag_rs==1)then
                        wld2(i,j)=ep(1)*wl(ii,jj)   

                    end if
                   
                    xs=gam(2)*x+(-hh)*(vx/vn(2))
                    ys=gam(2)*y
                    call flag_xy(xs,ys,ii,jj,flag_rs)                   
                    if(flag_rs==1)then
                        wud2(i,j)=ep(2)*wu(ii,jj)
                        !wud2(i,j)=ep(2)*(wu(ii,jj)+wld(ii,jj))
                    end if
                    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                end if
            end do
        end do
        ad(1)=get_ad(wld2)
        ad(2)=get_ad(wud2)
        !write(40,*)'ad2',ad
        return
    end subroutine bw_wd
    subroutine bw_vi(bw,wu,wl,wud,wld)
        implicit none
        real(8)::bw(2),wu(nr,np),wl(nr,np),wud(nr,np),wld(nr,np)
        real(8)::wu0(nr,np),wl0(nr,np),wud0(nr,np),wld0(nr,np)
        real(8)::wu2(nr,np),wl2(nr,np),wud2(nr,np),wld2(nr,np)
        real(8)::ad(2),ad2(2),ep(2),chi(2),res_w,residual_w,residual_ad
        real(8),parameter::ad_ep=0.0005,vi_ep=0.001
        integer::i,j,i_bw
        !write(*,*)'bw_vi',bw

        do i=1,nr
            do j=1,np
                wu(i,j)=0.0;wl(i,j)=0.0;wud(i,j)=2.0;wld(i,j)=2.0
                wu0(i,j)=0.0;wl0(i,j)=0.0;wud0(i,j)=0.0;wld0(i,j)=0.0
            end do
        end do
        !write(*,*)'bw_vi',i,j
        i_bw=0
        residual_w=1.0
        res_w=1.0
        residual_ad=1.0
        !
        do while((residual_w>vi_ep).and.(residual_ad>ad_ep))
            i_bw=i_bw+1
            wu0=wu;wl0=wl;wud0=wud;wld0=wld 
            ad(1)=get_ad(wld)
            ad(2)=get_ad(wud)
            !write(*,*)'ad        ',ad
            call get_epsilon(ad,bw,ep,chi)
            call bw_w(bw,wu,wl,wud,wld,chi)
            call bw_wd(bw,wu,wl,wud,wld,wud2,wld2,ep)
            ad2(1)=get_ad(wld2)
            ad2(2)=get_ad(wud2)
            residual_ad=max(abs(ad2(1)-ad(1)),abs(ad2(2)-ad(2)))
            !write(*,*)'res_bw_vi ad',i_bw,residual_ad
            wud=wud2
            wld=wld2
            residual_w=max(abs(wu0(1,1)-wu(1,1)),abs(wl0(1,1)-wl(1,1)))
            do i=1,nr
                do j=1,np
                    res_w=max(abs(wu0(i,j)-wu(i,j)),abs(wl0(i,j)-wl(i,j)))
                    if(res_w>residual_w) then
                        residual_w=res_w
                    end if
                end do
            end do
            if((i_bw>3).and.(residual_ad<0.005))then
                residual_w=0.0
            else if((i_bw>6).and.(residual_w<0.01))then
                residual_w=0.0
            else if((i_bw>10).and.(residual_w<0.05))then
                residual_w=0.0
            end if
            !write(*,*)'res bw_vi',i_bw,residual_w
        end do!do while 循环结束
        !write(*,*)'res bw_vi',i_bw,residual_w
        !write(43,*)'res bw_vi',i_bw,residual_w
        return
    end subroutine bw_vi
    function fun_w(bw,theta)
        implicit none
        real(8)::fun_w(2),bw(2),theta(6)
        real(8)::ad(2)
        real(8),allocatable::wu(:,:),wl(:,:),wud(:,:),wld(:,:)
  
        real(8)::psi,r_j,r_r
        real(8)::wx,wy,w_xy,dtheta,theta_w,betax,alphax,cy
        real(8)::ep(2),chi(2),momentum(2),factor_r
        integer::i,j
  
        allocate(wu(nr,np),wl(nr,np),wud(nr,np),wld(nr,np))
        call bw_vi(bw,wu,wl,wud,wld)
        ad(1)=get_ad(wld)
        ad(2)=get_ad(wud)
        fun_w=(/0.0,0.0/)
        momentum=(/0.0,0.0/)
        wx=0.0;wy=0.0;w_xy=0.0
        !write(*,*)'     start fun_w!'
        do i=1,nr
           do j=1,np
                psi=j*dpsi
                r_j=rr0+i*dr
                !r_r=r_j/rr
 
                wx=omega*(r_j)+(1.0)*vx*sind(psi)
                wy=vz+wu(I,J)+wld(I,J)
                w_xy=sqrt(wx**2+wy**2)
                dtheta=theta(1)+theta(2)*cosd(psi+gamma)+theta(3)*sind(psi+gamma)!修改了总距计算公式
                theta_w=phi70+dphi0*(r_j/rr-0.7)+dtheta
                if ((abs(wy)<0.0001).and.(abs(wx)<0.0001))then
                    betax=0.0
                else
                    betax=atand(wy/wx)
                end if
                alphax=theta_w-betax
                !if (wx>0) then
                !    cy=RCl(alphax,w_xy)     !????????????????????????????????
                !else
                !    cy=RCl(180.0-abs(alphax),w_xy)
                !end if

                !if(r_r>0.97) then
                !    cy=0.0
                !end if
                
                cc=2*(-0.15*(r_j)/rr+0.29)
                !fun_w(1)=fun_w(1)+(kk*rho/4.0/pi)*cc*(cy/cosd(betax))*(w_xy**2)*(2.0*pi/np)*(dr)
                fun_w(1)=fun_w(1)+(kk*rho/4.0/pi)*cc*alphax*(w_xy**2)*(2.0*pi/np)*(dr)!改用简单的（升力线斜率alcsrf*翼型迎角）方式计算截面升力系数
                psi=-psi        !!!!!!!!!!
                wx=omega*(r_j)+(-1.0)*vx*sind(psi)
                wy=vz+wl(I,J)+wud(I,J)
                w_xy=sqrt(wx**2+wy**2)
                dtheta=theta(4)+theta(5)*cosd(psi+gamma)+theta(6)*sind(psi+gamma)
                theta_w=phi70+dphi0*(r_j/rr-0.7)+dtheta
                if ((abs(wy)<0.0001).and.(abs(wx)<0.0001))then
                    betax=0.0
                else
                    betax=atand(wy/wx)
                end if
                alphax=theta_w-betax
                !if (wx>0) then
                !    cy=RCl(alphax,w_xy)
                !else
                !    cy=RCl(180.0-abs(alphax),w_xy)
                !end if
                !
                !if(r_r>0.97) then
                !    cy=0.0
                !end if
                !
                !fun_w(2)=fun_w(2)+(kk*rho/4.0/pi)*cc*(cy/cosd(betax))*(w_xy**2)*(2.0*pi/np)*(dr)
                fun_w(2)=fun_w(2)+(kk*rho/4.0/pi)*cc*alphax*(w_xy**2)*(2.0*pi/np)*(dr)!改用简单的（升力线斜率alcsrf*翼型迎角）方式计算截面升力系数
            end do
        end do
        !write(40,*)'ct fun',fun_w/(0.5*rho*pi*rr**2*omega_rr**2)
        call get_epsilon(ad,bw,ep,chi)
        !write(*,*)'     bw',bw
        momentum(1)=(2*rho*ad(1)*sqrt(vx**2+(vz+bw(1)+ep(1)*bw(2))**2)+2*rho*(aa-ad(1))*sqrt(vx**2+(vz+bw(1))**2))*bw(1)
        momentum(2)=(2*rho*ad(2)*sqrt(vx**2+(vz+bw(2)+ep(2)*bw(1))**2)+2*rho*(aa-ad(2))*sqrt(vx**2+(vz+bw(2))**2))*bw(2)
        !fun_w=fun_w-momentum
        !write(*,*)'     momentum',momentum
        !write(*,*)'     fun_w',fun_w
        fun_w(1)=fun_w(1)-momentum(1)
        fun_w(2)=fun_w(2)-momentum(2)
    end function fun_w
    
    
    subroutine get_vi_forward(theta,viu,vil,wi_bar)
        implicit none
        real(8)::theta(6),viu(nr,np),vil(nr,np),wi_bar
        !中间变量
        real(8)::wu(nr,np),wl(nr,np),wud(nr,np),wld(nr,np)

        real(8)::ad_final(2),bw_final(2)      !最终计算值
        real(8)::ep(2),chi(2),epsilon_ul(2)
        real(8)::y_f_goal(2)            !目标量（保持不变）
        !通用函数变量
        real(8)::x0(2),x1(2),dx(2),y0(2),y1(2)
        real(8)::jac(2,2)
        real(8)::residual_y(2),residual_vi,res_dx
        real(8),parameter::relax_vi=1.0,fun_ep=0.5,dx_ep=1E-6
        real(8)::ijac(2,2),xmkl(2),ymkl(2)
        integer::i,j,k,m,i_vi
        !
        !write(*,*)'Start get_vi_forward!'
        
        
               
        rr0=r_cut*rr
        dpsi=360.0/np
        dr=(rr-rr0)/nr
        aa=pi*(rr**2-rr0**2)
        vx=mu*omega*rr
        vz=mu*tand(alpha_s)*omega*rr
        
        ad_final=(/0.0,0.0/)
        x0=(/6.0,6.0/)!bw
        y_f_goal=(/30000,30000/)
        dx=(/1.0,1.0/)!d_bw
        i_vi=0
        residual_vi=10.0
        res_dx=10.0
        !do while((abs(residual_vi)>fun_ep).and.(abs(res_dx)>dx_ep))!fun_w中是没无量纲化的拉力差，因此收敛值可以比较大
        !    i_vi=i_vi+1
        !    do j=1,2,1
        !        do m=1,2,1
        !            if(m==j)then
        !                x1(m)=x0(m)+dx(m)
        !            else
        !                x1(m)=x0(m)
        !            end if
        !        end do
        !        !write(*,*)'fun_w'
        !        !write(*,*)x0
        !        !write(*,*)x1 
        !        !write(*,*)dx(j)
        !        y0=fun_w(x0,theta)
        !        y1=fun_w(x1,theta)             !??????????????????????????????
        !        do k=1,2,1
        !            jac(k,j)=(y1(k)-y0(k))/(dx(j))
        !        end do
        !    end do
        !    x1=x0+relax_vi*dx
        !    y1=fun_w(x1,theta)
        !    residual_y=y_f_goal-y1
        !    residual_y=abs(residual_y)
        !    residual_vi=max(residual_y(1),residual_y(2))
        !    call inverse(jac,ijac)
        !    ymkl=y_f_goal-y1
        !    call gemv(ijac,ymkl,xmkl)
        !    dx=xmkl
        !    res_dx=max(abs(dx(1)),abs(dx(2)))
        !    write(*,*)'res get_vi_forward',i_vi,residual_vi
        !    write(43,*)'res get_vi_forward',i_vi,residual_vi
        !    x0=x1
        !end do
        !write(*,*)'res get_vi_forward',i_vi,residual_vi
        !write(43,*)'res get_vi_forward',i_vi,residual_vi
        !!**********************************************************
        !bw_final=x1
        bw_final(1)=wi_bar
        bw_final(2)=wi_bar
        call bw_vi(bw_final,wu,wl,wud,wld)
        viu=wu+wld
        vil=wl+wud
        !write(30,*)'psi r   wu  wl  wud wld viu vil'
        !do i=1,np
        !    do j=1,nr
        !        write(30,"(4X,F14.9,\)")i*dpsi,(rr0+j*dr)/rr,wu(i,j),wl(i,j),wud(i,j),wld(i,j),wu(i,j)+wld(i,j),wl(i,j)+wud(i,j)
        !        write(30,*)'    '
        !    end do
        !end do
        ad_final(1)=get_ad(wld)
        ad_final(2)=get_ad(wud)
        call get_epsilon(ad_final,bw_final,ep,chi)
        epsilon_ul(1)=ep(2)*ad_final(2)/aa  
        epsilon_ul(2)=ep(1)*ad_final(1)/aa 
        !write(*,*)'ad_final*****',ad_final
        !write(*,*)'bw_final*****',x1
        !write(*,*)'mu  adl epu_l   adu epl_u'
        !!write(*,"(2X,F10.6,\)")mu,ad_final(2)/aa,epsilon_ul(1),ad_final(1)/aa,epsilon_ul(2)
        !write(*,"(2X,F10.6,\)")hh/rr,ad_final(2)/aa,epsilon_ul(1),ad_final(1)/aa,epsilon_ul(2)
        !write(*,*)'    '
        !!write(40,"(5X,F10.6,\)")mu,ad_final(2)/aa,epsilon_ul(1),ad_final(1)/aa,epsilon_ul(2)
        !write(40,"(5X,F10.6,\)")hh/rr,ad_final(2)/aa,epsilon_ul(1),ad_final(1)/aa,epsilon_ul(2)
        !write(40,*)'    '
        return
    end subroutine get_vi_forward
  !  function get_y_mid(theta_f) result(y_f)
  !
  !      implicit none
  !      real(8)::theta_f(6),y_f(6)  
  !      !real(8)::y_fun1(6),y_fun2(6),y_fun3(6)
  !      real(8)::y_u(6),y_l(6)
  !      real(8)::vi,vi2,viu(nr,np),vil(nr,np),viu_hover(nr),vil_hover(nr),omegaU,omegaL                                                !诱导速度
  !      real(8)::wu(nr,np),wl(nr,np),wud(nr,np),wld(nr,np)
  !      real(8)::phi,alphax,betax,phi2,alphax2,betax2        
  !      real(8)::dtheta,dtheta2                                  
  !      real(8)::w_xy,wx,wy,w_xy2,wx2,wy2
  !   
  !      real(8)::r_r,r_j,psi,psi2
  !      real(8)::dyy(6),yy(6)                                
  !      real(8)::dyy2(6),yy2(6)
  !      real(8)::cx,cy,dy,dx,cx2,cy2,dy2,dx2                        !翼型升、阻力系数，课本中间变量dY,dX
  !      real(8)::dt,dq,dt2,dq2,factor_r(2)
  !      integer::i,j,k   
  !      omegaU=omega
  !      omegaL=omega
  !      y_u=(/0.0,0.0,0.0,0.0,0.0,0.0/)
  !      y_l=(/0.0,0.0,0.0,0.0,0.0,0.0/)
  !      !初始化外层积分变量
  !      psi=0.0!初始化方位角积分变量
  !      psi2=0.0
  !      yy(:)=(/0.0,0.0,0.0,0.0,0.0,0.0/)
  !      yy2(:)=(/0.0,0.0,0.0,0.0,0.0,0.0/)
  !      !求解诱导速度
  !      if(mu<0.0)then
  !          write(*,*)'mu为负值，请检查输入的前进比！'
  !          pause
  !          stop
  !      end if
  !      if(abs(mu>0.000001))then
  !          call get_vi_forward(theta_f,viu,vil)            !??????????????????????????????????
  !      else
  !          call get_vi_hover(theta_f,viu_hover,vil_hover)
  !      end if
  !
  !      do i=1,np,1
  !          do j=1,nr,1
  !              psi=i*dpsi                                   
  !              psi2=i*dpsi
  !              r_j=rr0+j*dr                                  !注意：r默认是桨叶半径，保持不变，此处用新变量r_j来积分
  !              r_r=r_j/rr
  !
  !              if(mu>0.000001)then
  !                  vi=viu(i,j)
  !              else
  !                  vi=viu_hover(j)
  !              end if
  !              dtheta=theta_f(1)-theta_f(2)*cosd(psi+gamma)-theta_f(3)*sind(psi+gamma)      !周期变距***************************************************
  !              phi=phi70+dphi0*(r_r-0.7)+dtheta                   !!注意全部转化成弧度计算
  !
  !              wx=omegaU*(r_j)+vx*sind(psi)
  !              wy=vz+vi        
  !              if((abs(wy)<0.0001).and.(abs(wx)<0.0001))then
  !                  betax=0.0
  !              else
  !                  betax=atand(wy/wx)
  !              end if                                    !
  !              alphax=phi-betax
  !              w_xy=sqrt(wx**2+wy**2)               
  !              if (wx>0.0) then
  !                  cy=RCl(alphax,w_xy)
  !                  cx=RCD(alphax,w_xy)
  !              else
  !                  cy=RCl(180.0-abs(alphax),w_xy)
  !                  cx=RCD(180.0-abs(alphax),w_xy)
  !              end if
  !
  !              dy=(kk*rho/4.0/pi)*cy*cc*(wx**2+wy**2)*(2.0*pi/np)*dr
  !              dx=(kk*rho/4.0/pi)*cx*cc*(wx**2+wy**2)*(2.0*pi/np)*dr
  !              dt=dy*cosd(betax)-dx*sind(abs(betax))
  !              if(betax<0.0) then
  !                  dq=dy*sind(betax)-dx*cosd(betax)
  !              else
  !                  dq=dy*sind(betax)+dx*cosd(betax)
  !              end if
  !   
  !              dyy(1)=dt
  !              dyy(2)=0.0                                  
  !              dyy(3)=0.0                         
  !              dyy(4)=dq*r_j                             
  !              dyy(5)=0.0-dt*r_j*sind(psi)                
  !              dyy(6)=0.0-dt*r_j*cosd(psi)                  
  !
  !              psi2=-psi
  !              if(mu>0.000001)then
  !                  vi2=vil(i,j)
  !              else
  !                  vi2=vil_hover(j)
  !              end if
  !              dtheta2=theta_f(4)-theta_f(5)*cosd(psi2+gamma)-theta_f(6)*sind(psi2+gamma)   !周期变距********
  !              phi2=phi70+dphi0*(r_r-0.7)+dtheta2                   !注意全部转化成弧度计算
  !              !下面为求解来流速度
  !              wx2=omegaL*(r_j)+vx*sind(psi2)
  !              wy2=vz+vi2           
  !              if((abs(wy)<0.0001).and.(abs(wx)<0.0001))then
  !                  betax2=0.0
  !              else
  !                  betax2=atand(wy2/wx2)
  !              end if                                    !
  !              alphax2=phi2-betax2
  !              w_xy2=sqrt(wx2**2+wy2**2)               !注意乘以omega_rr的平方
  !              if (wx2>0.0) then
  !                  cy2=RCl(alphax2,w_xy2)
  !                  cx2=RCD(alphax2,w_xy2)
  !              else
  !                  cy2=RCl(180.0-abs(alphax2),w_xy2)
  !                  cx2=RCD(180.0-abs(alphax2),w_xy2)
  !              end if
  !              !factor_r(2)=(2.0/pi)*acos(exp(0.5*kk*(r_r-1.0)/(abs(wy2)/omega_rr)))!#######################
  !              !cy2=factor_r(2)*cy2
  !              if((r_r>0.97).and.(mu>0.000001)) then
  !                  cy2=0.0
  !              end if
  !              dy2=(kk*rho/4.0/pi)*cy2*cc*(wx2**2+wy2**2)*(2.0*pi/np)*dr
  !              dx2=(kk*rho/4.0/pi)*cx2*cc*(wx2**2+wy2**2)*(2.0*pi/np)*dr
  !              dt2=dy2*cosd(betax2)-dx2*sind(abs(betax2))
  !              !dq2=dx2*cosd(betax2)+dy2*sind(betax2)
  !              if(betax<0.0) then
  !                  dq2=dy2*sind(betax2)-dx2*cosd(betax2)
  !              else
  !                  dq2=dy2*sind(betax2)+dx2*cosd(betax2)
  !              end if
  !              !write(42,"(4X,F14.9,\)")r_j/rr,vi,vi2,betax,betax2,alphax,alphax2,cy,cy2,cx,cx2,dt,dt2,dq,dq2
  !              !write(42,*)'    '
  !              !
  !              dyy2(1)=dt2
  !              dyy2(2)=0.0          !
  !              dyy2(3)=0.0               !
  !              dyy2(4)=dq2*r_j          
  !              dyy2(5)=0.0-dt2*r_j*sind(-psi2)   
  !              dyy2(6)=0.0-dt2*r_j*cosd(-psi2)  
  !
  !              do k=1,6,1
  !                  yy(k)=yy(k)+dyy(k)!注意不可用dpsi,因为它是角度制的,不用再乘以dr,dyy中已经含有dr
  !                  yy2(k)=yy2(k)+dyy2(k)
  !              end do
  !          end do
  !      end do
  !      !close(20)
  !      !无量纲化
  !      do k=1,6,1
  !          if(k==1) then
  !              yy(k)=yy(k)/(rho*pi*rr**2*omega_rr**2)
  !              yy2(k)=yy2(k)/(rho*pi*rr**2*omega_rr**2)
  !          else
  !              yy(k)=yy(k)/(rho*pi*rr**3*omega_rr**2)
  !              yy2(k)=yy2(k)/(rho*pi*rr**3*omega_rr**2)        
  !          end if
  !      end do
  !
		!y_f(1)=yy(1)+yy2(1)
		!y_f(2)=yy(4)-yy2(4)
		!y_f(3)=abs(yy(5)-yy2(5))/abs(y_f(1))
		!y_f(4)=yy(5)+yy2(5)
		!y_f(5)=yy(6)+yy2(6)
		!y_f(6)=yy(6)-yy2(6)
		!!
  !      return
  !  end function get_y_mid
   ! function get_y_final(theta_orig)result(theta_final)
   !     implicit none
   !     real(8)::theta_orig(6),y_f_goal(6),theta_final(6)
   !     real(8)::y_final(6)
   !     real(8)::x0(6),x1(6),dx(6),y0(6),y1(6)
   !     real(8)::jac(6,6)
   !     real(8)::residual_y(6),residual
   !     real(8)::ijac(6,6),xmkl(6),ymkl(6)
   !     integer::i,j,k,m,i_final
   !     dx=(/1.0,0.5,0.5,1.0,0.5,0.5/)
   !     i_final=0
   !     y_f_goal=y_goal
   !     x0=theta_orig    !
   !     i_final=0
   !     residual=1.0
   !     do while(abs(residual)>epsilon_final)
   !         i_final=i_final+1
   !
   !         y0=get_y_mid(x0)
   !         do j=1,6,1
   !             do m=1,6,1
   !                 if(m==j)then
   !                     x1(m)=x0(m)+dx(m)
   !                 else
   !                     x1(m)=x0(m)
   !                 end if
   !             end do
   !             y1=get_y_mid(x1)
   !             do k=1,6,1
   !                 jac(k,j)=(y1(k)-y0(k))/(x1(j)-x0(j))
   !             end do
   !         end do
   !         x1=x0+dx
   !         y1=get_y_mid(x1)
   !         residual_y=y_f_goal-y1
   !         !
   !         residual_y=abs(residual_y)
   !         residual=max(residual_y(1),residual_y(2),residual_y(3),residual_y(4),residual_y(5),residual_y(6))
   !         !
   !         write(*,*)'res_jacobi',i_final,residual
   !         write(*,*)'x1'
   !         write(*,*)x1
   !         write(43,*)'res_jacobi',i_final,residual
   !         write(43,*)'x1'
   !         write(43,*)x1
   !         call inverse(jac,ijac)
   !         ymkl=y_f_goal-y1
   !         call gemv(ijac,ymkl,xmkl)
   !         dx=xmkl*relax_factor
   !
			!do i=1,6,1
			!	if(abs(dx(i))>max_dx)then
			!		dx(i)=sign(max_dx,dx(i))
			!	end if
			!end do
   !         !
   !         x0=x1
   !     end do
   !     y_final=y1
   !     theta_final=x1
   !     write(*,*)'y_final',y_final
   !     write(43,*)'y_final',y_final
   !     return
   ! end function get_y_final
end module get_vi_y