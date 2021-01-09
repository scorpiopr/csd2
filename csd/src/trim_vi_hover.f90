module vi_hover
    use LinearAlgebra
    use global
    use cl_cd
    !use mkl95_precision
    !use mkl95_blas
    implicit none
    
    REAL(8) :: pi=3.141592653
    private pi
    
    contains
!***************垂直飞行状态下*****************

function get_ep(z)
    implicit none
    real(8)::z,get_ep
    get_ep=1.0+z/sqrt(rr**2+z**2)
    return
end function get_ep

function get_flag(rs)
    implicit none
    integer::get_flag
    real(8)::rs
    if((rs>=rr0).and.(rs<rr))then
        get_flag=int((rs-rr0)/((rr-rr0)/nr))+1

    else if(rs==rr)then
        get_flag=nr
    else
        get_flag=0
    end if
    return
end function get_flag

subroutine get_wd(wu,wl,wud,wld)
    implicit none
    integer::i,f
    real(8)::r1,r2,r3
    real(8)::wu(nr),wl(nr),wld(nr),wud(nr)
    do i=1,nr
        r1=rr0+(rr-rr0)/nr*i
        r2=r1*sqrt((vz+wu(i)+wld(i))/(vz+get_ep(hh)*(wu(i)+wld(i))))
        f=get_flag(r2-(rr-rr0)/nr/2.0)
        if (f/=0) then
            wld(i)=get_ep(-hh)*wl(f)
        else
            wld(i)=0.0
        end if
        r3=r1*sqrt((vz+wl(i)+wud(i))/(vz+get_ep(-hh)*(wl(i)+wud(i))))
        f=get_flag(r3-(rr-rr0)/nr/2.0)
        if(f/=0) then
            wud(i)=get_ep(hh)*wu(f)
        else
            wud(i)=0.0
        end if
    end do
    return
end subroutine get_wd

function fun(i,theta_f,w,wd)
    implicit none
    integer::i
    real(8)::theta_f(3),w,wd,fun
    !中间变量
    real(8)::rs,phi,alphax,betax,cy
    real(8)::wx,wy,w_xy,y1,y2,r_r,factor_r
    factor_r=1.0
    !
    rs=(rr0+(rr-rr0)*i/nr)
    r_r=rs/rr
    !write(*,*)'w,wd,rr,Vz,omega,rho,theta7,dtheta,Cl,c,R'
    !write(*,*)w,wd,rs,Vz,omega,rho,theta_f(1),dphi0,cc,rr,rr0
    phi=phi70+dphi0*(r_r-0.7)+theta_f(1)
    wx=omega*rs
    wy=vz+w+wd        
    w_xy=sqrt(wx**2+wy**2)
    betax=atand(wy/wx)
    alphax=phi-betax
    cy=RCl(alphax,w_xy)

    factor_r=(2.0/pi)*acos(exp(0.5*kk*(r_r-1.0)/(abs(wy)/omega_rr)))
    !write(41,*)'r,factor_r,pi,kk,wy,omega_rr'
    !write(41,*),r_r,factor_r,pi,kk,wy,omega_rr
    !write(*,*)'cy',cy
    y1=4.0*pi*rho*rs*abs(wy)*w
    y2=factor_r*0.5*kk*rho*cc*(wx**2+wy**2)*cy/cosd(betax)
    fun=y1-y2
    return
end function fun
!运用动量—叶素方程计算自诱导速度和互诱导速度
subroutine get_vi_hover(theta,viu,vil)
    implicit none
    real(8)::theta(6),viu(nr),vil(nr),thetaU(3),thetaL(3)
    real(8)::wu_orig(nr),wl_orig(nr),wu(nr),wl(nr),wud(nr),wld(nr),ad(2),ep_ul(2)

    real(8)::a1,a2,a3,y1,y2,y3,theta_f(3)
    real(8)::res_vi,residual_vi
    real(8),parameter::vi_ep=0.001,fun_ep=0.05,relax_vi=0.9
    integer::i,j,i_vi
    !write(*,*)'Start get_vi_hover !'
    do i=1,nr,1
        wu(i)=1.0;wl(i)=1.0;wud(i)=0.0;wld(i)=0.0
        wu_orig(i)=0.0;wl_orig(i)=0.0;
    end do
    call get_wd(wu,wl,wud,wld)
    thetaU=(/theta(1),theta(2),theta(3)/)
    thetaL=(/theta(4),theta(5),theta(6)/)
    !
    i_vi=0
    residual_vi=10.0
    do while(residual_vi>vi_ep)
        i_vi=i_vi+1
        !write(*,*)'do res_vi',i_vi
        wu_orig=wu;wl_orig=wl
        !
        do i=1,nr,1
            !write(*,*)'nr',i
            !上旋翼
            theta_f=thetaU
            a1=-30.0
            y1=1.0;y2=1.0
            do while(y1*y2>0.0)
                a1=a1+3.0
                y1=fun(i,theta_f,a1,wld(i))
                y2=fun(i,theta_f,a1+3.0,wld(i))
            end do
            a2=a1+3.0
            !write(*,*)'y1,y2',fun(i,theta_f,a1,wld(i)),fun(i,theta_f,a2,wld(i))
            y1=10.0;y2=0.0
            do while((abs(y2-y1)>fun_ep).and.(abs(a2-a1)>vi_ep))
                a3=0.5*(a1+a2)
                y1=fun(i,theta_f,a1,wld(i))
                y2=fun(i,theta_f,a2,wld(i))
                y3=fun(i,theta_f,a3,wld(i))
                if(abs(y1)<fun_ep)then
                    a3=a1;exit
                else if(abs(y2)<fun_ep)then
                    a3=a2;exit
                else if(abs(y3)<fun_ep)then
                    exit!a3不变
                end if
                !以上表示如果找到解就跳出循环
                if(y1*y3<0.0)then
                    a2=a3
                else if(y2*y3<0.0)then
                    a1=a3
                end if
                !write(*,*)'a2-a1',a2-a1
            end do!(do while(y2-y1))
            wu(i)=relax_vi*a3+(1.0-relax_vi)*wu(i)
            call get_wd(wu,wl,wud,wld)
            !下旋翼
            theta_f=thetaL
            a1=-30.0
            y1=1.0;y2=1.0
            do while(y1*y2>0.0)
                a1=a1+3.0
                y1=fun(i,theta_f,a1,wud(i))
                y2=fun(i,theta_f,a1+3.0,wud(i))
            end do
            a2=a1+3
            !
            y1=10.0;y2=0.0
            a3=0.5*(a1+a2)
            do while((abs(y2-y1)>fun_ep).and.(abs(a2-a1)>vi_ep))
                a3=0.5*(a1+a2)
                y1=fun(i,theta_f,a1,wud(i))
                y2=fun(i,theta_f,a2,wud(i))
                y3=fun(i,theta_f,a3,wud(i))
                if(abs(y1)<fun_ep)then
                    a3=a1;exit
                else if(abs(y2)<fun_ep)then
                    a3=a2;exit
                else if(abs(y3)<fun_ep)then
                    exit!a3不变
                end if
                !以上表示如果找到解就跳出循环
                if(y1*y3<0.0)then
                    a2=a3
                else if(y2*y3<0.0)then
                    a1=a3
                end if
            end do!(do while(y2-y1))
            wl(i)=relax_vi*a3+(1.0-relax_vi)*wl(i)
            !call get_wd(wu,wl,wud,wld)
        end do!
        !
        residual_vi=max(abs(wu(1)-wu_orig(1)),abs(wl(1)-wl_orig(1)))
        do i=1,nr,1
            res_vi=max(abs(wu(i)-wu_orig(i)),abs(wl(i)-wl_orig(i)))
            if(res_vi>residual_vi)then
                residual_vi=res_vi
            end if
        end do
        !write(*,*)'res_vi_hover',i_vi,residual_vi
    end do
    write(*,*)'res_vi_hover',i_vi,residual_vi
    write(43,*)'res_vi_hover',i_vi,residual_vi
    !
    call get_wd(wu,wl,wud,wld)
    ad(1)=get_ad_hover(wld)
    ad(2)=get_ad_hover(wud)
    ep_ul(1)=get_ep(hh)*ad(2)/aa   
    ep_ul(2)=get_ep(-hh)*ad(1)/aa
    viu=wu+wld
    vil=wl+wud
    
    write(30,*)'r   wu  wl  wud wld viu vil'
    do j=1,nr
        write(30,"(4X,F14.9,\)")(rr0+j*dr)/rr,wu(j),wl(j),wud(j),wld(j),wu(j)+wld(j),wl(j)+wud(j)
        write(30,*)'    '
    end do
    return
end subroutine get_vi_hover

    function get_ad_hover(wd)
        implicit none
        real(8)::get_ad_hover,wd(nr)
        integer::j
        get_ad_hover=0.0
        do j=1,nr,1
            if(wd(j)/=0.0) then
                get_ad_hover=get_ad_hover+pi*((rr0+(rr-rr0)*j/nr)**2-(rr0+(rr-rr0)*(j-1)/nr)**2)
            end if
        end do
        return
    end function
end module vi_hover