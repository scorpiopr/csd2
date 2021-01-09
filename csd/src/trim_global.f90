module global                                               !全局变量
    integer::kk=3     !桨叶片数（kk）
    real(8)::rr,r_cc,cc,phi70,dphi0,omega_rr,omega,h_r,hh,r_cut,rr0,aa    !读入的桨叶参数
    real(8)::mu,alpha_s,gamma,lambd0,lambd02,vx,vz                !读入的旋翼参数：前进比，旋翼构造迎角,提前操纵角（暂时没用到）,流入比
    real(8),parameter::rho=1.139,g=9.80665   !恒定变量：pi,空气密度(kg/m3),重力加速度
    real(8)::ct,cmn,ee,cml_m,cmm_p,cmm_m
    real(8)::ct1,ct2,cml01,cmm01,cml02,cmm02
    real(8)::y_goal(6)           !目标量，保持不变
	integer::flag_goal
    real(8)::th0u,a1u,b1u,th0l,a1l,b1l
    integer,parameter::np=180,nr=50                   !积分分段数：R的分段，方位角ψ的分段
    real(8)::dr,dpsi
    real(8)::epsilon_final,relax_factor                     !迭代收敛控制 
	real(8)::max_dx=10.0
    contains
        !读入input10.dat的数据
        subroutine read_input()
            implicit none
            open(10,file='10input.dat')
            read(10,*)
            read(10,*)
            read(10,*)
            read(10,*)r_cc,cc,phi70,dphi0,kk,omega_rr,h_r,r_cut
            read(10,*)
            read(10,*)
            read(10,*)mu,alpha_s,gamma
            read(10,*)
            read(10,*)
            read(10,*)ct,cmn,ee,cml_m,cmm_p,cmm_m
            read(10,*)
            read(10,*)
            read(10,*)th0u,a1u,b1u,th0l,a1l,b1l
            read(10,*)
            read(10,*)
            read(10,*)epsilon_final,relax_factor 
            read(10,*)
            write(10,*)'end!**************************************************************************'
            close(10)
            return
        end subroutine read_input    
        !角度和弧度转化函数
        function angle_rad(angle)result(rad)!角度->弧度
            implicit none
            real(8)::rad,angle
            rad=(angle*3.141592653)/180.0
            return
        end function angle_rad
        function rad_angle(rad)result(angle)!弧度->角度
            implicit none
            real(8)::angle,rad
            angle=rad*180.0/3.141592653
            return
        end function rad_angle
end module global