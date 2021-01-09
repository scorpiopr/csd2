module global                                               !ȫ�ֱ���
    integer::kk=3     !��ҶƬ����kk��
    real(8)::rr,r_cc,cc,phi70,dphi0,omega_rr,omega,h_r,hh,r_cut,rr0,aa    !����Ľ�Ҷ����
    real(8)::mu,alpha_s,gamma,lambd0,lambd02,vx,vz                !��������������ǰ���ȣ�������ӭ��,��ǰ���ݽǣ���ʱû�õ���,�����
    real(8),parameter::rho=1.139,g=9.80665   !�㶨������pi,�����ܶ�(kg/m3),�������ٶ�
    real(8)::ct,cmn,ee,cml_m,cmm_p,cmm_m
    real(8)::ct1,ct2,cml01,cmm01,cml02,cmm02
    real(8)::y_goal(6)           !Ŀ���������ֲ���
	integer::flag_goal
    real(8)::th0u,a1u,b1u,th0l,a1l,b1l
    integer,parameter::np=180,nr=50                   !���ֶַ�����R�ķֶΣ���λ�Ǧ׵ķֶ�
    real(8)::dr,dpsi
    real(8)::epsilon_final,relax_factor                     !������������ 
	real(8)::max_dx=10.0
    contains
        !����input10.dat������
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
        !�ǶȺͻ���ת������
        function angle_rad(angle)result(rad)!�Ƕ�->����
            implicit none
            real(8)::rad,angle
            rad=(angle*3.141592653)/180.0
            return
        end function angle_rad
        function rad_angle(rad)result(angle)!����->�Ƕ�
            implicit none
            real(8)::angle,rad
            angle=rad*180.0/3.141592653
            return
        end function rad_angle
end module global