function [R_eigen,t_eigen] = eigen_solver(z_h,y_h,R_prio,K)

f_pix=K(1,1);  % focal length
ix_pix=2*K(1,3); % length of camera plane
u_pix=K(1,3);
iy_pix=2*K(2,3); % height of camera plane
v_pix=K(2,3);

z_h_n=normc(z_h);
y_h_n=normc(y_h);
point_num=size(z_h,2);
T_eigen=opengv2('eigensolver', 1:point_num, z_h_n  , y_h_n, R_prio);
R_eigen=T_eigen(:,1:3);
t_eigen=T_eigen(:,4)/norm(T_eigen(:,4));
E_eigen=cal_E(R_eigen,t_eigen);

intrinsics = cameraIntrinsics([f_pix f_pix],[u_pix v_pix],[ix_pix iy_pix]);
[R_eigen,t_eigen] = relativeCameraPose(E_eigen,intrinsics,f_pix*y_h(1:2,:)'+[u_pix v_pix],f_pix*z_h(1:2,:)'+[u_pix v_pix]);
R_eigen=R_eigen(:,:,1);
t_eigen=-R_eigen*t_eigen(1,:)';