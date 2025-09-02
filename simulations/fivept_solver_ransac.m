function [R,t] = fivept_solver (z_h,y_h,K)

f_pix=K(1,1);  % focal length
ix_pix=2*K(1,3); % length of camera plane
u_pix=K(1,3);
iy_pix=2*K(2,3); % height of camera plane
v_pix=K(2,3);

z_h_n=normc(z_h);
y_h_n=normc(y_h);
all_solution=opengv('fivept_nister_ransac', z_h_n  , y_h_n);
R=all_solution(:,1:3);
t=all_solution(:,4);
t=t/norm(t);



% size_all_solution=size(all_solution);
% % E=all_solution(:,:,size_all_solution(3));
% E=all_solution(:,:,1);
% 
% intrinsics = cameraIntrinsics([f_pix f_pix],[u_pix v_pix],[ix_pix iy_pix]);
% [R,t] = relativeCameraPose(E,intrinsics,f_pix*y_h(1:2,:)'+[u_pix v_pix],f_pix*z_h(1:2,:)'+[u_pix v_pix]);
% R=R(:,:,1);
% t=-R*t(1,:)';
