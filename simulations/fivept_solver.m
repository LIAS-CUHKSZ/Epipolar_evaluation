function [R,t] = fivept_solver (z_h,y_h,K)

f_pix=K(1,1);  % focal length
ix_pix=2*K(1,3); % length of camera plane
u_pix=K(1,3);
iy_pix=2*K(2,3); % height of camera plane
v_pix=K(2,3);

z_h_n=normc(z_h);
y_h_n=normc(y_h);
all_solution=opengv('fivept_nister', 1:5, z_h_n  , y_h_n);
size_all_solution=size(all_solution);
E_temp=cal_E(eye(3),[1 0 0]');
if size_all_solution(3)==0
solution_num=size_all_solution(3)
E_temp=[-0.270837610209851	-0.245617347464686	0.446824843594849;
    0.500121451210174	-0.217541002594527	-0.189440701970119;
    -0.229283841000323	0.463158350059213	-0.257384141624731];
end
error_temp=trace((z_h'*E_temp*y_h).^2);
for i=1:size_all_solution(3)
    if trace((z_h'*all_solution(:,:,i)*y_h).^2)<error_temp
        E_temp=all_solution(:,:,i);
        error_temp=trace((z_h'*all_solution(:,:,i)*y_h).^2);
    end
end
E=E_temp;
intrinsics = cameraIntrinsics([f_pix f_pix],[u_pix v_pix],[ix_pix iy_pix]);
[R,t] = relativeCameraPose(E,intrinsics,f_pix*y_h(1:2,:)'+[u_pix v_pix],f_pix*z_h(1:2,:)'+[u_pix v_pix]);
R=R(:,:,1);
t=-R*t(1,:)';
