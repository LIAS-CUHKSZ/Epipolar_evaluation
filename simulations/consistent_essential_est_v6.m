function [E_est,R_est,t_est,R_GN,t_GN,var_est] = consistent_essential_est_v6(z_h,y_h,K) % K is the intrinsic matrix
m=size(z_h,2); % number of point
f_pix=K(1,1);  % focal length
ix_pix=2*K(1,3); % length of camera plane
u_pix=K(1,3);
iy_pix=2*K(2,3); % height of camera plane
v_pix=K(2,3);
z=z_h(1:2,:);
y=y_h(1:2,:);
W=[1 0 0;0 1 0];
e3=[0 0 1]';
e1=[1 0 0]';
e2=[0 1 0]';
%% noise variance estimation
bar_A=zeros(m,9);
bar_Y=zeros(3,3);
for i=1:m
    bar_A(i,:)=kron(y_h(:,i)',z_h(:,i)');
    bar_Y=bar_Y+y_h(:,i)*y_h(:,i)'/m;
end
Q_m=bar_A'*bar_A/m;
S_m=kron(bar_Y,[W;0 0 0]);
var_est=1/max(eig(Q_m\S_m));

%% estimation via SVD
Q_BE=Q_m-var_est*S_m;

% format shortE
% eig(Q_BE)

[V,~]=eig(Q_BE);   
theta_est=V(:,4);    % will different matlab versions have different orders of eigevalues?
E_est=reshape(theta_est,3,3);
intrinsics = cameraIntrinsics([f_pix f_pix],[u_pix v_pix],[ix_pix iy_pix]);
[R_est,t_est] = relativeCameraPose(E_est,intrinsics,f_pix*y'+[u_pix v_pix],f_pix*z'+[u_pix v_pix]);
R_est=R_est(:,:,1);
t_est=-R_est*t_est(1,:)';  % we adopt z=Ry+t, while this toolbox adopts z=Ry-Rt !!!!!!

%%
%Gauss-newton iteration
gamma0=asin(t_est(3));
if t_est(1)>0
    theta0=atan(t_est(2)/t_est(1));
else
    theta0=atan(t_est(2)/t_est(1))+pi;
end
partial_t_theta=[-cos(gamma0)*sin(theta0) cos(gamma0)*cos(theta0) 0]';
partial_t_gamma=[-sin(gamma0)*cos(theta0) -sin(gamma0)*sin(theta0) cos(gamma0)]';

J=zeros(2*m,5);
z_predic=zeros(2*m,1);
phi = [0 0 0;0 0 1;0 -1 0;0 0 -1;0 0 0;1 0 0;0 1 0;-1 0 0;0 0 0];
triple_t=kron(eye(3),t_est); % 9 by 3
pxii=kron(y_h',R_est)*phi; 
Ryih_I=kron(R_est*y_h,eye(3));
Ryih_t=kron(R_est*y_h,t_est);
for i = 1:m  
    Ryih=R_est*y_h(:,i); % 3 by 1
    matrix1=[[-e3';0 0 0;e1'] [0 0 0;-e3';e2'] [zeros(2,2) z(:,i);-z(:,i)' 0]]; % 3 by 9
    matrix2=[[e3';0 0 0;-z(1,i)*e3'] [0 0 0;e3';-z(2,i)*e3'] [-eye(2) z(:,i);0 0 0]]; % 3 by 9
    denominator=t_est'*matrix2*triple_t*Ryih;
    numerator=Ryih'*matrix1*triple_t*Ryih;
    partial_ki_s=(denominator*(Ryih'*(matrix1*triple_t+triple_t'*matrix1')*pxii(3*i-2:3*i,:))-numerator*(t_est'*matrix2*triple_t*pxii(3*i-2:3*i,:)))/(denominator^2);
    
    
    partial_ki_theta=(denominator*(Ryih'*matrix1*Ryih_I(:,3*i-2:3*i)*partial_t_theta)-numerator*(t_est'*matrix2*Ryih_I(:,3*i-2:3*i)+Ryih_t(:,i)'*matrix2')*partial_t_theta)/(denominator^2);
    partial_ki_gamma=(denominator*(Ryih'*matrix1*Ryih_I(:,3*i-2:3*i)*partial_t_gamma)-numerator*(t_est'*matrix2*Ryih_I(:,3*i-2:3*i)+Ryih_t(:,i)'*matrix2')*partial_t_gamma)/(denominator^2);

    scale=numerator/denominator;
    g=W* (R_est*y_h(:,i)+scale*t_est);
    h=e3'* (R_est*y_h(:,i)+scale*t_est);
    z_predic(2*i-1:2*i)=g/h;
    J(2*i-1:2*i,1:3) = (h*W*(pxii(3*i-2:3*i,:)+t_est*partial_ki_s)-g*e3'*(pxii(3*i-2:3*i,:)+t_est*partial_ki_s))/(h^2);
    J(2*i-1:2*i,4) =(h*W*(scale*partial_t_theta+partial_ki_theta*t_est)-g*e3'*(scale*partial_t_theta+partial_ki_theta*t_est))/(h^2);
    J(2*i-1:2*i,5) =(h*W*(scale*partial_t_gamma+partial_ki_gamma*t_est)-g*e3'*(scale*partial_t_gamma+partial_ki_gamma*t_est))/(h^2);
end
initial = [0;0;0;theta0;gamma0];
results = initial + (J' * J)\(J' * (z(:) - z_predic));
X_GN = results(1:3);
t_GN = [cos(results(5))*cos(results(4)) cos(results(5))*sin(results(4)) sin(results(5))]';
Xhat = [0 -X_GN(3) X_GN(2); X_GN(3) 0 -X_GN(1); -X_GN(2) X_GN(1) 0];
R_GN = R_est * expm(Xhat);  