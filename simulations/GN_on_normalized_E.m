function [E_GN,R_GN,t_GN] = GN_on_normalized_E(z_h,y_h,E,K)
f_pix=K(1,1);  % focal length
ix_pix=2*K(1,3); % length of camera plane
u_pix=K(1,3);
iy_pix=2*K(2,3); % height of camera plane
v_pix=K(2,3);
m=size(z_h,2);
[U,E0,V]=svd(E);
if det(U)==-1
    U(:,3)=-U(:,3);
end
if det(V)==-1
    V(:,3)=-V(:,3);
end

bar_M=zeros(m,9);
for i=1:m
    bar_M(i,:)=vec(z_h(:,i)*y_h(:,i)')';
end
mathcal_M=(bar_M'*bar_M)/m;

Qx=[0 0 0;0 0 -1;0 1 0];
Qy=[0 0 1;0 0 0;-1 0 0];
Qz=[0 -1 0;1 0 0;0 0 0];
Q1=[vec(Qx/sqrt(2)) vec(Qy/sqrt(2)) vec(Qz/2) zeros(9,2)];
Q2=[zeros(9,2) vec(-Qz/2) vec(Qx/sqrt(2)) vec(Qy/sqrt(2))];
J=kron(V,U)*(kron(E0,eye(3))*Q1-kron(eye(3),E0)*Q2);
first_deriv=J'*mathcal_M*vec(E);
second_deriv=J'*mathcal_M*J;
if cond(second_deriv)>10^15
    x=-(second_deriv+10^(-8)*eye(5))\first_deriv;
    disp('near singular')
else
    x=-second_deriv\first_deriv;
end

omega1=[0 -x(3)/sqrt(2) x(2);x(3)/sqrt(2) 0 -x(1);-x(2) x(1) 0]/sqrt(2);
omega2=[0 x(3)/sqrt(2) x(5);-x(3)/sqrt(2) 0 -x(4);-x(5) x(4) 0]/sqrt(2);
U_GN=U*expm(omega1);
V_GN=V*expm(omega2);
E_GN=U_GN*E0*V_GN';

intrinsics = cameraIntrinsics([f_pix f_pix],[u_pix v_pix],[ix_pix iy_pix]);
[R_GN,t_GN] = relativeCameraPose(E_GN,intrinsics,f_pix*y_h(1:2,:)'+[u_pix v_pix],f_pix*z_h(1:2,:)'+[u_pix v_pix]);
R_GN=R_GN(:,:,1);
t_GN=-R_GN*t_GN(1,:)';
