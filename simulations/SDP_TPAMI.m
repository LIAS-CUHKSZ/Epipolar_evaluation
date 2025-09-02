function [E_est,R_est,t_est]=SDP_TPAMI(y_h,z_h,K)


f_pix=K(1,1);  % focal length
ix_pix=2*K(1,3); % length of camera plane
u_pix=K(1,3);
iy_pix=2*K(2,3); % height of camera plane
v_pix=K(2,3);

m=size(y_h,2);
C=zeros(9,9);
for i=1:m
    C=C+kron(y_h(:,i),z_h(:,i))*kron(y_h(:,i),z_h(:,i))';
end
C0=[C zeros(9,3); zeros(3,9) zeros(3,3)];
A1=diag([1 1 1 0 0 0 0 0 0 0 -1 -1]);
A2=diag([0 0 0 1 1 1 0 0 0 -1 0 -1]);
A3=diag([0 0 0 0 0 0 1 1 1 -1 -1 0]);
A4=zeros(12,12);
A4(1:3,4:6)=eye(3)/2;
A4(4:6,1:3)=eye(3)/2;
A4(10,11)=1/2;
A4(11,10)=1/2;
A5=zeros(12,12);
A5(1:3,7:9)=eye(3)/2;
A5(7:9,1:3)=eye(3)/2;
A5(10,12)=1/2;
A5(12,10)=1/2;
A6=zeros(12,12);
A6(7:9,4:6)=eye(3)/2;
A6(4:6,7:9)=eye(3)/2;
A6(12,11)=1/2;
A6(11,12)=1/2;
A7=diag([0 0 0 0 0 0 0 0 0 1 1 1]);

cvx_begin sdp quiet
variable X(12,12) symmetric 
minimize (trace(C0*X))
subject to
trace(A1*X)==0;
trace(A2*X)==0;
trace(A3*X)==0;
trace(A4*X)==0;
trace(A5*X)==0;
trace(A6*X)==0;
trace(A7*X)==1;
X>=0;
cvx_end

Xe=X(1:9,1:9);
[V,~]=eig(Xe);   
e_est=V(:,9);
E_est=[e_est(1:3)';e_est(4:6)';e_est(7:9)'];

%  [U,~,V] = svd(E_est);
% %  E_est=U*diag([1 1 0])*V'/sqrt(2);
% W=[0 -1 0;1 0 0;0 0 1];

intrinsics = cameraIntrinsics([f_pix f_pix],[u_pix v_pix],[ix_pix iy_pix]);
[R_est,t_est] = relativeCameraPose(E_est,intrinsics,f_pix*z_h(1:2,:)'+[u_pix v_pix],f_pix*y_h(1:2,:)'+[u_pix v_pix]);

R_est=R_est(:,:,1)';
t_est=t_est(1,:)';
% R_est=R_est(:,:,1);
% t_est=-R_est*t_est(1,:)';