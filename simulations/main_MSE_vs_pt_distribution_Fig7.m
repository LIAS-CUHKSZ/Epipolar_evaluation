% ki is a function of R and t_normal, and t_normal is constrained to a 2-sphere
% noise is only added to the second camera image
% compasion with the CRB, the CRB is derived under constraints of SO(3), 2-sphere, and ki being a function of R and t_normal
% the two-step estimator can asymptotically achieve the CRB
clc
clear
%%
%true pose and parameters
alpha=20;
beta=20;
gamma=20;
true_R=cal_rotation(alpha,beta,gamma);
true_t=[1 1 1]';
true_t_normal=true_t/norm(true_t);
true_E=[cross(true_t,true_R(:,1)) cross(true_t,true_R(:,2)) cross(true_t,true_R(:,3))];
true_E_normal=true_E/norm(true_E,'fro');
true_theta=vec(true_E_normal);
f=0.05*20; % focal length
f_pix=800;
ix=0.08*20; % length of camera plane
ix_pix=640;
u_pix=320;
iy=0.06*20; % height of camera plane
iy_pix=480;
v_pix=240;
max_dis=5*20; % max depth of point that cameras can see
min_dis=[5 4 3.5 3 2.5 2 1.5 1]*20; % min depth of point that cameras can see
K=[f_pix 0 u_pix;0 f_pix v_pix;0 0 1];
sigma_pix=0.5;
sigma=sigma_pix/f_pix;
m=1000;  % the number of points

%%
% initilization
MSE_R_LM=zeros(1,length(min_dis));
MSE_t_LM=zeros(1,length(min_dis));
MSE_R_eigen=zeros(1,length(min_dis));
MSE_t_eigen=zeros(1,length(min_dis));
MSE_R_SDP=zeros(1,length(min_dis));
MSE_t_SDP=zeros(1,length(min_dis));
MSE_R_GNE=zeros(1,length(min_dis));
MSE_t_GNE=zeros(1,length(min_dis));
MSE_R=zeros(1,length(min_dis));
MSE_t=zeros(1,length(min_dis));
MSE_t_GN=zeros(1,length(min_dis));
MSE_R_GN=zeros(1,length(min_dis));
CRB_R=zeros(1,length(min_dis));
CRB_t=zeros(1,length(min_dis));

min_eigen=zeros(1,length(min_dis));

monte_carlo=1;  % monte carlo number for calculating the MSE

%%
% test begins
for j=1:length(min_dis)
x=zeros(3,m); % 3D points
x_h=zeros(4,m); % 3D points
y=zeros(2,m);  % 2D points in the first camera plane
z=zeros(2,m);  % 2D points in the second camera plane
z_noisefree=zeros(2,m);

for k=1:monte_carlo
%% 
% generate 2D and 3D points
counter=0;
while counter<m
    y_temp=[rand*ix-ix/2 rand*iy-iy/2]';
    depth=min_dis(j)+rand*(max_dis-min_dis(j));
    x_temp=[y_temp*depth/f ; depth];
    x_temp2=true_R*x_temp+true_t;
    if (x_temp2(3)>f) && (abs(x_temp2(1))<ix/2*x_temp2(3)/f) && (abs(x_temp2(2))<iy/2*x_temp2(3)/f)
        counter=counter+1;
        x(:,counter)=x_temp;
        y(:,counter)=y_temp;
        z_noisefree(:,counter)=x_temp2(1:2)/x_temp2(3)*f;
        z(:,counter)=x_temp2(1:2)/x_temp2(3)*f+randn(2,1)*sigma;
    end
end
x_h=[x;ones(1,m)];
y_h=[y;ones(1,m)];
z_h=[z;ones(1,m)];
z_noisefree_h=[z_noisefree;ones(1,m)];

%%
% pose estimation
[R_5pt,t_5pt]=fivept_solver (z_h,y_h,K);
[R_eigen,t_eigen] = eigen_solver(z_h,y_h,R_5pt,K);
[~,R_SDP,t_SDP]=SDP_TPAMI(y_h,z_h,K);
E_5pt=cal_E(R_5pt,t_5pt);
[E_GNE,R_GNE,t_GNE] = GN_on_normalized_E(z_h,y_h,E_5pt,K);
[E_GNE,R_GNE,t_GNE] = GN_on_normalized_E(z_h,y_h,E_GNE,K);
[E_GNE,R_GNE,t_GNE] = GN_on_normalized_E(z_h,y_h,E_GNE,K);
[E_GNE,R_GNE,t_GNE] = GN_on_normalized_E(z_h,y_h,E_GNE,K);
[E_est,R_est,t_est,R_GN,t_GN,var_est] = consistent_essential_est(z_h,y_h,K);
[E_LM,R_LM,t_LM]=basic_LM(z_h,y_h,R_5pt,t_5pt);



%%
% calculate MSE and minimum eigenvalue
MSE_R_LM(j)=MSE_R_LM(j)+norm(R_LM-true_R,'fro')^2/monte_carlo;
MSE_t_LM(j)=MSE_t_LM(j)+norm(t_LM-true_t_normal)^2/monte_carlo;
MSE_R_eigen(j)=MSE_R_eigen(j)+norm(R_eigen-true_R,'fro')^2/monte_carlo;
MSE_t_eigen(j)=MSE_t_eigen(j)+norm(t_eigen-true_t_normal)^2/monte_carlo;
MSE_R_SDP(j)=MSE_R_SDP(j)+norm(R_SDP-true_R,'fro')^2/monte_carlo;
MSE_t_SDP(j)=MSE_t_SDP(j)+norm(t_SDP-true_t_normal)^2/monte_carlo;
MSE_R_GNE(j)=MSE_R_GNE(j)+norm(R_GNE-true_R,'fro')^2/monte_carlo;
MSE_t_GNE(j)=MSE_t_GNE(j)+norm(t_GNE-true_t_normal)^2/monte_carlo;
MSE_R(j)=MSE_R(j)+norm(R_est-true_R,'fro')^2/monte_carlo;
MSE_t(j)=MSE_t(j)+norm(t_est-true_t_normal)^2/monte_carlo;
MSE_R_GN(j)=MSE_R_GN(j)+norm(R_GN-true_R,'fro')^2/monte_carlo;
MSE_t_GN(j)=MSE_t_GN(j)+norm(t_GN-true_t_normal)^2/monte_carlo;


end
min_eigen(j)=min(eig(x_h*x_h'/m));

%%
% calculate CRB
% two remarks: 1. ki needs to be a function of R and t_normal when
% calculating the derivative of the likehood function; 2. cannot calculate
% the CRB of R and t_normal seperately by assuming the other is known. the
% fisher matrix should be 12 * 12 
W=[1 0 0;0 1 0];
e3=[0 0 1]';
e1=[1 0 0]';
e2=[0 1 0]';
F=zeros(12,12);
g_true=W* (true_R*y_h+true_t./x(3,:));
h_true=e3'* (true_R*y_h+true_t./x(3,:));
true_cov=sigma^2*eye(2);
for i=1:m
    L_i=kron(y_h(:,i)',eye(3));
    matrix1_former=[-true_t_normal(3)*eye(2);true_t_normal(1:2)'];
    matrix2_former=[zeros(2,2) z_noisefree(:,i);-z_noisefree(:,i)' 0];
    matrix3_former=[true_t_normal(3)*eye(2);-true_t_normal(3)*z_noisefree(:,i)'];
    matrix4_former=[-eye(2) z_noisefree(:,i);0 0 0];
    vector1=[y_h(:,i)'*true_R'*matrix1_former y_h(:,i)'*true_R'*matrix2_former*true_t_normal];
    vector2=[true_t_normal'*matrix3_former true_t_normal'*matrix4_former*true_t_normal];
    matrix1=[-e3';0 0 0;e1'];
    matrix2=[0 0 0;-e3';e2'];
    matrix3=[zeros(2,2) z_noisefree(:,i);-z_noisefree(:,i)' 0];
    matrix4=[e3';0 0 0;-z_noisefree(1,i)*e3'];
    matrix5=[0 0 0;e3';-z_noisefree(2,i)*e3'];
    matrix6=[-eye(2) z_noisefree(:,i);0 0 0];
%     pxii=kron(y_h(:,i)',R_est)*phi;
    Ryih=true_R*y_h(:,i);
    ki=(vector1*Ryih)/(vector2*Ryih);
    g_true(:,i)=W* (true_R*y_h(:,i)+ki*true_t_normal);
    h_true(i)=e3'* (true_R*y_h(:,i)+ki*true_t_normal);
    denominator=vector2*Ryih;
    numerator=vector1*Ryih;
    partial_ki_theta=(denominator*(vector1*L_i+Ryih'*[true_t_normal'*matrix1'*L_i;true_t_normal'*matrix2'*L_i;true_t_normal'*matrix3'*L_i])-numerator*(vector2*L_i))/(denominator^2);
    partial_ri_theta=(h_true(i)*(W*L_i+W*true_t_normal*partial_ki_theta)-g_true(:,i)*(e3'*L_i+e3'*true_t_normal*partial_ki_theta))'/(h_true(i)^2);

    partial_ki_t=(denominator*(Ryih'*[Ryih'*matrix1;Ryih'*matrix2;Ryih'*matrix3])-numerator*(Ryih'*[true_t_normal'*(matrix4+matrix4');true_t_normal'*(matrix5+matrix5');true_t_normal'*(matrix6+matrix6')]))/(denominator^2);
    partial_ri_t=(h_true(i)*(W*ki+W*true_t_normal*partial_ki_t)-g_true(:,i)*(e3'*ki+e3'*true_t_normal*partial_ki_t))'/(h_true(i)^2);
    partial_ri_thetat=[partial_ri_theta;partial_ri_t];
    F=F+partial_ri_thetat/true_cov*partial_ri_thetat';
end
O=zeros(3,1);
constraints=[2*true_R(:,1)' O' O' O';true_R(:,2)' true_R(:,1)' O' O';true_R(:,3)' O' true_R(:,1)' O';O' 2*true_R(:,2)' O' O';O' true_R(:,3)' true_R(:,2)' O';O' O' 2*true_R(:,3)' O';O' O' O' 2*true_t_normal'];
UU=null(constraints);
bar_F=UU*((UU'*F*UU)\UU');
CRB_R(j)=trace(bar_F(1:9,1:9));
CRB_t(j)=trace(bar_F(10:12,10:12));
end

%%
% save results and semilogy a pair of figures
% save result_MSE_vs_pt_distribution MSE_R_5pt MSE_t_5pt MSE_R_eigen MSE_t_eigen MSE_R_SDP MSE_t_SDP MSE_R_GNE MSE_t_GNE MSE_R MSE_t MSE_R_GN MSE_t_GN CRB_R CRB_t min_eigen

semilogy(min_eigen,MSE_R_LM,'-o','color',[0.8627 0.1961 0.1373],'linewidth',2.5,'markersize',8);
hold on
semilogy(min_eigen,MSE_R_eigen,'-s','color',[0.1882 0.2667 0.6275],'linewidth',2.5,'markersize',10);
semilogy(min_eigen,MSE_R_SDP,'-^','color',[0.3725 0.7843 0.2941],'linewidth',2.5,'markersize',8);
semilogy(min_eigen,MSE_R_GNE,'-x','color',[1.0000 0.6863 0],'linewidth',2.5,'markersize',10);
semilogy(min_eigen,MSE_R,'-d','color',[0.7255 0.2941 0.6863],'linewidth',2.5,'markersize',8);
semilogy(min_eigen,MSE_R_GN,'--d','color',[0.7255 0.2941 0.6863],'linewidth',2.5,'markersize',8);
semilogy(min_eigen,CRB_R,'-*','color',[0.3010 0.7450 0.9330],'linewidth',2.5,'markersize',10);
legend('LM','Eigen','SDP','GN-E','CECME1','CECME2','CRB');
grid on
% set(gca, 'XTick',min_eigen, 'XLim',[0 2]);

figure
semilogy(min_eigen,MSE_t_LM,'-o','color',[0.8627 0.1961 0.1373],'linewidth',2.5,'markersize',8);
hold on
semilogy(min_eigen,MSE_t_eigen,'-s','color',[0.1882 0.2667 0.6275],'linewidth',2.5,'markersize',10);
semilogy(min_eigen,MSE_t_SDP,'-^','color',[0.3725 0.7843 0.2941],'linewidth',2.5,'markersize',8);
semilogy(min_eigen,MSE_t_GNE,'-x','color',[1.0000 0.6863 0],'linewidth',2.5,'markersize',10);
semilogy(min_eigen,MSE_t,'-d','color',[0.7255 0.2941 0.6863],'linewidth',2.5,'markersize',8);
semilogy(min_eigen,MSE_t_GN,'--d','color',[0.7255 0.2941 0.6863],'linewidth',2.5,'markersize',8);
semilogy(min_eigen,CRB_t,'-*','color',[0.3010 0.7450 0.9330],'linewidth',2.5,'markersize',10);
legend('LM','Eigen','SDP','GN-E','CECME1','CECME2','CRB');
grid on
% set(gca, 'XTick',min_eigen, 'XLim',[0 2]);


