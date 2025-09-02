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
if norm(true_t)==0
    true_t_normal=[1 1 1]'/sqrt(3);
end
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
min_dis=1*20; % min depth of point that cameras can see
K=[f_pix 0 u_pix;0 f_pix v_pix;0 0 1];
sigma_pix=1;
sigma=sigma_pix/f_pix;
m=[10 30 10^2 300 1000 3000];  % the number of points

%%
% initilization
MSE_t_GN_v9=zeros(length(sigma),length(m));
MSE_R_GN_v9=zeros(length(sigma),length(m));
MSE_t_GN_v8=zeros(length(sigma),length(m));
MSE_R_GN_v8=zeros(length(sigma),length(m));
MSE_t_GN_v7=zeros(length(sigma),length(m));
MSE_R_GN_v7=zeros(length(sigma),length(m));
MSE_t_GN_v6=zeros(length(sigma),length(m));
MSE_R_GN_v6=zeros(length(sigma),length(m));
CRB_R=zeros(length(sigma),length(m));
CRB_t=zeros(length(sigma),length(m));
MSE_sigma=zeros(length(sigma),length(m));

monte_carlo=1;  % monte carlo number for calculating the MSE

%%
% test begins
for l=1:length(sigma)
for j=1:length(m)
x=zeros(3,m(j)); % 3D points
y=zeros(2,m(j));  % 2D points in the first camera plane
z=zeros(2,m(j));  % 2D points in the second camera plane
z_noisefree=zeros(2,m(j));
for k=1:monte_carlo
%% 
% generate 2D and 3D points
counter=0;
while counter<m(j)
    y_temp=[rand*ix-ix/2 rand*iy-iy/2]';
    depth=min_dis+rand*(max_dis-min_dis);
    x_temp=[y_temp*depth/f ; depth];
    x_temp2=true_R*x_temp+true_t;
    if (x_temp2(3)>f) && (abs(x_temp2(1))<ix/2*x_temp2(3)/f) && (abs(x_temp2(2))<iy/2*x_temp2(3)/f)
        counter=counter+1;
        x(:,counter)=x_temp;
        y(:,counter)=y_temp;
        z_noisefree(:,counter)=x_temp2(1:2)/x_temp2(3)*f;
        z(:,counter)=x_temp2(1:2)/x_temp2(3)*f+randn(2,1)*sigma(l);
    end
end
y_h=[y;ones(1,m(j))];
z_h=[z;ones(1,m(j))];
z_noisefree_h=[z_noisefree;ones(1,m(j))];

%%
% pose estimation
[~,~,~,R_GN_v9,t_GN_v9,var_est] = consistent_essential_est(z_h,y_h,K);
[~,~,~,R_GN_v8,t_GN_v8,var_est] = consistent_essential_est_v8(z_h,y_h,K);
[~,~,~,R_GN_v7,t_GN_v7,var_est] = consistent_essential_est_v7(z_h,y_h,K);
[~,~,~,R_GN_v6,t_GN_v6,var_est] = consistent_essential_est_v6(z_h,y_h,K);


%%
% calculate MSE and bias
MSE_R_GN_v9(l,j)=MSE_R_GN_v9(l,j)+norm(R_GN_v9-true_R,'fro')^2/monte_carlo;
MSE_t_GN_v9(l,j)=MSE_t_GN_v9(l,j)+norm(t_GN_v9-true_t_normal)^2/monte_carlo;
MSE_R_GN_v8(l,j)=MSE_R_GN_v8(l,j)+norm(R_GN_v8-true_R,'fro')^2/monte_carlo;
MSE_t_GN_v8(l,j)=MSE_t_GN_v8(l,j)+norm(t_GN_v8-true_t_normal)^2/monte_carlo;
MSE_R_GN_v7(l,j)=MSE_R_GN_v7(l,j)+norm(R_GN_v7-true_R,'fro')^2/monte_carlo;
MSE_t_GN_v7(l,j)=MSE_t_GN_v7(l,j)+norm(t_GN_v7-true_t_normal)^2/monte_carlo;
MSE_R_GN_v6(l,j)=MSE_R_GN_v6(l,j)+norm(R_GN_v6-true_R,'fro')^2/monte_carlo;
MSE_t_GN_v6(l,j)=MSE_t_GN_v6(l,j)+norm(t_GN_v6-true_t_normal)^2/monte_carlo;
MSE_sigma(l,j)=MSE_sigma(l,j)+norm(var_est-sigma(l)^2)^2/monte_carlo;

end
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
for i=1:m(j)
    true_cov=sigma(l)^2*eye(2);
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
CRB_R(l,j)=trace(bar_F(1:9,1:9));
CRB_t(l,j)=trace(bar_F(10:12,10:12));
end
end

%%
% save results and plot a pair of figures
% save MSE_R_5pt MSE_R_5pt
% save MSE_t_5pt MSE_t_5pt
% save MSE_R_eigen MSE_R_eigen
% save MSE_t_eigen MSE_t_eigen
% save MSE_R_SDP MSE_R_SDP
% save MSE_t_SDP MSE_t_SDP
% save MSE_R_GNE MSE_R_GNE
% save MSE_t_GNE MSE_t_GNE
% save MSE_R MSE_R
% save MSE_t MSE_t
% save MSE_R_GN MSE_R_GN
% save MSE_t_GN MSE_t_GN
% save CRB_R CRB_R
% save CRB_t CRB_t
% save MSE_sigma MSE_sigma
% 
% save bias_R_5pt bias_R_5pt
% save bias_t_5pt bias_t_5pt
% save bias_R_eigen bias_R_eigen
% save bias_t_eigen bias_t_eigen
% save bias_R_SDP bias_R_SDP
% save bias_t_SDP bias_t_SDP
% save bias_R_GNE bias_R_GNE
% save bias_t_GNE bias_t_GNE
% save bias_R bias_R
% save bias_t bias_t
% save bias_R_GN bias_R_GN
% save bias_t_GN bias_t_GN
% save result_MSE_bias_vs_m_sigma MSE_R_5pt MSE_t_5pt MSE_R_eigen MSE_t_eigen MSE_R_SDP MSE_t_SDP MSE_R_GNE MSE_t_GNE MSE_R MSE_t MSE_R_GN MSE_t_GN CRB_R CRB_t MSE_sigma bias_R_5pt bias_t_5pt bias_R_eigen bias_t_eigen bias_R_SDP bias_t_SDP bias_R_GNE bias_t_GNE bias_R bias_t bias_R_GN bias_t_GN

loglog(m,CRB_R(1,:),'-*','color',[0.3010 0.7450 0.9330],'linewidth',2.5,'markersize',10);
hold on 
loglog(m,MSE_R_GN_v9(1,:),'--d','color',[0.7255 0.2941 0.6863],'linewidth',2.5,'markersize',8);
loglog(m,MSE_R_GN_v8(1,:),'-o','color',[0.8627 0.1961 0.1373],'linewidth',2.5,'markersize',8);
loglog(m,MSE_R_GN_v7(1,:),'-s','color',[0.1882 0.2667 0.6275],'linewidth',2.5,'markersize',10);
loglog(m,MSE_R_GN_v6(1,:),'-^','color',[0.3725 0.7843 0.2941],'linewidth',2.5,'markersize',8);
legend('CRB','v9','v8','v7','v6','fontsize',10);
grid on
set(gca, 'XTick',m, 'XLim',[10 3000]);
xlabel('Number of points','fontsize',12)
ylabel('Rotation MSE','fontsize',12)

figure
loglog(m,CRB_t(1,:),'-*','color',[0.3010 0.7450 0.9330],'linewidth',2.5,'markersize',10);
hold on
loglog(m,MSE_t_GN_v9(1,:),'--d','color',[0.7255 0.2941 0.6863],'linewidth',2.5,'markersize',8);
loglog(m,MSE_t_GN_v8(1,:),'-o','color',[0.8627 0.1961 0.1373],'linewidth',2.5,'markersize',8);
loglog(m,MSE_t_GN_v7(1,:),'-s','color',[0.1882 0.2667 0.6275],'linewidth',2.5,'markersize',10);
loglog(m,MSE_t_GN_v6(1,:),'-^','color',[0.3725 0.7843 0.2941],'linewidth',2.5,'markersize',8);
legend('CRB','v9','v8','v7','v6','fontsize',10);
grid on
set(gca, 'XTick',m, 'XLim',[10 3000]);
xlabel('Number of points','fontsize',12)
ylabel('Translation MSE','fontsize',12)



