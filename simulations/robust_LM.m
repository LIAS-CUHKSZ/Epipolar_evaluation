function [E_LM,R_LM,t_LM]=robust_LM(z_h,y_h,R_init,t_init)

%%
% % TLS LM iteration
%     options = optimoptions('lsqnonlin', 'Display', 'off', ...
%                           'Algorithm', 'levenberg-marquardt');
%     % initial value
%     init_kexi = [0;0;0;t_init];
% 
%     % optimization
%     opt_kexi = lsqnonlin(@myfun, init_kexi, [], [], options);
%     
%     X_LM = opt_kexi(1:3);
%     t_LM = opt_kexi(4:6)/norm(opt_kexi(4:6));
%     Xhat = [0 -X_LM(3) X_LM(2); X_LM(3) 0 -X_LM(1); -X_LM(2) X_LM(1) 0];
%     R_LM = R_init * expm(Xhat);  
%     E_LM=cal_E(R_LM,t_LM);
% 
%     opt_cost=myfun(opt_kexi);
%     function g = myfun(kexi)
%         X = kexi(1:3);
%         t = kexi(4:6);
%         Xhat2 = [0 -X(3) X(2); X(3) 0 -X(1); -X(2) X(1) 0];
%         R = R_init * expm(Xhat2);
%         E=cal_E(R,t);
%         l=E*y_h;
%         d=sum(z_h.*l,1)./vecnorm(l(1:2,:));
%         delta=10^(-2);
%         g=truncate_loss(d,delta);
%     end

% TLS LM iteration
    options = optimoptions('lsqnonlin', 'Display', 'off', ...
                          'Algorithm', 'levenberg-marquardt');
    % initial value
    gamma0=asin(t_init(3));
if t_init(1)>0
    theta0=atan(t_init(2)/t_init(1));
else
    theta0=atan(t_init(2)/t_init(1))+pi;
end
    init_kexi = [0;0;0;theta0;gamma0];
    init_cost=myfun(init_kexi);

    % optimization
    opt_kexi = lsqnonlin(@myfun, init_kexi, [], [], options);
    
    X_LM = opt_kexi(1:3);
    t_LM = [cos(opt_kexi(5))*cos(opt_kexi(4)) cos(opt_kexi(5))*sin(opt_kexi(4)) sin(opt_kexi(5))]';
    Xhat = [0 -X_LM(3) X_LM(2); X_LM(3) 0 -X_LM(1); -X_LM(2) X_LM(1) 0];
    R_LM = R_init * expm(Xhat);  
    E_LM=cal_E(R_LM,t_LM);

    opt_cost=myfun(opt_kexi);
    function g = myfun(kexi)
        X = kexi(1:3);
        t = [cos(kexi(5))*cos(kexi(4)) cos(kexi(5))*sin(kexi(4)) sin(kexi(5))]';
        Xhat2 = [0 -X(3) X(2); X(3) 0 -X(1); -X(2) X(1) 0];
        R = R_init * expm(Xhat2);
        E=cal_E(R,t);
        l=E*y_h;
        d=sum(z_h.*l,1)./vecnorm(l(1:2,:));
        delta=5*10^(-3);
        g=truncate_loss(d,delta);

        % epipolar error on the other side
        l=E'*z_h;
        d=sum(y_h.*l,1)./vecnorm(l(1:2,:));
        delta=5*10^(-3);
        g=truncate_loss(d,delta);
    end



function L = truncate_loss(e, delta)
    abs_e = abs(e);
    L = min(abs_e, delta);
    L = L / delta;
end
end