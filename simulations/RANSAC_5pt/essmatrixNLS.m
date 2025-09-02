% ESSMATRIXNLS - refines essential matrix from N points
%
% Function refines the essential matrix from N matching points in
% a stereo pair of images. Non-linear least squares is used to minimize
% the average symmetric epipolar distance.
% Rotation is parameterized with the Rodrigues rotation vector
%
% Usage:   E = essmatrixNLS(E0, x1, x2)
%          E = essmatrixNLS(E0, x)
%
% Arguments:
%          E0     - Initial E estimate
%          x1, x2 - Two 3xN sets of corresponding homogeneous points
%         
%          x      - If a single argument is supplied it is assumed that it
%                   is in the form x = [x1; x2]
% Returns:
%          E  - A refined 3x3 essential matrix
%
% See Also: ESSMATRIX

% Requires marquardt.m

% Copyright (c) 2018 Manolis Lourakis
% Institute of Computer Science, Foundation for Research & Technology - Hellas
% Heraklion, Crete, Greece
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

% May 2018  - Original version.

function E = essmatrixNLS(varargin)
    
    [E0, x1, x2, npts] = checkargs(varargin(:));
    %Octave = exist('OCTAVE_VERSION', 'builtin') == 5; % Are we running under Octave

    if(npts < 6)
      error('essmatrixNLS: at least 6 point pairs expected');
    end
    
    % Normalise each set of points so that the origin 
    % is at centroid and mean distance from origin is sqrt(2). 
    % normalise2dpts also ensures the scale parameter is 1.
    %%%[x1, T1] = normalise2dpts(x1);
    %%%[x2, T2] = normalise2dpts(x2);
    
    % make point sets inhomogeneous: divide by third coord if non-unit
    if(sum(x1(3,:)~=1)>0), x1=x1(1:2,:)./repmat(x1(3,:), 2,1); else, x1=x1(1:2,:); end
    if(sum(x2(3,:)~=1)>0), x2=x2(1:2,:)./repmat(x2(3,:), 2,1); else, x2=x2(1:2,:); end

    rt=E2rt(E0);
    [rt info]=marquardt(@objfunc, {x1, x2}, rt, [1E-04 1E-14 1E-14 50]);
%     fprintf('\nessmatrixNLS: LM terminated after %d iterations, reason %d\n', info(5), info(6));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
    % Gauss-Newton
    e=zeros(npts, 1);
    J=zeros(npts, 6);
    prevnrm=inf;
    for j=1:20
      [e J]=objfunc(rt, {x1, x2});
      nrm=norm(e);
      if(1-nrm/prevnrm<1000*eps)
        break;
      end
      prevnrm=nrm;

      % d=(J'*J)\(J'*e); % normal eqs
      [Q, R]=qr(J); d=R\(Q'*e); % QR
      rt=rt - d;
    end
    fprintf('\nessmatrixNLS: GN terminated after %d iterations\n', j);
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    E=rt2E(rt);

    % Denormalise
    %%%E = T2'*E*T1;
    
%--------------------------------------------------------------------------
function [e, J] = objfunc(rt, parms)
  x1=parms{1};
  x2=parms{2};

  E=rt2E(rt);
  Ev=reshape(E.',1,[]); % flatten
  E_jac=calcEJac_rt(rt);

  npts=size(x1, 2);
  e=zeros(npts, 1);
  J=zeros(npts, 6);

  for i=1:npts
    p1=x1(:,i);
    p2=x2(:,i);
    e(i)=calcEpipDist(p1, p2, rt);
    %J(i,:)=calcEpipDistJac(p1, p2, rt); % direct calculation
    dist_grad=calcEpipDistJac_E(p1, p2, Ev);
    J(i,:)=dist_grad*E_jac; % chain rule: d dist(E) / d r,t = d dist(E)/ d E|_[E=rt2E(r,t)] * d E / d r,t
  end

%--------------------------------------------------------------------------
% finite difference Jacobian for verification
function J = fjac(rt, parms)
  x1=parms{1};
  x2=parms{2};

  delta=1e-6*(max(1,sqrt(norm(rt))));
  npts=size(x1, 2);
  m=6;
  J=zeros(npts, m);
  drt=zeros(m, 1);
  for i=1:npts
    p1=x1(:,i);
    p2=x2(:,i);
    e=calcEpipDist(p1, p2, rt);
    for j=1:m
      drt(j)=delta;
      J(i,j)=(calcEpipDist(p1, p2, rt+drt) - e)/delta; % fwd
      %J(i,j)=(calcEpipDist(p1, p2, rt+drt) - calcEpipDist(p1, p2, rt-drt))/(2*delta); % central
      drt(j)=0;
    end
  end

%--------------------------------------------------------------------------
% Function to check argument values and set defaults

function [E0, x1, x2, npts] = checkargs(arg);
    
    if length(arg) == 3
        E0 = arg{1};
        x1 = arg{2};
        x2 = arg{3};
        if ~all(size(x1)==size(x2))
            error('x1 and x2 must have the same size');
        elseif size(x1,1) ~= 3
            error('x1 and x2 must be 3xN');
        end
        
    elseif length(arg) == 2
        if size(arg{2},1) ~= 6
            error('Single argument x must be 6xN');
        else
            E0 = arg{1};
            x1 = arg{2}(1:3,:);
            x2 = arg{2}(4:6,:);
        end
    else
        error('Wrong number of arguments supplied');
    end
      
    npts = size(x1,2);
    if npts < 5
        error('At least 5 points are needed to compute the essential matrix');
    end
    
%--------------------------------------------------------------------------
% extract a 6 vector rt from E
function [rt] = E2rt(E)
  D=[0 1 0; -1 0 0; 0 0 1];
  [U,dummy,V]=svd(E);
  R=U*D*V'; % Ra
  t=U(:,3); % tu
  if(det(R)<0)
    R=-R;
    t=-t;
  end

  % convert R to rotation vector; adapted from rodrigues.m

  % project the rotation matrix to SO(3);
  %[U,S,V] = svd(R);
  %R = U*V';

  tr = (trace(R)-1)/2;
  theta = real(acos(tr));

  if sin(theta) >= 1E-5,
    vth = 1/(2*sin(theta));
    om1 = [R(3,2)-R(2,3), R(1,3)-R(3,1), R(2,1)-R(1,2)]';
    om = vth*om1;
    rv = om*theta;
  else
    if tr > 0;       % case norm(om)=0
      rv = [0 0 0]';
    else         % case norm(om)=pi
      rv = theta * (sqrt((diag(R)+1)/2).*[1;2*(R(1,2:3)>=0)'-1]);
    end
  end

  rt=[rv; t];

%--------------------------------------------------------------------------
% convert a 6 vector made of r, t to E (i.e., [t]x*R)
function [E] = rt2E(rt)
  t1 = rt(1) ^ 2;
  t2 = rt(2) ^ 2;
  t3 = rt(3) ^ 2;
  t4 = t1 + t2 + t3;
  t5 = sqrt(t4);
  t6 = sin(t5);
  t8 = t6 / t5;
  t9 = t8 * rt(3);
  t10 = cos(t5);
  t11 = 0.1e1 - t10;
  t12 = 0.1e1 / t4;
  t13 = t11 * t12;
  t15 = t13 * rt(2) * rt(1);
  t16 = t9 + t15;
  t18 = t8 * rt(2);
  t20 = t13 * rt(3) * rt(1);
  t21 = -t18 + t20;
  t24 = t12 * t3;
  t25 = t12 * t1;
  t28 = 0.1e1 + t11 * (-t24 - t25);
  t30 = t8 * rt(1);
  t32 = t13 * rt(3) * rt(2);
  t33 = t30 + t32;
  t36 = -t30 + t32;
  t38 = t12 * t2;
  t41 = 0.1e1 + t11 * (-t38 - t25);
  t46 = 0.1e1 + t11 * (-t24 - t38);
  t50 = -t9 + t15;
  t54 = t18 + t20;

  E= [-rt(6) * t16 + rt(5) * t21, -rt(6) * t28 + rt(5) * t33, -rt(6) * t36 + rt(5) * t41;...
       rt(6) * t46 - rt(4) * t21,  rt(6) * t50 - rt(4) * t33,  rt(6) * t54 - rt(4) * t41;...
      -rt(5) * t46 + rt(4) * t16, -rt(5) * t50 + rt(4) * t28, -rt(5) * t54 + rt(4) * t36];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Code automatically generated by maple

% epipolar distance from rt
function d = calcEpipDist(m1, m2, rt)
  t1 = m2(1);
  t2 = rt(6);
  t3 = rt(1);
  t4 = t3 ^ 2;
  t5 = rt(2);
  t6 = t5 ^ 2;
  t7 = rt(3);
  t8 = t7 ^ 2;
  t9 = t4 + t6 + t8;
  t10 = sqrt(t9);
  t11 = sin(t10);
  t13 = t11 / t10;
  t14 = t13 * t7;
  t15 = cos(t10);
  t16 = 0.1e1 - t15;
  t17 = 0.1e1 / t9;
  t18 = t16 * t17;
  t20 = t18 * t5 * t3;
  t21 = t14 + t20;
  t23 = rt(5);
  t24 = t13 * t5;
  t22 = t18 * t7;
  t26 = t22 * t3;
  t27 = -t24 + t26;
  t29 = -t2 * t21 + t23 * t27;
  t30 = m1(1);
  t32 = t17 * t8;
  t33 = t17 * t4;
  t36 = 0.1e1 + t16 * (-t32 - t33);
  t38 = t13 * t3;
  t40 = t22 * t5;
  t41 = t38 + t40;
  t43 = -t2 * t36 + t23 * t41;
  t44 = m1(2);
  t46 = -t38 + t40;
  t47 = t2 * t46;
  t48 = t17 * t6;
  t51 = 0.1e1 + t16 * (-t48 - t33);
  t52 = t23 * t51;
  t53 = t29 * t30 + t43 * t44 - t47 + t52;
  t55 = m2(2);
  t58 = 0.1e1 + t16 * (-t32 - t48);
  t60 = rt(4);
  t62 = t2 * t58 - t60 * t27;
  t64 = -t14 + t20;
  t67 = t2 * t64 - t60 * t41;
  t69 = t24 + t26;
  t70 = t2 * t69;
  t71 = t60 * t51;
  t72 = t62 * t30 + t67 * t44 + t70 - t71;
  t74 = t23 * t58;
  t75 = t60 * t21;
  t78 = t23 * t64;
  t79 = t60 * t36;
  t82 = t23 * t69;
  t83 = t60 * t46;
  t85 = (t53 * t1 + t55 * t72 + (-t74 + t75) * t30 + (-t78 + t79) * t44 - t82 + t83) ^ 2;
  t86 = t53 ^ 2;
  t87 = t72 ^ 2;
  t93 = t29 * t1 + t62 * t55 - t74 + t75;
  t97 = t43 * t1 + t67 * t55 - t78 + t79;
  t104 = (t30 * t93 + t44 * t97 + (-t47 + t52) * t1 + (t70 - t71) * t55 - t82 + t83) ^ 2;
  t105 = t93 ^ 2;
  t106 = t97 ^ 2;
  d = t85 / (t86 + t87) + t104 / (t105 + t106);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% gradient for epipolar distance from rt
function d_grad = calcEpipDistJac(m1, m2, rt)
  t1 = m2(1);
  t2 = rt(6);
  t3 = rt(1);
  t4 = t3 ^ 2;
  t5 = rt(2);
  t6 = t5 ^ 2;
  t7 = rt(3);
  t8 = t7 ^ 2;
  t9 = t4 + t6 + t8;
  t10 = sqrt(t9);
  t11 = sin(t10);
  t12 = 0.1e1 / t10;
  t13 = t11 * t12;
  t14 = t13 * t7;
  t15 = cos(t10);
  t16 = 0.1e1 - t15;
  t17 = 0.1e1 / t9;
  t18 = t16 * t17;
  t19 = t5 * t3;
  t20 = t18 * t19;
  t21 = t14 + t20;
  t23 = rt(5);
  t24 = t13 * t5;
  t25 = t7 * t3;
  t26 = t18 * t25;
  t27 = -t24 + t26;
  t29 = -t2 * t21 + t23 * t27;
  t30 = m1(1);
  t32 = t17 * t8;
  t33 = t17 * t4;
  t34 = -t32 - t33;
  t35 = t16 * t34;
  t36 = 0.1e1 + t35;
  t38 = t13 * t3;
  t39 = t7 * t5;
  t40 = t18 * t39;
  t41 = t38 + t40;
  t43 = -t2 * t36 + t23 * t41;
  t44 = m1(2);
  t46 = -t38 + t40;
  t47 = t2 * t46;
  t48 = t17 * t6;
  t49 = -t48 - t33;
  t50 = t16 * t49;
  t51 = 0.1e1 + t50;
  t52 = t23 * t51;
  t53 = t29 * t30 + t43 * t44 - t47 + t52;
  t55 = m2(2);
  t56 = -t32 - t48;
  t57 = t16 * t56;
  t58 = 0.1e1 + t57;
  t60 = rt(4);
  t62 = t2 * t58 - t60 * t27;
  t64 = -t14 + t20;
  t67 = t2 * t64 - t60 * t41;
  t69 = t24 + t26;
  t70 = t2 * t69;
  t71 = t60 * t51;
  t72 = t62 * t30 + t67 * t44 + t70 - t71;
  t74 = t23 * t58;
  t75 = t60 * t21;
  t78 = t23 * t64;
  t79 = t60 * t36;
  t82 = t23 * t69;
  t83 = t60 * t46;
  t84 = t53 * t1 + t55 * t72 + (-t74 + t75) * t30 + (-t78 + t79) * t44 - t82 + t83;
  t85 = t53 ^ 2;
  t86 = t72 ^ 2;
  t87 = t85 + t86;
  t89 = t84 / t87;
  t90 = t15 * t17;
  t91 = t90 * t25;
  t94 = t11 * t12 * t17;
  t95 = t94 * t25;
  t96 = t4 * t5;
  t97 = t94 * t96;
  t98 = t9 ^ 2;
  t99 = 0.1e1 / t98;
  t100 = t16 * t99;
  t102 = 0.2e1 * t100 * t96;
  t103 = t18 * t5;
  t104 = t91 - t95 + t97 - t102 + t103;
  t106 = t90 * t19;
  t107 = t94 * t19;
  t108 = t4 * t7;
  t109 = t94 * t108;
  t111 = 0.2e1 * t100 * t108;
  t112 = t18 * t7;
  t113 = -t106 + t107 + t109 - t111 + t112;
  t115 = -t2 * t104 + t23 * t113;
  t119 = t99 * t8;
  t120 = t119 * t3;
  t122 = t99 * t4 * t3;
  t123 = t17 * t3;
  t116 = t13 * t3;
  t127 = t116 * t34 + t16 * (0.2e1 * t120 + 0.2e1 * t122 - 0.2e1 * t123);
  t129 = t90 * t4;
  t130 = t94 * t4;
  t131 = t25 * t5;
  t132 = t94 * t131;
  t134 = 0.2e1 * t100 * t131;
  t135 = t129 - t130 + t13 + t132 - t134;
  t137 = -t2 * t127 + t23 * t135;
  t139 = -t129 + t130 - t13 + t132 - t134;
  t140 = t2 * t139;
  t143 = t99 * t6;
  t144 = t143 * t3;
  t148 = t116 * t49 + t16 * (0.2e1 * t144 + 0.2e1 * t122 - 0.2e1 * t123);
  t149 = t23 * t148;
  t150 = t30 * t115 + t137 * t44 - t140 + t149;
  t157 = t116 * t56 + t16 * (0.2e1 * t120 + 0.2e1 * t144);
  t160 = t2 * t157 - t60 * t113;
  t162 = -t91 + t95 + t97 - t102 + t103;
  t165 = t2 * t162 - t60 * t135;
  t167 = t106 - t107 + t109 - t111 + t112;
  t168 = t2 * t167;
  t169 = t60 * t148;
  t170 = t160 * t30 + t165 * t44 + t168 - t169;
  t172 = t23 * t157;
  t173 = t60 * t104;
  t176 = t23 * t162;
  t177 = t60 * t127;
  t180 = t23 * t167;
  t181 = t60 * t139;
  t185 = t84 ^ 2;
  t186 = t87 ^ 2;
  t188 = t185 / t186;
  t196 = t29 * t1 + t62 * t55 - t74 + t75;
  t200 = t43 * t1 + t67 * t55 - t78 + t79;
  t206 = t30 * t196 + t44 * t200 + (-t47 + t52) * t1 + (t70 - t71) * t55 - t82 + t83;
  t207 = t196 ^ 2;
  t208 = t200 ^ 2;
  t209 = t207 + t208;
  t211 = t206 / t209;
  t214 = t115 * t1 + t160 * t55 - t172 + t173;
  t218 = t137 * t1 + t165 * t55 - t176 + t177;
  t227 = t206 ^ 2;
  t228 = t209 ^ 2;
  t230 = t227 / t228;
  d_grad(1) = 0.2e1 * t89 * (t1 * t150 + t55 * t170 + (-t172 + t173) * t30 + (-t176 + t177) * t44 - t180 + t181) - t188 * (0.2e1 * t53 * t150 + 0.2e1 * t72 * t170) + 0.2e1 * t211 * (t30 * t214 + t44 * t218 + (-t140 + t149) * t1 + (t168 - t169) * t55 - t180 + t181) - t230 * (0.2e1 * t196 * t214 + 0.2e1 * t200 * t218);
  t236 = t90 * t39;
  t237 = t94 * t39;
  t238 = t6 * t3;
  t239 = t94 * t238;
  t241 = 0.2e1 * t100 * t238;
  t242 = t18 * t3;
  t243 = t236 - t237 + t239 - t241 + t242;
  t245 = t90 * t6;
  t246 = t94 * t6;
  t247 = -t245 + t246 - t13 + t132 - t134;
  t249 = -t2 * t243 + t23 * t247;
  t253 = t119 * t5;
  t254 = t99 * t4;
  t255 = t254 * t5;
  t240 = t13 * t5;
  t259 = t240 * t34 + t16 * (0.2e1 * t253 + 0.2e1 * t255);
  t261 = t6 * t7;
  t262 = t94 * t261;
  t264 = 0.2e1 * t100 * t261;
  t265 = t106 - t107 + t262 - t264 + t112;
  t267 = -t2 * t259 + t23 * t265;
  t269 = -t106 + t107 + t262 - t264 + t112;
  t270 = t2 * t269;
  t274 = t99 * t6 * t5;
  t275 = t17 * t5;
  t279 = t240 * t49 + t16 * (0.2e1 * t274 - 0.2e1 * t275 + 0.2e1 * t255);
  t280 = t23 * t279;
  t281 = t249 * t30 + t267 * t44 - t270 + t280;
  t288 = t240 * t56 + t16 * (0.2e1 * t253 + 0.2e1 * t274 - 0.2e1 * t275);
  t291 = t2 * t288 - t60 * t247;
  t293 = -t236 + t237 + t239 - t241 + t242;
  t296 = t2 * t293 - t60 * t265;
  t298 = t245 - t246 + t13 + t132 - t134;
  t299 = t2 * t298;
  t300 = t60 * t279;
  t301 = t291 * t30 + t296 * t44 + t299 - t300;
  t303 = t23 * t288;
  t304 = t60 * t243;
  t307 = t23 * t293;
  t308 = t60 * t259;
  t311 = t23 * t298;
  t312 = t60 * t269;
  t323 = t249 * t1 + t291 * t55 - t303 + t304;
  t327 = t267 * t1 + t296 * t55 - t307 + t308;
  d_grad(2) = 0.2e1 * t89 * (t1 * t281 + t55 * t301 + (-t303 + t304) * t30 + (-t307 + t308) * t44 - t311 + t312) - t188 * (0.2e1 * t53 * t281 + 0.2e1 * t72 * t301) + 0.2e1 * t211 * (t30 * t323 + t44 * t327 + (-t270 + t280) * t1 + (t299 - t300) * t55 - t311 + t312) - t230 * (0.2e1 * t196 * t323 + 0.2e1 * t200 * t327);
  t341 = t90 * t8;
  t342 = t94 * t8;
  t343 = t341 - t342 + t13 + t132 - t134;
  t345 = t8 * t3;
  t346 = t94 * t345;
  t348 = 0.2e1 * t100 * t345;
  t349 = -t236 + t237 + t346 - t348 + t242;
  t351 = -t2 * t343 + t23 * t349;
  t356 = t99 * t8 * t7;
  t357 = t17 * t7;
  t358 = t254 * t7;
  t340 = t13 * t7;
  t362 = t340 * t34 + t16 * (0.2e1 * t356 - 0.2e1 * t357 + 0.2e1 * t358);
  t364 = t8 * t5;
  t365 = t94 * t364;
  t367 = 0.2e1 * t100 * t364;
  t368 = t91 - t95 + t365 - t367 + t103;
  t370 = -t2 * t362 + t23 * t368;
  t372 = -t91 + t95 + t365 - t367 + t103;
  t373 = t2 * t372;
  t376 = t143 * t7;
  t380 = t340 * t49 + t16 * (0.2e1 * t376 + 0.2e1 * t358);
  t381 = t23 * t380;
  t382 = t351 * t30 + t370 * t44 - t373 + t381;
  t389 = t340 * t56 + t16 * (0.2e1 * t356 - 0.2e1 * t357 + 0.2e1 * t376);
  t392 = t2 * t389 - t60 * t349;
  t394 = -t341 + t342 - t13 + t132 - t134;
  t397 = t2 * t394 - t60 * t368;
  t399 = t236 - t237 + t346 - t348 + t242;
  t400 = t2 * t399;
  t401 = t60 * t380;
  t402 = t392 * t30 + t397 * t44 + t400 - t401;
  t404 = t23 * t389;
  t405 = t60 * t343;
  t408 = t23 * t394;
  t409 = t60 * t362;
  t412 = t23 * t399;
  t413 = t60 * t372;
  t424 = t351 * t1 + t392 * t55 - t404 + t405;
  t428 = t370 * t1 + t397 * t55 - t408 + t409;
  d_grad(3) = 0.2e1 * t89 * (t1 * t382 + t55 * t402 + (-t404 + t405) * t30 + (-t408 + t409) * t44 - t412 + t413) - t188 * (0.2e1 * t53 * t382 + 0.2e1 * t72 * t402) + 0.2e1 * t211 * (t30 * t424 + t44 * t428 + (-t373 + t381) * t1 + (t400 - t401) * t55 - t412 + t413) - t230 * (0.2e1 * t196 * t424 + 0.2e1 * t200 * t428);
  t442 = -t27;
  t444 = -t41;
  t446 = -0.1e1 + t442 * t30 + t444 * t44 - t50;
  t457 = t442 * t55 + t14 + t20;
  t460 = 0.1e1 + t444 * t55 + t35;
  d_grad(4) = 0.2e1 * t89 * (t55 * t446 + t21 * t30 + t36 * t44 - t38 + t40) - 0.2e1 * t188 * t72 * t446 + 0.2e1 * t211 * (t30 * t457 + t44 * t460 - t51 * t55 - t38 + t40) - t230 * (0.2e1 * t196 * t457 + 0.2e1 * t200 * t460);
  t474 = 0.1e1 + t27 * t30 + t41 * t44 + t50;
  t487 = -0.1e1 + t27 * t1 - t57;
  t490 = t41 * t1 + t14 - t20;
  t469 = t58 * t30;
  t470 = t64 * t44;
  d_grad(5) = 0.2e1 * t89 * (t1 * t474 - t469 - t470 - t24 - t26) - 0.2e1 * t188 * t53 * t474 + 0.2e1 * t211 * (t30 * t487 + t44 * t490 + t51 * t1 - t24 - t26) - t230 * (0.2e1 * t196 * t487 + 0.2e1 * t200 * t490);
  t501 = -t21;
  t503 = -t36;
  t505 = t30 * t501 + t503 * t44 + t38 - t40;
  t509 = t469 + t470 + t24 + t26;
  t521 = t501 * t1 + t58 * t55;
  t525 = t503 * t1 + t64 * t55;
  d_grad(6) = 0.2e1 * t89 * (t1 * t505 + t55 * t509) - t188 * (0.2e1 * t53 * t505 + 0.2e1 * t72 * t509) + 0.2e1 * t211 * (t30 * t521 + t44 * t525 - t46 * t1 + t69 * t55) - t230 * (0.2e1 * t196 * t521 + 0.2e1 * t200 * t525);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% epipolar distance from rt gradient using the chain rule (and the jacobian of E)
function d_grad = calcEpipDistJac_E(m1, m2, Evec)
  t1 = m2(1);
  t2 = Evec(1);
  t3 = m1(1);
  t5 = Evec(2);
  t6 = m1(2);
  t8 = Evec(3);
  t9 = t2 * t3 + t5 * t6 + t8;
  t11 = m2(2);
  t12 = Evec(4);
  t14 = Evec(5);
  t16 = Evec(6);
  t17 = t12 * t3 + t14 * t6 + t16;
  t19 = Evec(7);
  t21 = Evec(8);
  t23 = Evec(9);
  t24 = t1 * t9 + t11 * t17 + t19 * t3 + t21 * t6 + t23;
  t25 = t9 ^ 2;
  t26 = t17 ^ 2;
  t27 = t25 + t26;
  t29 = (t24 / t27);
  t30 = (t1 * t3);
  t32 = t24 ^ 2;
  t33 = t27 ^ 2;
  t35 = t32 / t33;
  t40 = t2 * t1 + t12 * t11 + t19;
  t44 = t5 * t1 + t14 * t11 + t21;
  t48 = t3 * t40 + t44 * t6 + t8 * t1 + t16 * t11 + t23;
  t49 = t40 ^ 2;
  t50 = t44 ^ 2;
  t51 = t49 + t50;
  t53 = (t48 / t51);
  t55 = t48 ^ 2;
  t56 = t51 ^ 2;
  t58 = t55 / t56;
  t52 = (t35 * t9);
  t59 = (t58 * t40);
  d_grad(1) = 2 * t29 * t30 - 2 * t52 * t3 + 2 * t53 * t30 - 2 * t59 * t1;
  t62 = t1 * t6;
  t66 = (t58 * t44);
  d_grad(2) = 2 * t29 * t62 - 2 * t52 * t6 + 2 * t53 * t62 - 2 * t66 * t1;
  d_grad(3) = 2 * t29 * t1 - 2 * t52 + 2 * t53 * t1;
  t76 = t11 * t3;
  t73 = (t35 * t17);
  d_grad(4) = 2 * t29 * t76 - 2 * t73 * t3 + 2 * t53 * t76 - 2 * t59 * t11;
  t84 = t11 * t6;
  d_grad(5) = 2 * t29 * t84 - 2 * t73 * t6 + 2 * t53 * t84 - 2 * t66 * t11;
  d_grad(6) = 2 * t29 * t11 - 2 * t73 + 2 * t53 * t11;
  d_grad(7) = 2 * t29 * t3 + 2 * t53 * t3 - 2 * t59;
  d_grad(8) = 2 * t29 * t6 + 2 * t53 * t6 - 2 * t66;
  d_grad(9) = 2 * t29 + 2 * t53;
  return;

% jacobian of E with respect to rt
function E_jac = calcEJac_rt(rt)
  t1 = rt(6);
  t2 = rt(1);
  t3 = t2 ^ 2;
  t4 = rt(2);
  t5 = t4 ^ 2;
  t6 = rt(3);
  t7 = t6 ^ 2;
  t8 = t3 + t5 + t7;
  t9 = sqrt(t8);
  t10 = cos(t9);
  t11 = 0.1e1 / t8;
  t12 = t10 * t11;
  t13 = t2 * t6;
  t14 = t13 * t12;
  t15 = sin(t9);
  t16 = 0.1e1 / t9;
  t18 = t15 * t16 * t11;
  t19 = t18 * t13;
  t20 = t3 * t4;
  t21 = t20 * t18;
  t22 = 0.1e1 - t10;
  t23 = t8 ^ 2;
  t24 = 0.1e1 / t23;
  t25 = t22 * t24;
  t27 = 0.2e1 * t25 * t20;
  t28 = t22 * t11;
  t29 = t4 * t28;
  t30 = t14 - t19 + t21 - t27 + t29;
  t32 = rt(5);
  t33 = t2 * t4;
  t34 = t12 * t33;
  t35 = t18 * t33;
  t36 = t6 * t3;
  t37 = t18 * t36;
  t39 = 0.2e1 * t36 * t25;
  t40 = t28 * t6;
  t41 = -t34 + t35 + t37 - t39 + t40;
  E_jac(1,1) = -t1 * t30 + t41 * t32;
  t44 = t16 * t15;
  t45 = t11 * t7;
  t46 = t11 * t3;
  t47 = -t45 - t46;
  t50 = t7 * t24;
  t51 = t50 * t2;
  t53 = t24 * t3 * t2;
  t54 = t11 * t2;
  t48 = t44 * t2;
  t58 = t48 * t47 + t22 * (0.2e1 * t51 + 0.2e1 * t53 - 0.2e1 * t54);
  t60 = t12 * t3;
  t61 = t3 * t18;
  t62 = t13 * t4;
  t63 = t18 * t62;
  t65 = 0.2e1 * t25 * t62;
  t66 = t60 - t61 + t44 + t63 - t65;
  E_jac(2,1) = -t1 * t58 + t32 * t66;
  t68 = -t60 + t61 - t44 + t63 - t65;
  t70 = t11 * t5;
  t71 = -t70 - t46;
  t74 = t24 * t5;
  t75 = t74 * t2;
  t79 = t48 * t71 + t22 * (0.2e1 * t75 + 0.2e1 * t53 - 0.2e1 * t54);
  E_jac(3,1) = -t1 * t68 + t32 * t79;
  t81 = -t45 - t70;
  t87 = t48 * t81 + t22 * (0.2e1 * t51 + 0.2e1 * t75);
  t89 = rt(4);
  E_jac(4,1) = t1 * t87 - t89 * t41;
  t91 = -t14 + t19 + t21 - t27 + t29;
  E_jac(5,1) = t1 * t91 - t89 * t66;
  t94 = t34 - t35 + t37 - t39 + t40;
  E_jac(6,1) = t1 * t94 - t89 * t79;
  E_jac(7,1) = -t32 * t87 + t89 * t30;
  E_jac(8,1) = -t32 * t91 + t89 * t58;
  E_jac(9,1) = -t32 * t94 + t89 * t68;
  t103 = t4 * t6;
  t104 = t12 * t103;
  t105 = t103 * t18;
  t106 = t5 * t2;
  t107 = t18 * t106;
  t109 = 0.2e1 * t25 * t106;
  t110 = t28 * t2;
  t111 = t104 - t105 + t107 - t109 + t110;
  t113 = t12 * t5;
  t114 = t18 * t5;
  t115 = -t113 + t114 - t44 + t63 - t65;
  E_jac(1,2) = -t1 * t111 + t32 * t115;
  t119 = t50 * t4;
  t120 = t24 * t3;
  t121 = t120 * t4;
  t108 = t44 * t4;
  t125 = t108 * t47 + t22 * (0.2e1 * t119 + 0.2e1 * t121);
  t127 = t5 * t6;
  t128 = t18 * t127;
  t130 = 0.2e1 * t25 * t127;
  t131 = t34 - t35 + t128 - t130 + t40;
  E_jac(2,2) = -t1 * t125 + t32 * t131;
  t133 = -t34 + t35 + t128 - t130 + t40;
  t138 = t24 * t5 * t4;
  t139 = t11 * t4;
  t143 = t108 * t71 + t22 * (0.2e1 * t138 - 0.2e1 * t139 + 0.2e1 * t121);
  E_jac(3,2) = -t1 * t133 + t32 * t143;
  t150 = t108 * t81 + t22 * (0.2e1 * t119 + 0.2e1 * t138 - 0.2e1 * t139);
  E_jac(4,2) = t1 * t150 - t89 * t115;
  t153 = -t104 + t105 + t107 - t109 + t110;
  E_jac(5,2) = t1 * t153 - t89 * t131;
  t156 = t113 - t114 + t44 + t63 - t65;
  E_jac(6,2) = t1 * t156 - t89 * t143;
  E_jac(7,2) = -t32 * t150 + t89 * t111;
  E_jac(8,2) = -t32 * t153 + t89 * t125;
  E_jac(9,2) = -t32 * t156 + t89 * t133;
  t165 = t12 * t7;
  t166 = t18 * t7;
  t167 = t165 - t166 + t44 + t63 - t65;
  t169 = t7 * t2;
  t170 = t18 * t169;
  t172 = 0.2e1 * t25 * t169;
  t173 = -t104 + t105 + t170 - t172 + t110;
  E_jac(1,3) = -t1 * t167 + t32 * t173;
  t178 = t24 * t7 * t6;
  t179 = t11 * t6;
  t180 = t6 * t120;
  t164 = t44 * t6;
  t184 = t164 * t47 + t22 * (0.2e1 * t178 - 0.2e1 * t179 + 0.2e1 * t180);
  t186 = t4 * t7;
  t187 = t18 * t186;
  t189 = 0.2e1 * t25 * t186;
  t190 = t14 - t19 + t187 - t189 + t29;
  E_jac(2,3) = -t1 * t184 + t32 * t190;
  t192 = -t14 + t19 + t187 - t189 + t29;
  t196 = t74 * t6;
  t200 = t164 * t71 + t22 * (0.2e1 * t196 + 0.2e1 * t180);
  E_jac(3,3) = -t1 * t192 + t32 * t200;
  t207 = t164 * t81 + t22 * (0.2e1 * t178 - 0.2e1 * t179 + 0.2e1 * t196);
  E_jac(4,3) = t1 * t207 - t89 * t173;
  t210 = -t165 + t166 - t44 + t63 - t65;
  E_jac(5,3) = t1 * t210 - t89 * t190;
  t213 = t104 - t105 + t170 - t172 + t110;
  E_jac(6,3) = t1 * t213 - t200 * t89;
  E_jac(7,3) = -t32 * t207 + t89 * t167;
  E_jac(8,3) = -t32 * t210 + t89 * t184;
  E_jac(9,3) = -t32 * t213 + t89 * t192;
  E_jac(1,4) = 0.0e0;
  E_jac(2,4) = 0.0e0;
  E_jac(3,4) = 0.0e0;
  t222 = t44 * t4;
  t223 = t28 * t13;
  E_jac(4,4) = t222 - t223;
  t224 = t44 * t2;
  t225 = t28 * t103;
  E_jac(5,4) = -t224 - t225;
  E_jac(6,4) = -0.1e1 - t22 * t71;
  t227 = t44 * t6;
  t228 = t28 * t33;
  E_jac(7,4) = t227 + t228;
  E_jac(8,4) = 0.1e1 + t22 * t47;
  E_jac(9,4) = -t224 + t225;
  E_jac(1,5) = -E_jac(4,4);
  E_jac(2,5) = -E_jac(5,4);
  E_jac(3,5) = -E_jac(6,4);
  E_jac(4,5) = 0.0e0;
  E_jac(5,5) = 0.0e0;
  E_jac(6,5) = 0.0e0;
  E_jac(7,5) = -0.1e1 - t22 * t81;
  E_jac(8,5) = t227 - t228;
  E_jac(9,5) = -t222 - t223;
  E_jac(1,6) = -E_jac(7,4);
  E_jac(2,6) = -E_jac(8,4);
  E_jac(3,6) = -E_jac(9,4);
  E_jac(4,6) = -E_jac(7,5);
  E_jac(5,6) = -E_jac(8,5);
  E_jac(6,6) = -E_jac(9,5);
  E_jac(7,6) = 0.0e0;
  E_jac(8,6) = 0.0e0;
  E_jac(9,6) = 0.0e0;
  return;
