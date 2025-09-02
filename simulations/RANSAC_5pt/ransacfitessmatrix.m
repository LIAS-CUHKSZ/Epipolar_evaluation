% RANSACFITESSMATRIX - fits essential matrix using RANSAC
%
% Usage:   [E, inliers] = ransacfitessmatrix(x1, x2, t)
%
% Arguments:
%          x1  - 2xN or 3xN set of homogeneous points.  If the data is
%                2xN it is assumed the homogeneous scale factor is 1.
%          x2  - 2xN or 3xN set of homogeneous points such that x1<->x2.
%                Note that x1, x2 are assumed normalized, i.e. the effect
%                of the camera intrinsics K has been removed.
%          t   - The distance threshold between data point and the model
%                used to decide whether a point is an inlier or not. 
%                Note that this distance applies to normalized point
%                coordinates, not pixels. The value of t should reflect this,
%                say in the range 0.0001 - 0.01  
%                If -1<t<0, then -t is interpreted as the outliers fraction
%                and the threshold is adaptively computed from the corresponding
%                percentile of all distances.
%
% Note that it is assumed that the matching of x1 and x2 are putative and it
% is expected that a percentage of matches will be wrong.
%
% Returns:
%          E       - The 3x3 essential matrix such that x2'Ex1 = 0.
%          inliers - An array of indices of the elements of x1, x2 that were
%                    the inliers for the best model.
%
% See Also: RANSAC, ESSMATRIX

% Copyright (c) 2018 Manolis Lourakis
% Institute of Computer Science, Foundation for Research & Technology - Hellas
% Heraklion, Crete, Greece
%
% Based on code from ransacfitfundmatrix.m
% Copyright (c) 2004-2005 Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% http://www.csse.uwa.edu.au/
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

% May 2018  Original version
% Jan 2020 - essdist: initialized bestInliers; added check for negative threshold in the single solution case; removed abs

function [E, inliers] = ransacfitessmatrix(x1, x2, t, feedback)

  if ~all(size(x1)==size(x2))
    error('Data sets x1 and x2 must have the same dimension');
  end
    
  if nargin == 3
    feedback = 0;
  end
    
  [rows,npts] = size(x1);
  if ~(rows==2 || rows==3)
    error('x1 and x2 must have 2 or 3 rows');
  end
    
  if rows == 2    % Pad data with homogeneous scale factor of 1
    x1 = [x1; ones(1,npts)];
    x2 = [x2; ones(1,npts)];        
  end
    
  % Normalise each set of points so that the origin is at centroid and
  % mean distance from origin is sqrt(2).  normalise2dpts also ensures the
  % scale parameter is 1.  Note that 'essmatrix' will also call
  % 'normalise2dpts' but the code in 'ransac' that calls the distance
  % function will not - so it is best that we normalise beforehand.
  %%%[x1, T1] = normalise2dpts(x1);
  %%%[x2, T2] = normalise2dpts(x2);

  s = 5;  % Number of points needed to fit an essential matrix.
    
  fittingfn = @essmatrix;
  distfn    = @essdist;
  degenfn   = @isdegenerate;
  % x1 and x2 are 'stacked' to create a 6xN array for ransac
  [E, inliers] = ransac([x1; x2], fittingfn, distfn, degenfn, s, t, feedback);

  if isempty(E)  % ransac failed to find a solution.
    return;     % Do not attempt to do a final non-linear least squares fit
  end
    
  % Now do a final non-linear least squares fit on the data points considered to
  % be inliers.
  E = essmatrixNLS(E, x1(:,inliers), x2(:,inliers));
    
  % Denormalise
  %%%E = T2'*E*T1;
    
%--------------------------------------------------------------------------
% Function to evaluate the first order approximation of the geometric error
% (Sampson distance) of the fit of an essential matrix with respect to a
% set of matched points as needed by RANSAC.  See: Hartley and Zisserman,
% 'Multiple View Geometry in Computer Vision', 2nd Ed. page 287.
%
% Note that this code allows for E being a cell array of essential matrices of
% which we have to pick the best one. (A 5 point solution can return up to 10
% solutions, less on average)

function [bestInliers, bestE] = essdist(E, x, t);
    
  x1 = x(1:3,:);    % Extract x1 and x2 from x
  x2 = x(4:6,:);
    
    
  if iscell(E)  % We have several solutions each of which must be tested
		  
    nE = length(E);   % Number of solutions to test
    bestE = E{1};     % Initial allocation of best solution
    bestInliers = []; % Initial inliers
    ninliers = 0;     % Number of inliers
	
    for k = 1:nE
      Ex1 = E{k}*x1;
      Etx2 = E{k}'*x2;     

      %x2tEx1 = zeros(1,length(x1));
      %for n = 1:length(x1)
      %  x2tEx1(n) = x2(:,n)'*E{k}*x1(:,n);
      %end
      % faster calculation of x2tEx1
      x2tEx1 = sum(x2.*Ex1);
	    
      % Evaluate distances
      d =  x2tEx1.^2 ./ ...
          (Ex1(1,:).^2 + Ex1(2,:).^2 + Etx2(1,:).^2 + Etx2(2,:).^2);

      if t<0
        sorted = sort(d(:));
        if t>-1 % interpret t as outlier fraction
          pcent = -t;
        else
          pcent = 0.3;
        end
        t = sorted(uint32((1-pcent)*length(sorted))); % ignore distances > (100-pcent*100)% of all distances 
        fprintf('ransacfitessmatrix: threshold at %g%% outliers is %g\n', pcent*100, t);
        % alternative: t = prctile(d(:), (1-pcent)*100);
      end
	    
      % d>0, hence abs(d)==d
      inliers = find(d < t);     % Indices of inlying points
	    
      if length(inliers) > ninliers   % Record best solution
        ninliers = length(inliers);
        bestE = E{k};
        bestInliers = inliers;
      end
  end
    
  else     % We just have one solution
    Ex1 = E*x1;
    Etx2 = E'*x2;     
	
    %x2tEx1 = zeros(1,length(x1));
    %for n = 1:length(x1)
    %  x2tEx1(n) = x2(:,n)'*E*x1(:,n);
    %end
    % faster calculation of x2tEx1
    x2tEx1 = sum(x2.*Ex1);
	
    % Evaluate distances
    d =  x2tEx1.^2 ./ ...
        (Ex1(1,:).^2 + Ex1(2,:).^2 + Etx2(1,:).^2 + Etx2(2,:).^2);
	
    if t<0
      sorted = sort(d(:));
      if t>-1 % interpret t as outlier fraction
        pcent = -t;
      else
        pcent = 0.3;
      end
      t = sorted(uint32((1-pcent)*length(sorted))); % ignore distances > (100-pcent*100)% of all distances 
      fprintf('ransacfitessmatrix: threshold at %g%% outliers is %g\n', pcent*100, t);
      % alternative: t = prctile(d(:), (1-pcent)*100);
    end
	
    bestInliers = find(d < t);     % Indices of inlying points
    bestE = E;                          % Copy E directly to bestE
	
  end
	


%----------------------------------------------------------------------
% (Degenerate!) function to determine if a set of matched points will result
% in a degeneracy in the calculation of an essential matrix as needed by
% RANSAC.  This function assumes this cannot happen...
     
function r = isdegenerate(x)
  r = 0;    
 
