% ESSMATRIX - computes essential matrix from 5 points
%
% Function computes the essential matrix from 5 matching points in
% a stereo pair of images. A minimal, 5 point algorithm is used.
%
% Usage:   Ecell = essmatrix(x1, x2)
%          Ecell = essmatrix(x)
%
% Arguments:
%          x1, x2 - Two 3xN sets of corresponding homogeneous points
%         
%          x      - If a single argument is supplied it is assumed that it
%                   is in the form x = [x1; x2]
% Returns:
%          Ecell  - A cell array of possible 3x3 essential matrices Ei such that x2'*Ei*x1 = 0.
%
% See Also: Kovesi's FUNDMATRIX

% Requires calibrated_fivepoint.m

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

function Ecell = essmatrix(varargin)
    
    [x1, x2, npts] = checkargs(varargin(:));
    %Octave = exist('OCTAVE_VERSION', 'builtin') == 5; % Are we running under Octave

    if(npts ~= 5)
      error('essmatrix: exactly 5 point pairs expected');
    end
    
    % Normalise each set of points so that the origin 
    % is at centroid and mean distance from origin is sqrt(2). 
    % normalise2dpts also ensures the scale parameter is 1.
    %%%[x1, T1] = normalise2dpts(x1);
    %%%[x2, T2] = normalise2dpts(x2);
    
    %%%%%%%%%%% Select a 5-point solver by uncommenting the corresponding code fragment %%%%%%%%%%% 

    % Use Stewenius' ISPRS Grobner basis solver
    Evec = calibrated_fivepoint(x1, x2);

    % Use Stewenius/Engels' ISPRS non-Grobner basis solver
    %{
    [SOLS, EE] = calibrated_fivepoint_non_gb(x1, x2);
    Evec = EE*[SOLS ones(size(SOLS,1),1) ]';
    Evec = Evec./(ones(9,1)*sqrt(sum(Evec.^2)));
    I = find(not(imag(Evec(1,:))));
    Evec = Evec(:,I);
    %}

    % Use Kukelova's BMVC polyeig solver
    %{
    Emat = peig5pt(x1, x2);
    Evec = reshape(Emat', 9, size(Emat,1)/3);
    Evec = Evec./(ones(9,1)*sqrt(sum(Evec.^2)));
    %}

    nsol = size(Evec,2);
    Emat = permute(reshape(Evec, 3, 3, nsol), [2,1,3]);
    Ecell = mat2cell(Emat, 3, 3, ones(1, nsol));

    % Denormalise
    %%%E = T2'*E*T1;
    
%--------------------------------------------------------------------------
% Function to check argument values and set defaults

function [x1, x2, npts] = checkargs(arg);
    
    if length(arg) == 2
        x1 = arg{1};
        x2 = arg{2};
        if ~all(size(x1)==size(x2))
            error('x1 and x2 must have the same size');
        elseif size(x1,1) ~= 3
            error('x1 and x2 must be 3xN');
        end
        
    elseif length(arg) == 1
        if size(arg{1},1) ~= 6
            error('Single argument x must be 6xN');
        else
            x1 = arg{1}(1:3,:);
            x2 = arg{1}(4:6,:);
        end
    else
        error('Wrong number of arguments supplied');
    end
      
    npts = size(x1,2);
    if npts < 5
        error('At least 5 points are needed to compute the essential matrix');
    end
    
