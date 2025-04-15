function [alpha,idxPoints] = cvxdec(obj,point)

% CVXDEC computes a convex decomposition corresponding to a given point
%
% Use: 
%   alpha = CVXDEC(set,point)
%   [alpha,idx] = CVXDEC(set,point)
%
% Input:
%   set:    pset.Gral object
%   point:  parameter value
%
% Output:
%   alpha:  coefficients of the convex decomposition, such that  
%               sum(alpha(:)) = 1, all(alpha) >= 0
%   idx:    indices corresponding to the nearest points to point
%
% The convex decomposition is computed as  
%
% 	point = alpha(1)*Vx(idx(1)) + ... + alpha(nv)*Vx(idx(nv)), or 
%   point = Vx(:,idx)*alpha
%
% where Vx are the points defining the set.
%
% See also pset.Gral

% fbianchi - 2021-03-30

if (nargin < 2)
    error('PSET:GRAL:CVXDEC:inputError','Insufficient input arguments');
end

if isempty(obj)
    error('PSET:GRAL:CVXDEC:inputError','CVXDEC is not available for empty sets');
end

% checking the inputs
if (length(point) ~= size(obj.points,1))
    error('PSET:GRAL:CVXDC:inputError','POINT must be a vector of lenght %2.0f',...
        size(obj.points,1))
end

% # of parameters
[np,nv] = size(obj.points);

% check if point is in the set and find the nearest points in the
% definition
if (nv == 1)
    if any(point ~= obj.points)
        error('PSET:GRAL:CVXDC:inputError','POINT is not in the parameter set')
    else
        idxPoints = 1;
        alpha = 1;
    end
else
    
    % use triangulation
    simplex = tsearchn(obj.points',obj.simplices,point');
    if isnan(simplex)
        error('PSET:GRAL:CVXDC:inputError','POINT is not in the parameter set')
    end
    idxPoints = obj.simplices(simplex,:);
    
    % decomposition
    n = length(idxPoints);
    nearestPoints = obj.points(:,idxPoints);
    alpha = [nearestPoints;ones(1,n)]\[point;1];
end

if (nargout == 1)
   aux = zeros(nv,1);
   aux(idxPoints) = alpha;
   alpha = aux; 
end


