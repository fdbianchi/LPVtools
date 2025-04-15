function [alpha,idxPoints] = cvxdec(obj,point)

% CVXDEC computes a convex decomposition corresponding to a give point
%
% Use: 
%   alpha = CVXDEC(set,point)
%   [alpha,idx] = CVXDEC(set,point)
%
% Input:
%   set:    pset.Box object
%   point:  parameter value
%
% Output:
%   alpha:  coefficients of the convex decomposition, such that  
%               sum(alpha(:)) = 1, all(alpha) >= 0
%   idx:    vertex indices
%
% The convex decomposition is computed as  
%
% 	point = alpha(1)*Vx(idx(1)) + ... + alpha(nv)*Vx(idx(nv)), or 
%   point = Vx(:,idx)*alpha
%
% where Vx are the set vertices.
%
% See also PSET.BOX

% fbianchi - 2021-03-29

if (nargin < 2)
    error('PSET:BOX:CVXDEC:inputError','Insufficient input arguments');
end

if isempty(obj)
    error('PSET:BOX:CVXDEC:inputError','CVXDEC is not available for empty sets');
end

% check if point is consistent with the parameter set
[bool,msg] = checkval(obj,point);
if (bool <= 0)
    error('PSET:BOX:CVXDEC:inputError',msg);
end

prange = obj.range;
np = 1;

% initialization
if (diff(prange(1,:)) == 0)
    alpha = 1;
    idx = 1;

else
    alpha_i = (prange(1,2) - point(1))/(prange(1,2) - prange(1,1));
    alpha = [alpha_i (1 - alpha_i)]';
    idx = [1 2];

end
% rest of vertices
for ii = 2:size(prange,1)
    if (diff(prange(ii,:)) > 0)
        alpha_i = (prange(ii,2) - point(ii))/(prange(ii,2) - prange(ii,1));
        alpha = [alpha*alpha_i; alpha*(1 - alpha_i)];
        np = np + 1;
    end
end

% indices corresponding to the nearest points in the set to point 
% (for coherence with other pset objects)
if (np == 1)
    idxPoints = idx;
    
elseif (np > 1)
    idx = num2cell(pgrid(repmat(idx,np,1)),2);
    idxPoints = sub2ind(2*ones(1,np),idx{:});
    
end

nv = size(obj.points,2);
if (nargout == 1)
   aux = zeros(nv,1);
   aux(idxPoints) = alpha;
   alpha = aux; 
end