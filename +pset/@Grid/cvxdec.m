function [alpha,idxPoints] = cvxdec(obj,point)

% CVXDEC computes a convex decomposition corresponding to a given point
%
% Use: 
%   alpha = CVXDEC(set,point)
%   [alpha,idx] = CVXDEC(set,point)
%
% Input:
%   set:    pset.Grid object
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
% See also PSET.GRID

% fbianchi - 2021-03-31

if (nargin < 2)
    error('PSET:GRID:CVXDEC:inputError','Insufficient input arguments');
end

if isempty(obj)
    error('PSET:GRID:CVXDEC:inputError','CVXDEC is not available for empty sets');
end

% check if point is consistent with the parameter set
[bool,msg] = checkval(obj,point);
if (bool <= 0)
    error('PSET:GRID:CVXDEC:inputError',msg);
end

% computing the convex combination and the nearest points
np = size(obj.points,1);
gDims(np) = 0; idx = [];
for ii = 1:np
    
    % remove duplicated value for ii dimension
    aux = unique(obj.points(ii,:));
    gDims(ii) = length(aux);
    if (gDims(ii) == 1) 
        if (ii == 1)
            alpha = 1;
        end
        if (np == 1)
            idx = 1;
        end
        
    else
        [~,auxIdx] = sort(abs(aux-point(ii)));
        auxIdx = sort(auxIdx(1:2));
        id1 = auxIdx(1);
        id2 = auxIdx(2);
        
        alpha_i = (aux(id2) - point(ii))/(aux(id2) - aux(id1));
        
        % create the outputs
        if (ii == 1)
            alpha = [alpha_i (1-alpha_i)]';
            idx = [id1 id2];
        else
            alpha = [alpha*alpha_i; alpha*(1-alpha_i)];
            idx = [idx; id1 id2];
        end
    end
end

% indices corresponding to the nearest points in the set to point 
if (np == 1)
    idxPoints = idx;
elseif np > 1
    idx = num2cell(pgrid(idx),2);
    idxPoints = sub2ind(gDims,idx{:});
end

if (nargout == 1)
    nv = size(obj.points,2);
    aux = zeros(nv,1);
    aux(idxPoints) = alpha;
    alpha = aux;
end