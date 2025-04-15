function [bool,msg] = checkval(obj,point)

% CHECKVAL checks if a parameter value point is in the set 
%
% Use:
%   [bool,msg] = CHECKVAL(set,point)
%
% Inputs:
%   set:    pset.Box object
%   point:  parameter value
%
% Output:
%   bool:    0 = out of range, 1 = pval in the set, 0 < invalid value
%   msg:     error message
%
% See also PSET.BOX

% fbianchi - 2021-03-29


if (nargin < 2)
    
    bool = -3;
    msg = 'insufficient input arguments';
    
else
    if isnumeric(point) && ismatrix(point)% && ~isempty(point))
        
        if isempty(point) && isempty(obj)
            % trivial case
            bool = true;
            msg = 'point is in the parameter set';
        
        elseif (size(point,1) == size(obj.points,1))
            
            % just to allow some tolerance in the comparisons
            tol = 1e-3;
            prange = obj.range;
            ck_max = (max(point,[],2) - prange(:,2)) <=  tol*abs(prange(:,2));
            ck_min = (min(point,[],2) - prange(:,1)) >= -tol*abs(prange(:,1));
            
            if (all(ck_max) && all(ck_min))
                bool = true;
                msg = 'point is in the parameter set';
            else
                bool = false;
                msg = 'point is not in the parameter set';
            end
            
        else
            bool = -1;
            msg = sprintf('point must be a vector or a matrix with %g rows',size(obj.points,1));
        end
        
    else
        bool = -2;
        msg = 'point must be a numeric vector or matrix';
    end
    
end