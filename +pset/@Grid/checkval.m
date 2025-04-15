function [bool,msg] = checkval(obj,point)

% CHECKVAL checks if a parameter value point is in the set 
%
% Use:
%   [bool,msg] = CHECKVAL(set,point)
%
% Inputs:
%   set:    pset.Grid object
%   point:  parameter value
%
% Output:
%   bool:    0 = out of range, 1 = pval in the set, 0 < invalid value
%   msg:     error message
%
% See also PSET.GRID

% fbianchi - 2021-03-30

tol = 1e-3;

if (nargin < 2)
    bool = -3;
    msg  = 'insufficient input arguments';
    
else
    
    if isnumeric(point) && ismatrix(point)
        
        if isempty(point) && isempty(obj)
            % trivial case
            bool = true;
            msg = 'POINT is in the parameter set';

        elseif (size(point,1) == size(obj.points,1))
            
            % just to allow some tolerance in the comparisons
            prange = obj.range;
            ck_max = (max(point,[],2) - prange(:,2)) <=  tol*abs(prange(:,2));
            ck_min = (min(point,[],2) - prange(:,1)) >= -tol*abs(prange(:,1));
            
            if (all(ck_max) && all(ck_min))
                bool = true;
                msg = 'POINT is in the parameter set';
            else
                bool = false;
                msg = 'POINT is out of range';
            end
            
        else
            bool = -1;
            msg = sprintf('POINT must be a vector or a matrix with %g rows',size(obj.points,1));
        end
        
    else
        bool = -2;
        msg = 'POINT must be a numeric vector or matrix';
    end
    
end