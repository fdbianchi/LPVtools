function [bool,msg] = checkval(obj,point)

% CHECKVAL checks if a parameter value point is in the set 
%
% Use:
%   [bool,msg] = CHECKVAL(set,pval)
%
% Inputs:
%   set:    pset.Hull object
%   point:  parameter value
%
% Output:
%   bool:    0 = out of range, 1 = pval in the set, 0 < invalid value
%   msg:     string corresponding to the error message
%
% See also PSET.HULL

% fbianchi - 2021-03-30

if (nargin < 2)
    bool = -3;
    msg  = 'insufficient input arguments';
    
else
    
    if (isnumeric(point) && ismatrix(point))
        
        [np,nv] = size(obj.points);
        
        if isempty(point) && isempty(obj)
            % trivial case
            bool = true;
            msg = 'POINT is in the parameter set';
            
        elseif (size(point,1) == np)
            
            if (nv == 1)
                if any(point ~= obj.points)
                    simplex = NaN;
                else
                    simplex = 1;
                end
            else
                % find a simplex in which the point could be
                simplex = tsearchn(obj.points',obj.simplices,point');
            end
            
            if isnan(simplex)
                bool = false;
                msg = 'POINT is not in the parameter set';
            else
                bool = true;
                msg = 'POINT is in the parameter set';
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