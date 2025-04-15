function [bool,msg] = checkval(obj,point)

% CHECKVAL checks if a parameter value point is in the set 
%
% Use:
%   [bool,msg] = CHECKVAL(set,point)
%
% Inputs:
%   set:    pset.Gral object
%   point:  parameter value
%
% Output:
%   bool:    0 = out of range, 1 = pval in the set, 0 < invalid value
%   msg:     string corresponding to the error message
%
% See also pset.Gral

% fbianchi - 2021-03-30

if (nargin < 2)
    bool = -3;
    msg = 'insufficient input arguments';
    
else
    
    if isnumeric(point) && ismatrix(point)
        
        if isempty(point) && isempty(obj)
            % trivial case
            bool = true;
            msg = 'point is in the parameter set';
            
        elseif (size(point,1) == size(obj.points,1))
            
            if (size(obj.points,2) == 1)
                if isequal(point,obj.points)
                    bool = true;
                    msg = 'POINT is consistent with the parameter set';
                else
                    bool = false;
                    msg = 'POINT is out of range';
                end
                
            else
                % find a simplex in which the point could be
                simplex = tsearchn(obj.points',obj.simplices,point');
                
                if isnan(simplex)
                    bool = false;
                    msg = 'POINT is out of range';
                else
                    bool = true;
                    msg = 'POINT is consistent with the parameter set';
                end
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
