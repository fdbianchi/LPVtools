function varargout  = size(obj,dim)

% SIZE returns the dimensions of a PSET object
%
% Use:
%   [np,nv] = size(set)
%   np = size(set)
%
% where
%   np: number of parameters
%   nv: number of points
%
% See also PSET.Box, PSET.Grid, PSET.Hull, PSET.Gral

% fbianchi - 2021-03-29

[np,nv] = size(obj.points);

if (nargout == 0)
    % if no output argument, then display set info
    disp(obj);
    
else
    
   if (nargin == 1)
        % if 1 input argument:
        if (nargout == 1)
            % then returns # of parameter
            varargout  = {np};
        else
            % or, then returns # of parameter and # of points
            varargout  = {np,nv};
        end
        
    elseif (nargin == 2)
        % if 2 input argument:
        if (nargout == 1)
            if isnumeric(dim) && isscalar(dim)
                switch dim
                    case 1
                        varargout  = {np};
                        
                    case 2
                        varargout  = {nv};
                        
                    otherwise
                        error('PSET:SIZE:inputError','DIMS must be a positive scalar lower than 2')
                end
            end
            
        else
            error('PSET:SIZE:outputError','Too many output arguments')
            
        end
    end
end