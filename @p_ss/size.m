function varargout = size(obj,dim)

% SIZE returns the dimensions of the LPV model pdG
%
% Use:
%   [ny,nu,ns,np,nv] = SIZE(pdG)
%   ny = SIZE(pdG,1)
%   nu = SIZE(pdG,2)
%   ns = SIZE(pdG,3)
%   np = SIZE(pdG,4)
%   nv = SIZE(pdG,5)
%
% where
%   - ny: number of outputs
%   - nu: number of inputs
%   - ns: number of states
%   - np: number of parameters
%   - nv: number of systems
%
% See also pass, ppss, pgss, pcss

% fbianchi - 2021-03-31


if (nargout == 0)
    disp(obj);
    
else
    
    if (nargin == 1)
        ns      = order(obj);
        [ny,nu] = iosize(obj);
        np      = npar(obj);
        nv      = nsys(obj);
        if (nargout == 1)
            varargout = {[ny,nu,ns,np,nv]};
        else
            varargout = {ny,nu,ns,np,nv};
        end
        
    elseif (nargin == 2) && isnumeric(dim) && isscalar(dim)
        if (nargout > 1)
            error('P_SS:SIZE:outputError','Wrong number of output arguments')
        end
        
        switch dim
            case 1
                varargout{1} = size(obj.D(:,:,1),1);
                
            case 2
                varargout{1} = size(obj.D(:,:,1),2);
                
            case 3
                varargout{1} = order(obj);
                
            case 4
                varargout{1} = npar(obj);
                
            case 5
                varargout{1} = nsys(obj);
            otherwise
                error('P_SS:SIZE:inputError',...
                    'DIMS must be a positive scalar lower than 5')
        end
    end
    
end
    