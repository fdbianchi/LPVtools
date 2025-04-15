function lyapStruct = createLyapFcn(varargin)

% CREATELYAPFCN creates a struct with the Lyapunov function information
%
% CREATELYAPFCN returns a struct with the information needed by lpvsyn and
% lpvanalysis to design or check controllers using parameter dependent
% Lyapunov functions
%
% For lpvanalysis use:
%       lyapStruct = createLyapFcn('cl',X,dX)
%
% For lpvsyn (state-feedback) use:
%       lyapStruct = createLyapFcn('sf',Y,dY)
%
% For lpvsyn (output-feedback) use:
%       lyapStruct = createLyapFcn(X,dX,Y,dY)
% (X is mandatory and the rest of the arguments are optionals)
%
% X, dX, Y and dX can be numeric scalar/matrices or anonymous functions
%
% CREATELYAPFCN will check consistency of X, dX, Y and dX (e.g. is X=0, 
% dX is set to 0)
%
% See also lpvsyn, lpvanalysis

% fbianchi - 2021-07-06


if ischar(varargin{1}) && strcmp(varargin{1},'cl')
    
    % lyapStruct = createLyapFcn('cl',X,dX)
    if checkFcn(varargin{2})
        lyapStruct.parfcn  = varargin{2};
    else
        error('CREATELYAPFCN:InputError',...
            'Invalid Lyapunov Function')
    end
    if (nargin > 2)
        if checkFcn(varargin{3})
            if isnumeric(lyapStruct.parfcn)
                lyapStruct.dparfcn = 0;
            else
                lyapStruct.dparfcn  = varargin{3};
            end
        else
            error('CREATELYAPFCN:InputError',...
                'Invalid Lyapunov Function')
        end
    else
        lyapStruct.dparfcn  = 0;
    end
    
elseif ischar(varargin{1}) && strcmp(varargin{1},'sf')
    
    % lyapStruct = createLyapFcn('sf',Y,dY)
    if checkFcn(varargin{2})
        lyapStruct.parfcnY  = varargin{2};
    else
        error('CREATELYAPFCN:InputError',...
            'Invalid Lyapunov Function')
    end
    if (nargin > 2)
        if checkFcn(varargin{3})
            if isnumeric(lyapStruct.parfcnY)
                lyapStruct.dparfcnY = 0;
            else
                lyapStruct.dparfcnY  = varargin{3};
            end
        else
            error('CREATELYAPFCN:InputError',...
                'Invalid Lyapunov Function')
        end
    else
        lyapStruct.dparfcnY  = 0;
    end
    
else
    if (nargin == 1)
        
        % lyapStruct = createLyapFcn(X)
        if checkFcn(varargin{1})
            lyapStruct.parfcnX  = varargin{1};
        else
            error('CREATELYAPFCN:InputError',...
                'Invalid Lyapunov Function')
        end
        lyapStruct.dparfcnX = 0;
        lyapStruct.parfcnY  = 0;
        lyapStruct.dparfcnY = 0;
        
    elseif (nargin == 2)
        % lyapStruct = createLyapFcn(X,dX)
        if checkFcn(varargin{1})
            lyapStruct.parfcnX  = varargin{1};
        else
            error('CREATELYAPFCN:InputError',...
                'Invalid Lyapunov Function')
        end
        if checkFcn(varargin{2})
            if isnumeric(lyapStruct.parfcnX)
                lyapStruct.dparfcnX = 0;
            else
                lyapStruct.dparfcnX  = varargin{2};
            end
        else
            error('CREATELYAPFCN:InputError',...
                'Invalid Lyapunov Function')
        end
        lyapStruct.parfcnY  = 0;
        lyapStruct.dparfcnY = 0;
        
    elseif (nargin == 3)
        % lyapStruct = createLyapFcn(X,dX,Y)
        if checkFcn(varargin{1})
            lyapStruct.parfcnX  = varargin{1};
        else
            error('CREATELYAPFCN:InputError',...
                'Invalid Lyapunov Function')
        end
        if checkFcn(varargin{2})
            if isnumeric(lyapStruct.parfcnX)
                lyapStruct.dparfcnX = 0;
            else
                lyapStruct.dparfcnX  = varargin{2};
            end
        else
            error('CREATELYAPFCN:InputError',...
                'Invalid Lyapunov Function')
        end
        if checkFcn(varargin{3})
            lyapStruct.parfcnY  = varargin{3};
        else
            error('CREATELYAPFCN:InputError',...
                'Invalid Lyapunov Function')
        end
        lyapStruct.dparfcnY = 0;
        
    elseif (nargin == 4)
        % lyapStruct = createLyapFcn(X,dX,Y,dY)
        if checkFcn(varargin{1})
            lyapStruct.parfcnX  = varargin{1};
        else
            error('CREATELYAPFCN:InputError',...
                'Invalid Lyapunov Function')
        end
        if checkFcn(varargin{2})
            if isnumeric(lyapStruct.parfcnX)
                lyapStruct.dparfcnX = 0;
            else
                lyapStruct.dparfcnX  = varargin{2};
            end
        else
            error('CREATELYAPFCN:InputError',...
                'Invalid Lyapunov Function')
        end
        if checkFcn(varargin{3})
            lyapStruct.parfcnY  = varargin{3};
        else
            error('CREATELYAPFCN:InputError',...
                'Invalid Lyapunov Function')
        end
        if checkFcn(varargin{4})
            if isnumeric(lyapStruct.parfcnY)
                lyapStruct.dparfcnY = 0;
            else
                lyapStruct.dparfcnY  = varargin{4};
            end
        else
            error('CREATELYAPFCN:InputError',...
                'Invalid Lyapunov Function')
        end
        
    end
end


% internal functions
function bool = checkFcn(fcn)

bool = isnumeric(fcn) || isa(fcn,'function_handle');


