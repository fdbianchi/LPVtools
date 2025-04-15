function newOpts = lpvsettings(varargin)

% LPVSETTINGS allows to set the optimization settings
%
% Use:
%   newOpts = LPVSETTINGS(Name, Value, ...)
%
% where opts is an struct with the some of the following fields
%
%   opts.eigtol  = 1e-6;             minimun value to decide positive definitiveness
%   opts.dkbnd   = 10000;            bound on Dk matrix
%   opts.penalty = 1e-6;             penalty on the norm(X) and norm(Y)
%   opts.XYpenalty = 1;              penalty on [X E; E Y], with E =
%                                    XYpenalty*I, this helps with
%                                    numerical issues
%   opts.Ctrlpenalty = 0;            penalty on controller matrices values
%
%   opts.subOpt  = false;            if true the controller is redesigned
%   opts.subOptBnd = 0.05;           for (1 + opts.subOptBnd)*gamma. This
%                                    also helps to get a controller with
%                                    better numerical characteristics
%   opts.verb    = true;             verbose mode
%   opts.debug   = false;            returns debug info
%
%   opts.solver  = 'sedumi';         solver: sedumi, sdpt3 or mosek
%   opts.solTol  = 1e-6;             solver tolerance
%   opts.maxiter = 100;              max solver iterations
%   opts.dualize = false;            check YALMIP help
%   opts.removeequalities = true;    check YALMIP help
%
% calling the function without arguments returns the default settings

% fbianchi - 2023-06-26


if (nargin == 0)
   % no arguments => return default settings
   
    newOpts.eigtol  = 1e-6;             % min eigenvalue to decide positive definitiveness
    newOpts.dkbnd   = 10000;            % bound on Dk matrix
    newOpts.penalty = 1e-6;             % penalty on the norm(X) and norm(Y)
    newOpts.XYpenalty = 1;              % penalty on [X I; I Y]
    newOpts.Ctrlpenalty = 0;            % penalty on controller matrices values
    
    newOpts.subOpt  = false;
    newOpts.subOptBnd = 0.05;
    newOpts.verb    = true;             % verbose mode
    newOpts.debug   = false;             % returns debug info

    newOpts.solver  = 'sedumi';         % solver: sedumi or sdpt3
    newOpts.solTol  = 1e-6;             
    newOpts.maxiter = 100;              % max solver iterations
    newOpts.dualize = false;            % check YALMIP help
    newOpts.removeequalities = true;    % check YALMIP help
    
else
    
    % new setting struct with default values
    newOpts = lpvsettings;
    
    for ii = 1:2:length(varargin)
        
        if isfield(newOpts, varargin{ii})
            
            % check for valid names
            checkField(varargin{ii}, varargin{ii+1});
            
            % changing values
            newOpts.(varargin{ii}) = varargin{ii+1};
            
        else
            error('LPVSYNSETTINGS:inputError',...
                '%s is not an acceptable setting', varargin{ii})
            
        end
        
    end
        
    
end


% ----------------------------------------------------------------------

function checkField(field, val)

    % checking for valid name-value pairs

    bool = false;
    switch field
        case {'eigtol', 'dkbnd', 'penalty', 'maxiter', 'solTol',...
                'subOptBnd', 'XYpenalty', 'Ctrlpenalty'}
            
            if ~isnumeric(val)
                msg = 'a numeric value';
                bool = true;
            end
            
        case {'subOpt','verb','debug'}
            if ~(isnumeric(val) || islogical(val))
                msg = 'logical value';
                bool = true;
            end
            
        case {'solver'}
            
            if ~ismember(val,{'sedumi', 'sdpt3', 'mosek'})
                msg = 'sedumi, sdpt3 or mosek';
                bool = true;
            end
            
        case {'removeequalities'}
            
            if ~ismember(val,{-1, 0, 1, 2})
                msg = '-1, 0, 1, or 2, check sdpsettings for help';
                bool = true;
            end
            
        case 'dualize'
            
            if ~(islogical(val) || (val == 1) || (val == 0))
                msg = 'must be true o false';
                bool = true;
            end
            
        otherwise
            error('Unexpected field')    
        
    end

    if bool
        error('LPVSYNSETTINGS:inputError','%s must be %s',...
            field, msg)
    end


