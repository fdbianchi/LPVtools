classdef pset
    
    % PSET is a superclass to specify a parameter set
    %
    % See also PSET.Box, PSET.Grid, PSET.Hull, PSET.Gral
    
    % fbianchi - 2021-03-25

    properties (Dependent, SetAccess = private)
        range  = [];      % (np x 2) matrix with the minimum and maximum 
                          % values of the parameters
    end
    
    properties (SetAccess = protected)
        points = [];      % (np x nv) matrix with the polytope vertices 
                          % or grid points
    end
    
    properties (SetAccess = public)
        rate   = [];      % (np x 2) matrix with the minimum and maximum 
                          % rate values of the parameters
    end
    
    properties
        ParameterNames = [];   % parameter names
    end
    
    methods
        
        function obj = pset(varargin)
            
            % PSET is a superclass to specify a parameter set
            %
            % See also PSET.Box, PSET.Grid, PSET.Hull, PSET.Gral
            
            % ------------------------------------------------------------
            % empty object
            if (nargin == 0)
                return;
                
            % ------------------------------------------------------------
            % to manage copies of objects
            elseif isa(varargin{1},'pset')
                obj = varargin{1};
                
            else
                error('PSET:inputError','Input argument not valid')
                
            end
            
        end
        
        % ----------------------------------------------------------------
        % GET functions
        %
        % GET <- range (implementation is necessary for dependent
        % properties)
        function rng = get.range(obj)
            rng = [min(obj.points,[],2), max(obj.points,[],2)];
            
        end
        %
        % GET <- ParameterNames
        function names = get.ParameterNames(obj)
            names = obj.ParameterNames;
            
        end
       
        % ----------------------------------------------------------------
        % SET
        %
        % SET -> rate
        function obj = set.rate(obj,newrate)
            np = size(obj);
            if (~isnumeric(newrate))
                error('PSET:inputError','NEWRATE must be numeric array')
            end
            [r,c] = size(newrate);
            if (np == 0)
                obj.rate = newrate;
            elseif ((r == np) || (r == 0)) && ((c == 2) || (c == 0))
                obj.rate = newrate;
            elseif ((r == 1) || (r == 0)) && ((c == 2) || (c == 0))
                obj.rate = repmat(newrate,np,1);
            else
                error('PSET:inputError',...
                    'NEWRATE must be a matrix with two columns and less than %g rows',np)
            end
        end
        
        % SET -> parameter names
        function obj = set.ParameterNames(obj,newnames)
            np = size(obj);
            if (np > 0 && iscellstr(newnames) && length(newnames) ~= np)
                error('PSET:inputError','NEWNAMES must be a cell of string with lenght %g',np)
            else
                obj.ParameterNames = newnames;
            end
        end
        
    end
    
end


