classdef Grid < pset
    
    % PSET.GRID creates a grid parameter set
    %
    % Use:
    %   set = PSET.GRID(range,npoints,[rate],[names])
    %   set = PSET.GRID(subsets,[rate],[names])
    %
    % Inputs:
    %   range:   (np x 2) matrix with minimun and maximun values,
    %                     range(j,1) <= pj  <= range(j,2)
    %
    %   npoints: scalar or (1 x np) vector with the number of
    %            points in each dimension,
    %   subsets: cell with np vectors of dimension 1 x n_j, the set
    %            is
    %               subset_1 x ... x subset_np
    %
    %   rate:    (np x 2) matrix with minimun and maximun rate values,
    %                     rate(j,1)  <= d(pj)/dt  <= rate(j,2)
    %            (1 X 2)  vector defines the same minimun and
    %                     maximun rate values for all parameters
    %            (optional)
    %   names:   string cell with parameter names (optional)
    %
    % PSET.GRID Properties:
    %   range       - (np x 2) matrix with the minimum and maximum values of 
    %                 the parameters
    %   rate        - (np x 2) matrix with the minimum and maximum rate values 
    %                 of the parameters
    %   points      - (np x nv) matrix with the points defining the set
    %   ParameterNames - cell array with the parameter names
    %
    % PSET.GRID Methods:
    %   PSET.GRID   - class constructor
    %   isempty     - returns a nonzero value if the set is empty
    %   size        - returns the dimensions of the set
    %   checkval    - checks if a given parameter value is in the set
    %   cvxdec      - convex decomposition corresponding to a given point
    %   plot        - parameter set illustration (1, 2 or 3 parameters)
    %   subsref     - returns a subset of the parameter set, use subscripts
    %
    % See also pass, ppss
            
    % fbianchi - 2021-03-26

    methods
        
        function obj = Grid(varargin)
            
            % PSET.GRID creates a grid parameter set
            %
            % Use:
            %   set = PSET.GRID(range,npoints,[rate],[names])
            %   set = PSET.GRID(subsets,[rate],[names])
            %
            % Inputs:
            %   range:   (np x 2) matrix with minimun and maximun values,
            %                     range(j,1) <= pj  <= range(j,2)
            %
            %   npoints: scalar or (1 x np) vector with the number of
            %            points in each dimension,
            %   subsets: cell with np vectors of dimension 1 x n_j, the set
            %            is
            %               subset_1 x ... x subset_np
            %
            %   rate:    (np x 2) matrix with minimun and maximun rate values,
            %                     rate(j,1)  <= d(pj)/dt  <= rate(j,2)
            %            (1 X 2)  vector defines the same minimun and
            %                     maximun rate values for all parameters
            %            (optional)
            %   names:   string cell with parameter names (optional)
            
            % ------------------------------------------------------------
            % empty object
            if (nargin == 0)
                return;
                
            % ------------------------------------------------------------
            % to manage copies of objects
            elseif isa(varargin{1},'Grid')
                obj = varargin{1};
                
            % ------------------------------------------------------------
            % case: set = pset(range,npoints,[rate],[names])
            elseif isnumeric(varargin{1})
                
                if (nargin < 1)
                    error('PSET:GRID:inputError','Insufficient input arguments')
                end
                    
                if (size(varargin{1},2) >= 2)
                    range  = sort(varargin{1},2);
                else
                    error('PSET:GRID:inputError','RANGE must be a numeric matrix of dimension np x 2')
                end
                
                % points of the grid
                if (nargin < 2)
                    obj.points = pgrid(range,2);
                    
                elseif ~isempty(varargin{2}) && isnumeric(varargin{2})
                    obj.points = pgrid(range,varargin{2});
                else
                    error('PSET:GRID:inputError','Invalid number of points')
                end
                np = size(obj.points,1);

                % parameter rate (optional)
                if (nargin > 2) && ~isempty(varargin{3})
                    if (size(varargin{3},1) == np) || (size(varargin{3},1) == 1)
                        obj.rate  = varargin{3};
                    else
                        error('PSET:GRID:inputError','RATE must be a numeric matrix of dimension %g x 2',np)
                    end
                end
                
                % parameter names (optional)
                if (nargin > 3)
                    if iscellstr(varargin{4})
                        obj.ParameterNames = varargin{4};
                    elseif ischar(varargin{4})
                        obj.ParameterNames = varargin(4);
                    else
                        error('PSET:GRID:inputError','Names must be char or cell array')
                    end
                else
                    % default names
                    straux = sprintf('p%d,',1:np);
                    names = strsplit(straux,',');
                    obj.ParameterNames = names(1:end-1);
                end
                     
             % ------------------------------------------------------------
             % case: set = pset(subSets,[rate],[names])
             elseif iscell(varargin{1})
                 
                 % generate the points from the subsets
                 subSets = varargin{1};
                 subSets = subSets(:);
                 np = length(subSets);
                 
                 [auxPoints{1:np}] = ndgrid(subSets{:});
                 for ii = np:-1:1
                     obj.points(ii,:) = auxPoints{ii}(:)';
                 end
                 
                 % parameter rate (optional)
                 if (nargin > 1) && ~isempty(varargin{2})
                     if (size(varargin{2},1) == np || size(varargin{2},1) == 1)
                         obj.rate = varargin{2};
                     else
                         error('PSET:GRID:inputError','RATE must be a numeric matrix of dimension %g x 2',np)
                     end
                 end
                 
                 % parameter names (optional)
                 if (nargin > 2 && iscellstr(varargin{3}))
                     obj.ParameterNames = varargin{3};
                 else
                     straux = sprintf('p%d,',1:np);
                     names = strsplit(straux,',');
                     obj.ParameterNames = names(1:end-1);
                 end
                 
            else
                error('PSET:GRID:inputError','Input argument not valid')
            end
        end
        
    end
    
end
    

