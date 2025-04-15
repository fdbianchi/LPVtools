classdef Box < pset
    
    % PSET.BOX creates a box parameter set
    %
    % Use:
    %   set = PSET.BOX(range,[rate],[names])
    %
    % Inputs:
    %   range:   (np x 2) matrix with minimun and maximun values,
    %                     range(j,1) <= pj  <= range(j,2)
    %   rate:    (np x 2) matrix with minimun and maximun rate values,
    %                     rate(j,1)  <= d(pj)/dt  <= rate(j,2)
    %            (1 X 2)  vector defines the same minimun and
    %                     maximun rate values for all parameters
    %                     [optional]
    %   names:   string cell with parameter names [optional]
    %
    % Alternativally, the set can be constructed from a pvec object (lmitool)
    % In this case, use
    %
    %   set = PSET.BOX('pvec',pvec_obj)
    % 
    % where:
    %   pvec_obj: parameter set object from lmitool
    %
    % PSET.BOX Properties:
    %   range       - (np x 2) matrix with the minimum and maximum values of 
    %                 the parameters
    %   rate        - (np x 2) matrix with the minimum and maximum rate values 
    %                 of the parameters
    %   points      - (np x nv) matrix with points in the description
    %   ParameterNames - cell array with the parameter names 
    %
    % PSET.BOX Methods:
    %   pset.Box    - class constructor
    %   checkval    - checks if a given parameter value is in the set
    %   cvxdec      - convex decomposition corresponding to a given point
    %   pvec        - converts a pvec object into a pset one
    %   plot        - parameter set illustration (1, 2 or 3 parameters)
    %   subsref     - returns a subset of the parameter set, use subscripts
    %
    % See also pass, ppss, pgss
            
    % fbianchi - 2021-03-25

    methods
        
        function obj = Box(varargin)
            
            % PSET.BOX creates a box parameter set
            %
            % Use:
            %   set = PSET.BOX(range,[rate],[names])
            %
            % Inputs:
            %   range:   (np x 2) matrix with minimun and maximun values,
            %                     range(j,1) <= pj  <= range(j,2)
            %   rate:    (np x 2) matrix with minimun and maximun rate values,
            %                     rate(j,1)  <= d(pj)/dt  <= rate(j,2)
            %            (1 X 2)  vector defines the same minimun and
            %                     maximun rate values for all parameters
            %                     [optional]
            %   names:   string cell with parameter names [optional]
            %
            % Alternativally, the set can be constructed from a pvec object from
            % LMItool. In this case, use
            %
            %   set = PSET.BOX('pvec',pvec_obj)
            %
            % where:
            %   pvec_obj: parameter set object from lmitool
            %
            % See also pass, ppss, pgss
            
            % ------------------------------------------------------------
            % empty object
            if (nargin == 0)
                return;
                
            % ------------------------------------------------------------
            % to manage copies of objects
            elseif isa(varargin{1},'Box')
                obj = varargin{1};
                
            % ------------------------------------------------------------
            % case: set = pset(range,[rate],[names])
            elseif isnumeric(varargin{1})
                
                % range
                if  (size(varargin{1},2) >= 2)
                    range  = sort(varargin{1},2);
                else
                    error('PSET:BOX:inputError','RANGE must be a numeric matrix of dimension np x 2')
                end
                
                % points correspond to the vertices of the box
                obj.points = pgrid(range,2);
                np = size(obj.points,1);

                % parameter rate (optional)
                if (nargin > 1) && ~isempty(varargin{2})
                    if (size(varargin{2},1) == np) || (size(varargin{2},1) == 1)
                        obj.rate  = varargin{2};
                    else
                        error('PSET:BOX:inputError','RATE must be a numeric matrix of dimension %g x 2',np)
                    end
                end
                
                % parameter names (optional)
                if (nargin > 2)
                    if iscellstr(varargin{3})
                        obj.ParameterNames = varargin{3};
                    elseif ischar(varargin{3})
                        obj.ParameterNames = varargin(3);
                    else
                        error('PSET:BOX:inputError','Names must be char or cell array')
                    end
                else
                    % default names
                    straux = sprintf('p%d,',1:np);
                    names = strsplit(straux,',');
                    obj.ParameterNames = names(1:end-1);
                end
                    
            % ------------------------------------------------------------
            % case: set = pset('pvec',pvec_obj)
            elseif (strcmp(varargin{1},'pvec') && isnumeric(varargin{2}))
                
                data = varargin{2};
                [typ,np] = pvinfo(data);
                if (strcmp(typ,'box'))
                    if all(all(data(1:np,[4 5])))
                        % only range
                        obj = pset.Box(data(:,[2 3]));
                    else
                        % with rate information
                        obj = pset.Box(data(1:np,[2 3]),data(1:np,[4 5]));
                    end
                else
                    error('PSET:BOX:inputError','The pvec set must be type box')
                end
                
            else
                error('PSET:BOX:inputError','Input argument not valid')
            end
        end
        
    end
    
end
    

