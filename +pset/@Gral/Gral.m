classdef Gral < pset
    
    % PSET.GRAL creates a general parameter set defined by discrete points
    %
    % Use:
    %   set = PSET.GRAL(points,[rate],[names])
    %
    % Inputs:
    %   points:  (np x nv) matrix with the points in the set,
    %                     points = [p1, p2,...,pn]
    %   rate:    (np x 2) matrix with minimun and maximun rate values,
    %                     rate(j,1)  <= d(pj)/dt  <= rate(j,2)
    %            (1 X 2)  vector defines the same minimun and
    %                     maximun rate values for all parameters (optional)
    %   names:   string cell with parameter names (optional)
    %
    % Alternativally, the set can be constructed from a pvec object (lmitool)
    % In this case, use
    %
    %   set = PSET.GRAL('pvec',pvec_obj)
    % 
    % where:
    %   pvec_obj: parameter set lmitool object of type pol
    %
    % PSET.GRAL Properties:
    %   range       - (np x 2) matrix with the minimum and maximum values of
    %                 the parameters
    %   rate        - (np x 2) matrix with the minimum and maximum rate values
    %                 of the parameters
    %   points      - (np x nv) matrix with the polytope vertices or grid
    %                 points
    %   ParameterNames - cell array with the parameter names
    %
    % PSET.GRAL Methods:
    %   pset.Gral   - class constructor
    %   isempty     - returns a nonzero value if the set is empty
    %   size        - returns dimensions of the set
    %   checkval    - checks if a given parameter value is in the set
    %   cvxdec      - computes a convex decomposition
    %   pvec        - converts a pvec object into a pset one
    %   plot        - parameter set illustration (1, 2 or 3 parameters)
    %
    % See also pass, ppss, pgss
    
    % fbianchi - 2021-03-26

    properties (SetAccess = private)
        hullIndex = []; % point indeces corresponding to the convex hull
        simplices = []; % simplecis stored for speeding up computations
    end
    
    methods
        
        function obj = Gral(varargin)
            
            % PSET.GRAL creates a general parameter set defined by discrete points
            %
            % Use:
            %   set = PSET.GRAL(points,[rate],[names])
            %
            % Inputs:
            %   points:  (np x nv) matrix the with points in the set,
            %                     points = [p1, p2,...,pn]
            %   rate:    (np x 2) matrix with minimun and maximun rate values,
            %                     rate(j,1)  <= d(pj)/dt  <= rate(j,2)
            %            (1 X 2)  vector defines the same minimun and
            %                     maximun rate values for all parameters (optional)
            %   names:   string cell with parameter names (optional)
            %
            % Alternativally, the set can be constructed from a pvec object (lmitool)
            % In this case, use
            %
            %   set = PSET.GRAL('pvec',pvec_obj)
            %
            % where:
            %   pvec_obj: parameter set lmitool object of type pol
            
            % ------------------------------------------------------------
            % empty object
            if (nargin == 0)
                return;
                
            % ------------------------------------------------------------
            % to manage copies of objects
            elseif isa(varargin{1},'Gral')
                obj = varargin{1};
                
            % ------------------------------------------------------------
            % case: set = pset(points,[rate],[names])
            elseif isnumeric(varargin{1})

                % points in the defintion
                auxPoints = varargin{1};
                [np,nv] = size(auxPoints);
                
                % removing dimensions that don't changes
                if (nv > 1)
                    idx = all(~bsxfun(@minus, auxPoints, auxPoints(:,1)),2);
                    auxPoints(idx) = [];
                end
                
                % get the convex hull
                if (np == 1)
                    % 1 parameter
                    [~,idx(1,1)] = min(auxPoints);
                    [~,idx(2,1)] = max(auxPoints);
                    obj.hullIndex = idx;
                    if length(auxPoints) > 1
                        obj.simplices = delaunayn(auxPoints');
                    else
                        obj.simplices = 1;
                    end
                    
                elseif (np == 2) || (np == 3)
                    % better functions for 2 & 3 parameters
                    if (nv <= np)
                        obj.hullIndex = 1:nv;
                        obj.simplices = 1:nv;
                    else
                        idx = convhull(auxPoints');
                        obj.hullIndex = idx;%unique(idx(:));
                        obj.simplices = delaunay(auxPoints');
                    end
                else
                    if (nv <= np)
                        obj.hullIndex = 1:nv;
                        obj.simplices = 1:nv;
                    else
                        idx = convhulln(auxPoints');
                        obj.hullIndex = unique(idx(:));
                        obj.simplices = delaunayn(auxPoints');
                    end
                end
                obj.points = varargin{1};
                
                % parameter rate (optional)
                if (nargin > 1) && ~isempty(varargin{2})
                    if (size(varargin{2},1) == np || size(varargin{2},1) == 1)
                        obj.rate = varargin{2};
                    else
                        error('PSET:GRAL:inputError','RATE must be a numeric matrix of dimension %g x 2',np)
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
                
            % ------------------------------------------------------------
            % case: set = pset('pvec',pvec_obj)
            elseif (strcmp(varargin{1},'pvec') && isnumeric(varargin{2}))
                
                data = varargin{2};
                typ = pvinfo(data);
                
                if (strcmp(typ,'pol'))
                    % polytopic case
                    obj = pset.Gral(data(:,2:end));
                else
                    % box case
                    error('PSET:GRAL:inputError','PVEC is type box, use pset.Box')
                end
                
            else
                error('PSET:GRAL:inputError','Input argument not valid')
            end
        end
        
        
    end

    
end
    

