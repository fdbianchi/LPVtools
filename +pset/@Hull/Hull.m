classdef Hull < pset
    
    % PSET.HULL creates a convex hull parameter set
    %
    % Use:
    %   set = PSET.HULL(points,[rate],[names])
    %
    % Inputs:
    %   points:  (np x nv) matrix with point in the set,
    %                     points = [p1, p2,...,pn]
    %   rate:    (np x 2) matrix with minimun and maximun rate values,
    %                     rate(j,1)  <= d(pj)/dt  <= rate(j,2)
    %            (1 X 2)  vector defines the same minimun and
    %                     maximun rate values for all parameters
    %            (optional) 
    %   names:   string cell with the parameter names (optional)
    %
    % The function automatically finds the convex hull and discards
    % the interior points
    %
    % PSET.HULL Properties:
    %   range       - (np x 2) matrix with the minimum and maximum values of
    %                 the parameters
    %   rate        - (np x 2) matrix with the minimum and maximum rate values
    %                 of the parameters
    %   points      - (np x nv) matrix with the polytope vertices or grid
    %                 points
    %   ParameterNames - cell array with the parameter names
    %
    % PSET.HULL Methods:
    %   pset.Hull   - class constructor
    %   isempty     - returns a nonzero value if the set is empty
    %   size        - returns dimensions of the set
    %   checkval    - checks if a given parameter value is in the set
    %   cvxdec      - computes a convex decomposition
    %   plot        - parameter set illustration (1, 2 or 3 parameters)
    %
    % See also pass, ppss
    
    % fbianchi - 2021-03-26

    properties (SetAccess = private)
        simplices = []; % store the simplices for speeding up computations
        
    end
    
    methods
        
        function obj = Hull(varargin)
            
            % PSET.HULL creates a convex hull parameter set
            %
            % Use:
            %   set = PSET.HULL(points,[rate],[names])
            %
            % Inputs:
            %   points:  (np x nv) matrix with point in the set,
            %                     points = [p1, p2,...,pn]
            %   rate:    (np x 2) matrix with minimun and maximun rate values,
            %                     rate(j,1)  <= d(pj)/dt  <= rate(j,2)
            %            (1 X 2)  vector defines the same minimun and
            %                     maximun rate values for all parameters
            %   names:   string cell with the parameter names
            %
            % The function automatically finds the convex hull and discards
            % the interior points
            
            % ------------------------------------------------------------
            % empty object
            if (nargin == 0)
                return;
                
            % ------------------------------------------------------------
            % to manage copies of objects
            elseif isa(varargin{1},'Hull')
                obj = varargin{1};
                
            % ------------------------------------------------------------
            % case: set = pset(points,[rate],[names])
            elseif isnumeric(varargin{1})

                auxPoints = unique(varargin{1}','stable','rows')';
                [np,nv] = size(auxPoints);
                
                
                % get the convex hull
                if (nv == 1)
                    obj.points = auxPoints;
                    obj.simplices = 1;
                    
                else
                    if (np == 1)
                        % 1 parameter
                        if all(~(auxPoints - auxPoints(:,1)))
                            % when all points are the same
                            obj.points = auxPoints(1,1);
                            obj.simplices = 1;
                        else
                            [~,idx(1,1)] = min(auxPoints);
                            [~,idx(2,1)] = max(auxPoints);
                            obj.points = auxPoints(idx);
                            obj.simplices = [1 2];
                        end
                        
                    else
                        
                        % removing dimensions that don't changes
                        idx = all(~bsxfun(@minus, auxPoints, auxPoints(:,1)),2);
                        auxPoints_red = auxPoints(~idx,:);
                        
                        if (np == 2) || (np == 3)
                            % for 2d an 3d better to use these functions
                            if (nv <= np)
                                obj.simplices = 1:nv;
                                obj.points = auxPoints;
                            else
                                idx = convhull(auxPoints_red');
                                obj.points = auxPoints(:,unique(idx(:),'stable'));
                                obj.simplices = delaunay(obj.points');
                            end
                            
                        else
                            % for more than 3 parameters
                            idx = convhulln(auxPoints_red');
                            obj.points = auxPoints(:,unique(idx(:),'stable'));
                            obj.simplices = delaunayn(obj.points');
                        end
                    end
                end
                
                % parameter rate (optional)
                if (nargin > 1) && ~isempty(varargin{2})
                    if (size(varargin{2},1) == np || size(varargin{2},1) == 1)
                        obj.rate = varargin{2};
                    else
                        error('PSET:HULL:inputError','RATE must be a numeric matrix of dimension %g x 2',np)
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
                error('PSET:HULL:inputError','Input argument not valid')
                
            end
        end
        
    end

end
    

