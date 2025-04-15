classdef (InferiorClasses={?ss,?tf,?zpk,?frd,?pmat,?ureal,?ucomplex,...
        ?ucomplexm,?ultidyn,?udyn,?umat,?uss,?ufrd}) p_ss < ss
    
    % P_SS is superclass to describe LPV models
    %
    % See also pass, ppss, pgss
    
    % fbianchi - 2021-03-31
    
    properties
        parset = [];    % parameter set
    end
    
    methods
        function obj = p_ss(varargin)

            % create a ss object (superclass)
            obj@ss();
            
            % ------------------------------------------------------------
            % empty object
            if (nargin == 0)
                return

            % ------------------------------------------------------------
            % to manage the copy of objects
            elseif isa(varargin{1},'p_ss')
                obj = varargin{1};
                
            end
        end
    end
 end
