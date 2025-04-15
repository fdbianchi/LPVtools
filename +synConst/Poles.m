classdef Poles < synConst.synConst
    
    % SYNCONST.POLES specifies a constraint on the closed-loop pole location
    %
    % Use:
    %   const = synConst.Poles(property,values,...)
    %
    % where:
    %  - MinDecay:      positive scalar value
    %  - MinDamping:    positive scalar value between 0 and 1
    %  - MaxFreq:       positive scalar value
    %  - region:        LMI region
    
    % fbianchi - 2021-06-30
    
    properties
        region     = [];    % LMI region (matrices L and M)
        MinDecay   = 0;     % Minimum decay rate   
        MinDamping = 0;     % Minimum damping
        MaxFreq    = inf;   % Maximum frequency
        
    end
    
    methods
        
        function obj = Poles(varargin)
            
            % SYNCONST.POLES specifies a constraint on the closed-loop pole location
            %
            % Use:
            %   const = synConst.Poles(property,values,...)
            %
            % where:
            %  - MinDecay: positive scalar value
            %  - MinDamping: positive scalar value between 0 and 1
            %  - MaxFreq: positive scalar value
            %  - region: LMI region
                                   
            if (nargin == 0) || isempty(varargin{1})
                % empty object
                return;
                
            elseif isa(varargin{1},'Poles')
                % object copies
                obj = varargin{1};
                
            else 
                % general constructor
                %
                for ii = 1:2:length(varargin)
                    jj = fix(ii/2) + 1; 
                    switch varargin{ii}
                        case 'MinDecay'
                            % half plane
                            obj.MinDecay = varargin{ii+1};
                            if ~(isnumeric(obj.MinDecay) && isscalar(obj.MinDecay))...
                                    || (obj.MinDecay < 0) || isinf(obj.MinDecay)
                                error('SYNCONST_POLES:InputError',...
                                    'MinDecay must a positive scalar lower than inf')
                                
                            elseif (obj.MinDecay >= 0) && (obj.MinDecay < inf)
                                obj.region.L{jj} =  2*obj.MinDecay;
                                obj.region.M{jj} =  1;
                                
                            end
                            
                        case 'MinDamping'
                            % conic sector
                            obj.MinDamping = varargin{ii+1};
                            if ~(isnumeric(obj.MinDamping) && isscalar(obj.MinDamping))...
                                    || (obj.MinDamping > 1) || (obj.MinDamping < 0)
                                error('SYNCONST_POLES:InputError',...
                                    'MinDamping must a numeric scalar between 0 and 1')
                                
                            elseif (obj.MinDamping > 0) && (obj.MinDamping <= 1)
                                beta = atan(sqrt(1-obj.MinDamping^2)/obj.MinDamping);
                                theta = 2*beta;
                                obj.region.L{jj} = zeros(2);
                                obj.region.M{jj} = [sin(theta/2) -cos(theta/2);
                                                    cos(theta/2) sin(theta/2)];
                                                
                            end
                            
                        case 'MaxFreq'
                            % disk with center in s=0
                            obj.MaxFreq = varargin{ii+1};
                            if ~(isnumeric(obj.MaxFreq) && isscalar(obj.MaxFreq))...
                                    || (obj.MaxFreq < 0)
                                error('SYNCONST_POLES:InputError',...
                                    'MaxFreq must a positive scalar lower than inf')
                                
                            elseif (obj.MaxFreq >= 0) && (obj.MaxFreq < inf)
                                obj.region.L{jj} = -obj.MaxFreq*eye(2);
                                obj.region.M{jj} = [0 0;1 0];
                                
                            end
                            
                        
                        case 'Region'
                            obj.region = varargin{ii};
                            
                        otherwise
                            error('keyword not valid')
                    end
                end
            end
        end    
        
        function varargout = char(obj)
            
            % CHAR method for synConst.Gain class

            if isempty(obj)
                varargout{1} = sprintf('\tEmpty synConst.Poles constraint.\n');
                
            else
                
                str1 = 's.t.';
                str2 = sprintf('MinDecay > %4.3f, MinDamping > %4.3f, MaxFreq < %d',...
                    obj.MinDecay,obj.MinDamping,obj.MaxFreq);
                str3 = '';
                
                if (nargout < 2)
%                     varargout{1} = sprintf('\tPoles constraint: %s, with %s\n',str2,str3);
                    varargout{1} = sprintf('\tPoles constraint: %s\n',str2);
                else
                    varargout = {str1,str2,str3};
                end
            end
        end
        
    end
end