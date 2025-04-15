classdef GainH2g < synConst.synConst
    
    % SYNCONST.GAINH2G specifies a constraint on H2 generalized performance
    %
    %   factor*||Wout*G(outputs,inputs)*Win||_2g <= bound
    %
    % Use:
    %   const = synConst.Gain(inputs, outputs, [prop], [values],...)
    %
    % where:
    %  - inputs  = vector or cell array with the inputs
    %  - outputs = vector or cell array with the outputs
    %  - prop/values:
    %       factor: numeric value (default 1)
    %       bound: numeric value (defualt inf)
    %       Win: lti or numeric value
    %       Wout: lti or numeric value
    %
    % This specification states a constraint on inf-norm of the output when 
    % the input is a signal with 2-norm equal to 1
    
    % fbianchi - 2021-06-30
    
    properties
        factor  = 1;     % constant weight multiplying the norm
        bound   = inf;   % upper bound for the norm
        inputs  = [];    % input definition (vector or cell array)
        outputs = [];    % output definition (vector or cell array)
        Win     = [];    % input weight function
        Wout    = [];    % output weight function
        
    end
    
    methods
        
        function obj = GainH2g(varargin)
            
            % SYNCONST.GAIN specifies a constraint on H2 generalized performance
            %
            %   factor*||Wout*G(outputs,inputs)*Win||_2 <= bound
            %
            % Use:
            %   const = synConst.Gain(outputs, inputs, [prop], [values],...)
            %
            % where:
            %  - inputs  = vector or cell array with the inputs
            %  - outputs = vector or cell array with the outputs
            %  - prop/values: 
            %       factor: numeric value (default 1)
            %       bound: numeric value (defualt inf)
            %       Win: lti or numeric value
            %       Wout: lti or numeric value
            %
            % This specification states a constraint on inf-norm of the 
            % output when the input is a signal with 2-norm equal to 1
                                  
            if (nargin == 0)
                % empty object
                return;
                
            elseif isa(varargin{1},'GainH2g')
                % for object copies
                obj = varargin{1};
                
            else 
                % general constructor
                %
                % i/o map
                if iscellstr(varargin{1}) || ischar(varargin{1}) || isvector(varargin{1})
                    obj.inputs  = varargin{1};
                else
                    error('SYNCONST_GAINH2g:InputError',...
                          'Input must be a numeric vector or cellstr')
                end
                if iscellstr(varargin{2}) || ischar(varargin{2}) || isvector(varargin{2})
                    obj.outputs  = varargin{2};
                else
                    error('SYNCONST_GAINH2g:InputError',...
                          'Output must be a numeric vector or cellstr')
                end
                
                % property/value
                for ii = 3:2:length(varargin)
                    prop  = varargin{ii};
                    value = varargin{ii+1};
                    
                    switch prop
                        case 'factor'
                            if ~(isnumeric(value) && isscalar(value))...
                                    || (value < 0) || (value > 1)
                                error('SYNCONST_GAINH2g:InputError',...
                                    'Factor must be a positive scalar between 0 and 1')
                                
                            elseif (value >= 0) && (value < 1)
                                obj.factor = value;
                                
                            end
                            
                        case 'bound'
                            if ~(isnumeric(value) && isscalar(value))...
                                    || (value < 0)
                                error('SYNCONST_GAINH2g:InputError',...
                                    'Factor must be a non negative scalar')
                                
                            elseif (value >= 0) && (value < inf)
                                obj.bound = value;
                                
                            end
                            
                        case 'Win'
                            if isa(value,'lti') || isnumeric(value)
                                if isscalar(value) || ...
                                        (length(value.u) == 1 && length(value.y) == 1)
                                    obj.Win = value*eye(length(obj.inputs));
                                else
                                    obj.Win  = value;
                                end
                            elseif iscell(value)
                                obj.Win  = append(value{:});
                                
                            else
                                error('Invalid Win')
                            end
                            if (length(obj.inputs) ~= size(obj.Win,2))
                                error('Win with invalid dimensions')
                            end
                            
                        case 'Wout'
                            if isa(value,'lti') || isnumeric(value)
                                if isscalar(value) || ...
                                        (length(value.u) == 1 && length(value.y) == 1)
                                    obj.Wout = value*eye(length(obj.outputs));
                                else
                                    obj.Wout  = value;
                                end
                            elseif iscell(value)
                                obj.Wout  = append(value{:});
                                
                            else
                                error('Invalid Wout')
                            end
                            if (length(obj.outputs) ~= size(obj.Wout,2))
                                error('Wout with invalid dimensions')
                            end
                            
                        otherwise
                            error('Invalidad parameter')
                    end
                end
            end
        end 
        
        function varargout = char(obj)
            
            % CHAR method for synConst.GainH2g class

            if isempty(obj)
                varargout{1} = sprintf('\tEmpty synConst.GainH2g constraint.\n');
                
            else
                
                % inputs
                if isnumeric(obj.inputs)
                    strAux = num2str(obj.inputs);
                    
                elseif iscell(obj.inputs)
                    strAux = strjoin(obj.inputs);
                    
                else
                    strAux = obj.inputs;
                end
                strIns = sprintf('[%s]',strAux);
                
                % outputs
                if isnumeric(obj.outputs)
                    strAux = num2str(obj.outputs);
                    
                elseif iscell(obj.outputs)
                    strAux = strjoin(obj.outputs);
                    
                else
                    strAux = obj.outputs;
                end
                strOuts = sprintf('[%s]',strAux);
                
                if isa(obj.bound,'sdpvar') || isinf(obj.bound)
                    str1 = 'min';
                    str2 = sprintf('%4.2f*||Tzw||_2g ',obj.factor);
                    str3 = sprintf('z=%s, w=%s',strOuts,strIns);
                else                
                    str1 = 's.t.';
                    str2 = sprintf('%4.2f*||Tzw||_2g < %4.3f',obj.factor,obj.bound);
                    str3 = sprintf('z=%s, w=%s',strOuts,strIns);
                end
                
                if (nargout < 2)
                    if strcmp(str1,'min')
                        varargout{1} = sprintf('\tConstraint: min %s, with %s\n',str2,str3);
                    else
                        varargout{1} = sprintf('\tConstraint: %s, with %s\n',str2,str3);
                    end
                else
                    varargout = {str1,str2,str3};
                end
            end
        end

    end
end

