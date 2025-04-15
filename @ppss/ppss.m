classdef (InferiorClasses={?ss,?tf,?zpk,?frd,?pmat,?ureal,?ucomplex,...
        ?ucomplexm,?ultidyn,?udyn,?umat,?uss,?ufrd}) ppss < p_ss
    
    % PPSS creates a PPSS object describing a polytopic LPV model as
    %
    %   pdG = sum_i (alpha_i(p)*sys(:,:,i))
    %
    % where alphi_i >= 0 and all(alpha_i(p)) = 1
    %
    % Use:
    %   pdG = PPSS(A,B,C,D,set,'prop',value)
    %            A,B,C,D: 3D matrices corresponding to each terms
    %            set: parameter set description, can be:
    %                 - a pset object
    %                 - (np x nv) matrix (np # of parameters and
    %                      nv # of points (size(sys,3))
    %            prop: properties: InputName, StateName, OutputName, 
    %                 (see ss for more details)      
    %   pdG = PPSS(sys,set): 
    %            sys: 3D ss object with the elementary models
    %
    %   pdG = PPSS(): returns an empty object
    %   pdG = PPSS(pdGaff): conversion from a pass obj
    %   pdG = PPSS(psys): conversion from a psys obj (lmitool)
    %
    % The parameter set must have as many points as elements in the third
    % dimension of sys or system matrices A, B, C, D
    %
    % PPSS Properties:
    %   a, b, c, d and other from ss object
    %   parset    - pvset object with the parameter set information
    %
    % PPSS Methods:
    %   ppss       - class constructor
    %   iosize     - returns the number of input and outputs
    %   order      - returns the number of states
    %   npar       - returns the number of parameters
    %   nsys       - returns the number of LTI models used to describe 
    %                the LPV model
    %   size       - returns the model dimensions
    %   ispd       - checks if the model is parameter dependent
    %   isempty    - checks if a empty object
    %   subs       - evaluates an LPV model at a frozen parameter values
    %   ss         - converts into ss-object (tf & zpk also available)
    %   psys       - converts into psys object (lmitool)
    %
    % other overloaded functions: 
    %   eig, pole, tzero, dcgain
    %   series, parallel, connect, etc.
    %   bode, bodemag, setp, etc.
    %
    % See also pass, ss
    
    % fbianchi - 2021-03-31
    
   
    methods
        
        function obj = ppss(varargin)
 
            % PPSS creates a PPSS object describing a polytopic LPV model as
            %
            %   pdG = sum_i (alpha_i(p)*sys(:,:,i))
            %
            % where alphi_i >= 0 and all(alpha_i(p)) = 1
            %
            % Use:
            %   pdG = PPSS(A,B,C,D,set,'prop',value)
            %            A,B,C,D: 3D matrices corresponding to each terms
            %            set: parameter set description, can be:
            %                 - a pset object
            %                 - (np x nv) matrix (np # of parameters and
            %                      nv # of points (size(sys,3))
            %            prop: properties: InputName, StateName, OutputName,
            %                 (see ss for more details)
            %   pdG = PPSS(sys,set):
            %            sys: 3D ss object with the elementary models
            %
            %   pdG = PPSS(): returns an empty object
            %   pdG = PPSS(pdGaff): conversion from a pass obj
            %   pdG = PPSS(psys): conversion from a psys obj (lmitool)
            %
            % The parameter set must have as many points as elements in the third
            % dimension of sys or system matrices A, B, C, D

            obj@p_ss();
            
            % ------------------------------------------------------------
            % obj = PPSS() --> empty object
            if (nargin == 0)
                return
                
            % ------------------------------------------------------------
            % to manage the copy of objects
            elseif isa(varargin{1},'ppss')
                if (nargin == 1)
                    obj = varargin{1};
                end
                
            % ------------------------------------------------------------
            % pdG = PPSS(psys) --> conversion from a psys obj (lmitool)
            elseif ispsys(varargin{1})
                
                % psys object
                pdG = varargin{1};
                [typ,nv,ns,ni,no] = psinfo(pdG);
                if strcmp(typ,'aff')
                    pdG = aff2pol(pdG);
                end
                if (nargin == 2)
                    pv = varargin{2};
                else
                    pv = psinfo(pdG,'par');
                end
                % system data
                as(ns,ns,nv) = 0;
                bs(ns,ni,nv) = 0;
                cs(no,ns,nv) = 0;
                ds(no,ni,nv) = 0;
                for ii=1:nv
                    [as(:,:,ii),bs(:,:,ii),cs(:,:,ii),ds(:,:,ii)] = ltiss(psinfo(pdG,'sys',ii));
                end
                obj.A = as;
                obj.B = bs;
                obj.C = cs;
                obj.D = ds;
                % pvset construction
                if isempty(pv)
                    warning('The parameter set is not defined')
                else
                    typ = pvinfo(pv);
                    if strcmp(typ,'box')
                        obj.parset = pset.Box('pvec',pv);
                    else
                        obj.parset = pset.Gral('pvec',pv);
                    end
                end
                
            % ------------------------------------------------------------
            % obj = PPSS(A,B,C,D,set,'prop',value) --> standard syntax
            elseif isnumeric(varargin{1}) && isnumeric(varargin{2}) && ...
                   isnumeric(varargin{3}) && isnumeric(varargin{4}) 
                % 
                if (nargin < 4)
                    % missing parset
                    error('PPSS:PPSS:notEnoughInputs',...
                        'Not enough input arguments')
                end
                
                % system matrices (to manage matrices w/different 3rd dim)
                nvs = cellfun('size',varargin(1:4),3);
                nv  = max(nvs);
                if (nvs(1) == 1)
                    extraA = repmat(varargin{1},1,1,nv-nvs(1));
                    obj.A = cat(3,varargin{1},extraA);
                elseif (nvs(1) < nv)
                    error('PPSS:PPSS:inputError',...
                        'The 3rd dimension of A must 1 or %1.0f',nv)
                else
                    obj.A = varargin{1};
                end
                if (nvs(2) == 1)
                    extraB = repmat(varargin{2},1,1,nv-nvs(2));
                    obj.B = cat(3,varargin{2},extraB);
                elseif (nvs(2) < nv)
                    error('PPSS:PPSS:inputError',...
                        'The 3rd dimension of B must 1 or %1.0f',nv)
                else
                    obj.B = varargin{2};
                end
                if (nvs(3) == 1)
                    extraC = repmat(varargin{3},1,1,nv-nvs(3));
                    obj.C = cat(3,varargin{3},extraC);
                elseif (nvs(3) < nv)
                    error('PPSS:PPSS:inputError',...
                        'The 3rd dimension of C must 1 or %1.0f',nv)
                else
                    obj.C = varargin{3};
                end
                if (nvs(4) == 1)
                    extraD = repmat(varargin{4},1,1,nv-nvs(4));
                    obj.D = cat(3,varargin{4},extraD);
                elseif (nvs(4) < nv)
                    error('PPSS:PPSS:inputError',...
                        'The 3rd dimension of D must 1 or %1.0f',nv)
                else
                    obj.D = varargin{4};
                end
                
                % rest of input arguments
                if (nargin > 5)
                    % property/value
                    na = length(varargin) - 5;
                    if (rem(na,2) ~= 0)
                        error('PPSS:PPSS:inputError',...
                            'A pair Property/Valus must be provided')
                    else
                        for jj = 1:2:na
                            prop  = varargin{jj+5};
                            value = varargin{jj+6};
                            switch prop
                                case 'StateName'
                                    obj.StateName = value;
                                    
                                case 'StateUnit'
                                    obj.StateUnit = value;
                                    
                                case 'InputName'
                                    obj.InputName = value;
                                    
                                case 'InputUnit'
                                    obj.InputUnit = value;
                                    
                                case 'OutputName'
                                    obj.OutputName = value;
                                    
                                case 'OutputUnit'
                                    obj.OutputUnit = value;
                                    
                                case 'InputGroup'
                                    obj.InputGroup = value;
                                    
                                case 'OutputGroup'
                                    obj.OutputGroup = value;
                                    
                                otherwise
                                    error('PPSS:PPSS:inputError',...
                                        'Invalid property/value input argument')
                                    
                            end
                        end
                    end
                end
                
                % parameter set
                if (nargin < 5) && (nv == 1)
                    % case parameter independent
                    obj.parset = pset.Gral;
                    
                elseif (nargin > 4) && isa(varargin{5},'pset') && ...
                        (size(varargin{5},2) == nv)
                    % from pset object
                    obj.parset = varargin{5};
                    
                elseif (nargin > 4) &&  isnumeric(varargin{5}) && ...
                        (size(varargin{5},2) == nv) 
                    % from np x nv matrix
                    obj.parset = pset.Gral(varargin{5});
                    
                else
                    error('PPSS:PPSS:inputError',...
                        'The parameter set must be a pset object with %.0f points\n or a npx%.0f matrix',...
                        nv,nv)
                end
                
                
            % ------------------------------------------------------------
            % pdG = PPSS(pdGaff) --> conversion from a pass obj
            % [must be before the general constructor because pass is also ss]
            elseif isa(varargin{1},'pass')
                sys = ss(varargin{1});
                obj = ppss(sys,varargin{1}.parset);
           
            % ------------------------------------------------------------
            % obj = PPSS(sys,set)
            elseif isa(varargin{1},'ss')
                
                if (nargin < 2) 
                    if ismatrix(varargin{1})
                        varargin{2} = pset.Gral(1); 
                    
                    else
                        error('PPSS:PPSS:inputError',...
                              'Insufficient input arguments')
                    end
                end
                
                % system matrices
                obj.A  = varargin{1}.a;
                obj.B  = varargin{1}.b;
                obj.C  = varargin{1}.c;
                obj.D  = varargin{1}.d;
                
                % parameter set
                pv = varargin{2};
                % - from pset object
                if isa(pv,'pset') && (size(pv.points,2) == size(obj.A,3))
                    obj.parset = pv;
                % - from npxnv matrix (vertices)  
                elseif (isnumeric(pv) && size(pv,2) == size(obj.A,3))
                    obj.parset = pset.Gral(pv);
                    
                else
                    error('PPSS:PPSS:inputError',...
                          'SET must be a pset object or a npx%.0f matrix',...
                          size(obj.A,3))
                end
                
                % state, input and output names
                warning('off','Control:ltiobject:RepeatedChannelNames')
                obj.StateName   = varargin{1}.StateName;
                obj.InputName   = varargin{1}.InputName;
                obj.OutputName  = varargin{1}.OutputName;
                obj.InputGroup  = varargin{1}.InputGroup;
                obj.OutputGroup = varargin{1}.OutputGroup;
                warning('on','Control:ltiobject:RepeatedChannelNames')
                
            end
        end
        
     end
end
