classdef (InferiorClasses={?ss,?tf,?zpk,?frd,?pmat,?ureal,?ucomplex,...
         ?ucomplexm,?ultidyn,?udyn,?umat,?uss,?ufrd}) pass < p_ss
    
    % PASS creates a PASS object decribing an affine LPV model as
    %
    %   pdG = sys(:,:,1) + sum_i p_i*sys(:,:,i)
    %
    % Use:
    %   pdG = PASS(): returns an empty object
    %   pdG = PASS(A,B,C,D,set,'prop',value)
    %            A,B,C,D: 3D matrices corresponding to each terms
    %            set: parameter set description, a pset object or
    %                 (np x 2) numeric matrix
    %            prop: properties: InputName, StateName, OutputName, 
    %                 (see ss for more details)      
    %   pdG = PASS(sys,set)
    %            sys: 3D ss object with the elementary models
    %   pdG = PASS(psys) conversion from a PSYS object (lmitool)
    %    
    % PASS Properties:
    %   a, b, c, d and others from ss object
    %   parset     - pset object with the parameter set information
    %
    % PASS Methods:
    %   pass       - class constructor
    %   iosize     - returns the number of input and outputs
    %   order      - returns the number of states
    %   npar       - returns the number of parameters
    %   nsys       - returns the number of LTI models used to describe 
    %                the LPV model
    %   size       - returns the model dimensions
    %   ispd       - checks if the model is parameter dependent
    %   isempty    - checks if a empty object
    %   subs       - evaluates an LPV model at a frozen parameter values
    %   ss         - converts to ss-object (tf & zpk also available)
    %   uss        - converts to uss-object
    %   psys       - converts to psys object (lmitool)
    %
    % other overloaded functions: 
    %   eig, pole, tzero, dcgain
    %   series, parallel, connect, etc.
    %   bode, bodemag, step, etc.
    %
    % See also ss, ppss
    
    % fbianchi - 2021-03-30
    
    methods

        function obj = pass(varargin)

            % PASS creates a PASS object decribing an affine LPV model as
            %
            %   pdG = sys(:,:,1) + sum_i p_i*sys(:,:,i)
            %
            % Use:
            %   pdG = PASS(): returns an empty object
            %   pdG = PASS(A,B,C,D,set,'prop',value)
            %            A,B,C,D: 3D matrices corresponding to each terms
            %            set: parameter set description, a pset object or
            %                 (np x 2) numeric matrix
            %            prop: properties: InputName, StateName, OutputName,
            %                 (see ss for more details)
            %   pdG = PASS(sys,set)
            %            sys: 3D ss object with the elementary models
            %   pdG = PASS(psys) conversion from a PSYS object (lmitool)
            
            % create a ss object (superclass)
            obj@p_ss();
            
            % ------------------------------------------------------------
            % empty object
            if (nargin == 0)
                return

            % ------------------------------------------------------------
            % to manage the copy of objects
            elseif isa(varargin{1},'pass')
                obj = varargin{1};
                
            % ------------------------------------------------------------
            % obj = PASS(psys) --> conversion from a PSYS object (lmitool)
            elseif ispsys(varargin{1})
                
                % psys object
                pdG = varargin{1};
                [typ,nv,ns,ni,no] = psinfo(pdG);
                if strcmp(typ,'pol')
                    error('PASS:PASS:inputError','PSYS must be affine, see ppss')
                end
                pv  = psinfo(pdG,'par');
                
                % system data
                as(ns,ns,nv) = 0;
                bs(ns,ni,nv) = 0;
                cs(no,ns,nv) = 0;
                ds(no,ni,nv) = 0;
                for ii = 1:nv
                    [as(:,:,ii),bs(:,:,ii),cs(:,:,ii),ds(:,:,ii),~] = ltiss(psinfo(pdG,'sys',ii));
                end
                obj.A = as;
                obj.B = bs;
                obj.C = cs;
                obj.D = ds;
                
                % pvset construction
                if isempty(pv)
                    error('PASS:PASS:inputError','PSYS must have a parameter set')
                    
                else
                    typpv = pvinfo(pv);
                    if strcmp(typpv,'box')
                        obj.parset = pset.Box('pvec',pv);
                    elseif strcmp(typpv,'pol')
                        obj.parset = pset.Gral('pvec',pv);
                    else
                        error('PASS:PASS:inputError','Invalid parameter set')
                    end
                end
            
            % ------------------------------------------------------------
            % obj = PASS(A,B,C,D,set,'prop',value) --> standard syntax
            elseif isnumeric(varargin{1}) && isnumeric(varargin{2}) && ...
                   isnumeric(varargin{3}) && isnumeric(varargin{4}) 
                % 
                if (nargin < 4)
                    error('PASS:PASS:notEnoughInputs',...
                        'Not enough input arguments')
                end
                
                % system matrices (to manage matrices w/different 3rd dim)
                nvs = cellfun('size',varargin(1:4),3);
                nv  = max(nvs);
                dims = size(varargin{1});
                extraA = repmat(zeros(dims(1:2)),1,1,nv-nvs(1));
                dims = size(varargin{2});
                extraB = repmat(zeros(dims(1:2)),1,1,nv-nvs(2));
                dims = size(varargin{3});
                extraC = repmat(zeros(dims(1:2)),1,1,nv-nvs(3));
                dims = size(varargin{4});
                extraD = repmat(zeros(dims(1:2)),1,1,nv-nvs(4));
                
                obj.A = cat(3,varargin{1},extraA);
                obj.B = cat(3,varargin{2},extraB);
                obj.C = cat(3,varargin{3},extraC);
                obj.D = cat(3,varargin{4},extraD);
                
                % rest of input arguments
                if (nargin > 5)
                    if ischar(varargin{6})
                        % property/value
                        na = length(varargin) - 5;
                        if (rem(na,2) ~= 0)
                            error('PASS:PASS:inputError',...
                                'A pair Property/Valus must be given')
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
                                        error('PASS:PASS:inputError',...
                                            'Invalid property/value input argument')
                                        
                                end
                            end
                        end
                    else
                        error('PASS:PASS:inputError','Invalid input argument')
                    end
                end
                
                % parameter set
                if (nargin < 5) && (nv == 1)
                    % case parameter independent
                    obj.parset = pset.Box;
                elseif (nargin > 4) && isa(varargin{5},'pset') && ...
                        (size(varargin{5}) == nv-1)
                    % from pset object
                    obj.parset = varargin{5};
                % - from np x 2 matrix
                elseif (nargin > 4) &&  (isnumeric(varargin{5}) && ...
                        size(varargin{5},2) == 2 && size(varargin{5},1) == nv-1) 
                    obj.parset = pset.Box(varargin{5});
                else
                    error('PASS:PASS:inputError',...
                        'The parameter set must be a pset object of %.0f parameters\n or a %.0fx2 matrix',...
                        nv-1,nv-1)
                end
                
            % ------------------------------------------------------------
            % obj = PASS(sys,set) 
            elseif isa(varargin{1},'ss')
                % 
                if (nargin < 2) 
                    if ismatrix(varargin{1})
                        varargin{2} = pset.Box; 
                    
                    else
                        error('PASS:PASS:inputError',...
                              'Insufficient input arguments')
                    end
                end
                % matrices info from ss 3d object
                sys    = varargin{1};
                % system matrices
                obj.A  = sys.a;
                obj.B  = sys.b;
                obj.C  = sys.c;
                obj.D  = sys.d;
                
                % parameter set
                % - from pset object
                if isa(varargin{2},'pset')
                    obj.parset = varargin{2};
                % - from np x 2 matrix
                elseif (isnumeric(varargin{2}) && size(varargin{2},2) == 2 ...
                            && size(varargin{2},1) == size(obj.A,3)-1) 
                    obj.parset = pset.Box(varargin{2});
                else
                    error('PASS:PASS:inputError',...
                        'The parameter set must be a pset.pset object or a %.0fx2 matrix',...
                        size(obj.A,3))
                end                
                
                % state, input and output names
                obj.StateName   = sys.StateName;
                obj.InputName   = sys.InputName;
                obj.OutputName  = sys.OutputName;
                obj.InputGroup  = sys.InputGroup;
                obj.OutputGroup = sys.OutputGroup;

            % ------------------------------------------------------------
            % obj = PASS(usys,pdG) --> conversion from a USS object
            elseif isa(varargin{1},'uss')
                usys   = varargin{1};
                
                % system matrices: reconstruction from the umat object
                parName = fields(usys.Uncertainty); 
                nv = length(parName);
                ns = order(usys);
                [no,ni] = iosize(usys);
                % parameter set:
                if (nargin > 1)
                    % taken the parameter from pdG and forcing the same
                    % parameter order
                    if ~isa(varargin{2},'pass');
                        error('PASS:PASS:inputError',...
                            'pdG must be a pass object');
                    else
                        pdG = varargin{2};
                    end
                    if (nv+1 ~= nsys(pdG))
                        error('PASS:PASS:inputError',...
                            'USS object must have the same parameter number than PASS object');
                    end
                    parNamePdG = pdG.parset.ParameterNames;
                    if ~isempty(setdiff(parName,parNamePdG));
                        error('PASS:PASS:inputError',...
                            'USS and PASS objects must have the same parameter names');
                    end
                    % takes parset from other pass object
                    obj.parset = pdG.parset;
                else
                    % built from uss info
                    parRange(nv,2) = 0;
                    parNamePdG = parName;
                    for ii = 1:nv
                        parRange(ii,:) = eval(sprintf('usys.Uncertainty.%s.range',parNamePdG{ii}));
                    end
                    obj.parset = pset.Box(parRange,[],parNamePdG);
                end
                
                as = zeros(ns,ns,nv+1);
                bs = zeros(ns,ni,nv+1);
                cs = zeros(no,ns,nv+1);
                ds = zeros(no,ni,nv+1);
                for ii = 1:nv+1
                    for jj = 1:nv
                        if (ii == 1 || jj ~= ii-1)
                            vpar.(parNamePdG{jj}) = 0;
                        else
                            vpar.(parNamePdG{jj}) = 1;
                        end
                    end
                    if (ii == 1)
                        [as(:,:,ii),bs(:,:,ii),cs(:,:,ii),ds(:,:,ii)] = ssdata(usubs(usys,vpar));
                    else
                        [atmp,btmp,ctmp,dtmp] = ssdata(usubs(usys,vpar));
                        as(:,:,ii) = atmp - as(:,:,1);
                        bs(:,:,ii) = btmp - bs(:,:,1);
                        cs(:,:,ii) = ctmp - cs(:,:,1);
                        ds(:,:,ii) = dtmp - ds(:,:,1);
                    end
                end
                obj.A = as;
                obj.B = bs;
                obj.C = cs;
                obj.D = ds;

                % state, input and output names
                obj.StateName   = usys.StateName;
                obj.InputName   = usys.InputName;
                obj.OutputName  = usys.OutputName;
                obj.InputGroup  = usys.InputGroup;
                obj.OutputGroup = usys.OutputGroup;
            
            end 
        end
        
     end
end
