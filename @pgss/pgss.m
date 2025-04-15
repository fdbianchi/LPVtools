classdef (InferiorClasses={?ss,?tf,?zpk,?frd,?pmat,?ureal,?ucomplex,...
        ?ucomplexm,?ultidyn,?udyn,?umat,?uss,?ufrd}) pgss < p_ss

    % PGSS creates a PGSS object decribing an LPV model as
    %
    %   pdG = sys(:,:,1) + sum_i (fcn_i(p)*sys(:,:,i))
    %
    % Use:
    %   pdG = PGSS(): returns an empty object
    %   pdG = PGSS(A,B,C,D,set,fnc,'prop',value)
    %            A,B,C,D: 3D matrices corresponding to each terms
    %            set: parameter set description, a pset object or
    %                 (np x 2) numeric matrix
    %            fnc: anonymous function, 
    %            prop: properties: InputName, StateName, OutputName, 
    %                 (see ss for more details)      
    %   pdG = PGSS(sys,set,fnc)
    %            sys: 3D ss object with the elementary models
    %    
    % PASS Properties:
    %   a, b, c, d and others from ss object
    %   parset     - pset object with the parameter set information
    %   parfcn     - anonymous function
    %
    % PASS Methods:
    %   pgss       - class constructor
    %   iosize     - returns the number of input and outputs
    %   order      - returns the number of states
    %   npar       - returns the number of parameters
    %   nsys       - returns the number of LTI models used to describe 
    %                the LPV model
    %   size       - returns the model dimensions
    %   ispd       - checks if the model is parameter dependent
    %   isempty    - checks if a empty object
    %   subs       - evaluates an LPV model at a frozen parameter values
    %   ss         - converts into a ss-object (tf & zpk also available)
    %
    % other overloaded functions: 
    %   eig, pole, tzero, dcgain
    %   series, parallel, connect, etc.
    %   bode, bodemag, setp, etc.
    %
    % See also ss, pass, ppss
    
    % fbianchi - 2021-03-31
    
    properties (SetAccess = private)
        parfcn = NaN;
    end
    
    methods
        
        function obj = pgss(varargin)

            % PGSS creates a PGSS object decribing an LPV model as
            %
            %   pdG = sys(:,:,1) + sum_i (fcn_i(p)*sys(:,:,i))
            %
            % Use:
            %   pdG = PGSS(): returns an empty object
            %   pdG = PGSS(A,B,C,D,set,fnc,'prop',value)
            %            A,B,C,D: 3D matrices corresponding to each terms
            %            set: parameter set description, a pset object or
            %                 (np x 2) numeric matrix
            %            fnc: anonymous function,
            %            prop: properties: InputName, StateName, OutputName,
            %                 (see ss for more details)
            %   pdG = PGSS(sys,set,fnc)
            %            sys: 3D ss object with the elementary models
            
            % create a ss object (superclass)
            obj@p_ss();
            
            % ------------------------------------------------------------
            % empty object
            if (nargin == 0)
                return

            % ------------------------------------------------------------
            % to manage the copy of objects
            elseif isa(varargin{1},'pgss')
                if (nargin == 1)
                    obj = varargin{1};
                    
                end
                
            % ------------------------------------------------------------
            % obj = pgss(A,B,C,D,set,fnc,'prop',value) --> standard syntax
            elseif isnumeric(varargin{1}) && isnumeric(varargin{2}) && ...
                   isnumeric(varargin{3}) && isnumeric(varargin{4}) 
                % 
                if (nargin < 6)
                    % missing parset
                    error('PGSS:PGSS:notEnoughInputs',...
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
                
                % parameter functions fcn_i(p)
                if isa(varargin{6},'function_handle')
                    obj.parfcn = varargin{6};
                else
                    error('PGSS:PGSS:inputError','fnc must be anonymous function')
                end
                
                % rest of input arguments
                if (nargin > 7)
                    if ischar(varargin{7})
                        % property/value
                        na = length(varargin) - 6;
                        if (rem(na,2) ~= 0)
                            error('PGSS:PGSS:inputError',...
                                'A pair Property/Value must be provided')
                        else
                            for jj = 1:2:na
                                prop  = varargin{jj+6};
                                value = varargin{jj+7};
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
                                        error('PGSS:PGSS:inputError',...
                                            'Invalid property/value input argument')
                                        
                                end
                            end
                        end
                    else
                        error('PGSS:PGSS:inputError','Invalid input argument')
                    end
                end
                
                % parameter set
                % - from pset object
                if isa(varargin{5},'pset')
                    obj.parset = varargin{5};
                % - from np x 2 matrix
                elseif (isnumeric(varargin{5}) && size(varargin{5},2) == 2 ...
                            && size(varargin{5},1) == nv-1) 
                    obj.parset = pset.Box(varargin{5});
                else
                    error('PGSS:PGSS:inputError',...
                        'SET must be a box pset object or a %.0fx2 matrix',...
                        size(obj.A,3))
                end
                
                % checking dimensions of fcn_i
                if isa(obj.parfcn,'function_handle')
                    try
                        out = obj.parfcn(obj.parset.points(:,1));
                    catch err
                        error('PGSS:PGSS:inputError',...
                            'Function definition does not correspond with the parameter set')
                    end
                    if (length(out) ~= nv-1)
                        error('PGSS:PGSS:inputError',...
                            'FNC must be an anonymous function with %.0f output(s)',...
                            size(obj.A,3)-1)
                    end
                end

            % ------------------------------------------------------------
            % obj = pgss(sys,set,fnc) 
            elseif isa(varargin{1},'ss')
                % 
                if (nargin < 3)
                    error('PGSS:PGSS:notEnoughInputs', 'Not enough input arguments')
                end
                % matrices info from ss 3d object
                sys    = varargin{1};
                % system matrices
                obj.A  = sys.a;
                obj.B  = sys.b;
                obj.C  = sys.c;
                obj.D  = sys.d;
                
                % parameter functions fcn_i(p)
                if isa(varargin{3},'function_handle') || isnan(varargin{3})
                    obj.parfcn = varargin{3};
                else
                    error('PGSS:PGSS:inputError','FNC must be an anonymous function')
                end
                
                % parameter set
                % - from pset object
                if isa(varargin{2},'pset')
                    obj.parset = varargin{2};
                % - from np x 2 matrix
                elseif (isnumeric(varargin{2}) && size(varargin{2},2) == 2 ...
                            && size(varargin{2},1) == size(obj.A,3)-1) 
                    obj.parset = pset.Box(varargin{2});
                else
                    error('PGSS:PGSS:inputError',...
                        'SET must be a box pvset object or a %.0fx2 matrix',...
                        size(obj.A,3))
                end
                
                % checking dimensions of fcn_i
                if isa(obj.parfcn,'function_handle')
                    try
                        out = obj.parfcn(obj.parset.points(:,1));
                    catch err
                        error('PGSS:PGSS:inputError',...
                            'Function definition does not correspond with the parameter set')
                    end
                    if (length(out) ~= size(obj.A,3)-1)
                        error('PGSS:PGSS:inputError',...
                            'FNC must be an anonymous function of %.0f variable(s)',...
                            size(obj.a,3)-1)
                    end
                end
                
                % state, input and output names
                obj.StateName   = sys.StateName;
                obj.InputName   = sys.InputName;
                obj.OutputName  = sys.OutputName;
                obj.InputGroup  = sys.InputGroup;
                obj.OutputGroup = sys.OutputGroup;

            % ------------------------------------------------------------
            % obj = pgss(usys,pdG) --> conversion from a USS object
            elseif isa(varargin{1},'uss')
                usys   = varargin{1};
                
                % Getting fcn_i
                if (nargin > 1)
                    % gral model (copy function from other pass obj)
                    obj.parfcn = varargin{2}.parfcn;
                end
                
                % system matrices: reconstruction from the umat object
                ParName = fields(usys.Uncertainty); 
                nv = length(ParName);
                if (nargin > 1),
                   if (nv+1 ~= nsys(varargin{2}))
                       error('PASS:PASS:inputError',...
                           'USS object must have the same parameter number than PASS object');
                   end
                end

                ns = order(usys);
                [no,ni] = iosize(usys);
                as = zeros(ns,ns,nv+1);
                bs = zeros(ns,ni,nv+1);
                cs = zeros(no,ns,nv+1);
                ds = zeros(no,ni,nv+1);
                for ii = 1:nv+1
                    for jj = 1:nv
                        if (ii == 1 || jj ~= ii-1)
                            vpar.(ParName{jj}) = 0;
                        else
                            vpar.(ParName{jj}) = 1;
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

                % parameter set:
                if (nargin > 1)
                    % takes parset from other pass object
                    obj.parset = varargin{2}.parset;
                else
                    % built from uss info
                    parRange(size(as,2),2) = 0;
                    for ii = 1:size(as,2)
                        parRange(ii,:) = eval(sprintf('usys.Uncertainty.p%d.range',ii));
                    end
                    obj.parset = pset(parRange);
                end
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
