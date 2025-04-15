function varargout = lpvanalysis(varargin)

% LPVANALYSIS: General constraint analysis for LPV models
%
% LPVANALYSIS checks is the plant pdG satisfies the problem:
%	
%		min sum_i softConstraints_i(pdG)
%		s.t. hardConstraints_j
%		
% where:
%   - softConstraint_i: elements with soft constraints (e.g. || ||_2 < gamma)
%   - hardConstraint_j: elements with hard constraints (e.g. eigenvalue locations)
%
% each constraint includes the performance output, the disturbance inputs, 
% and factor, e.g: 
%   f1*||pdG(z1,w1)||_inf + f2*||pdG(z2,w2)||_2
% (see synConst.Gain, synConst.Poles for more details)
%	
% Use:
%   [constOut,obj,X] = lpvanalysis(pdG,constraints,[lyapfcn],[points],[options])
%
% Inputs:
%   - pdG:      the plant model that can be pass, ppss, pgss, psys, ss, tf, zpk
%   - constraint: synConst object specifying type and dist/perf signals
%   - lyapfcn:  Lyapunov function dependance (constant if not defined)
%               The following shortcuts can be used:
%               - cte: constant function  
%               - aff/affdX: for affine function derivative zero or nonzero
%               - pwa/pwadX: for pwa function derivative zero or nonzero
%               adding the last letter "i", the algorithm used the inverse
%               function, in some case with better numerical conditions
%               Otherwise a struct produce by createLyapFcn must be used 
%               with the function description
%   - points:   a (np x nv) matrices of points in which the constraints are
%               evaluated
%   - options:  a struct with settings produced by lpvsettings
%
% Outputs:
%   - constraint: synConst object with resulting constraints
%   - obj:      objective function value
%   - X:        Lyapunov matrix function

% fbianchi - 2021-07-05
% fbianchi - 2024-01-19 - rev

% ========================================================================
% Processing input arguments
%
if (nargin < 2)
    error('LPVANALYSIS:inputerror','Not enough input arguments.')
end
points = [];

switch nargin
    
    case 2
        % [constOut,obj] = lpvanalysis(pdG,constraints)
        [pdG,ftype] = checkSys(varargin{1}); 
        constraints = varargin{2};
        if isempty(pdG)
            error('LPVANALYSIS:inputerror','Invalid plant')
        end
        if ~isa(constraints,'synConst.synConst')
            error('LPVANALYSIS:inputerror','The second arguments must be a constraints')
        end
        lyapfcn = [];
        options = lpvsettings;
        
    case 3
        % [constOut,obj] = lpvanalysis(pdG,constraints,lyapFcn)
        [pdG,ftype] = checkSys(varargin{1}); 
        constraints = varargin{2};
        lyapfcn = varargin{3};
        if isempty(pdG)
            error('LPVANALYSIS:inputerror','Invalid plant')
        end
        if ~isa(constraints,'synConst.synConst')
            error('LPVANALYSIS:inputerror','The second arguments must be a constraints')
        end
        if strcmp(ftype,'lti') 
            if (size(pdG,3) == 1)
                warning('For LTI plants, the argument lyapFcn is ignored');
                lyapfcn = 'cte';
            else 
                error('LPVANALYSIS:inputerror',...
                    'In case of 3D plants, the parameter points must be provided')
            end
        end
        options = lpvsettings;

    case 4
        % [constOut,obj] = lpvanalysis(pdG,constraints,lyapFcn,points)
        [pdG,ftype] = checkSys(varargin{1}); 
        constraints = varargin{2};
        lyapfcn = varargin{3};
        points = varargin{4};
        if isempty(pdG)
            error('LPVANALYSIS:inputerror','Invalid plant')
        end
        if ~isa(constraints,'synConst.synConst')
            error('LPVANALYSIS:inputerror','The second arguments must be a constraints')
        end
        if strcmp(ftype,'lti') 
            if (size(pdG,3) == 1)
                warning('For 1D LTI plants, the argument lyapFcn and points are ignored');
                lyapfcn = 'cte';
            else
                error('LPVANALYSIS:inputerror',...
                    'In case of 3D plants, the number of cols in points must be equal to the number of plants')
            end
        end
        if strcmp(ftype,'lpv') && ischar(lyapfcn) && (strcmp(lyapfcn(1:3),'aff') || strcmp(lyapfcn(1:3),'pwa'))
            warning('For affine or PWA Lyapunov, the argument POINTS is ignored');
            points = [];
        end
        options = lpvsettings;
       
    case 5 
        % [constOut,obj] =  lpvanalysis(pdG,constraints,lyapFcn,points,options)
        [pdG,ftype] = checkSys(varargin{1}); 
        constraints = varargin{2};
        lyapfcn = varargin{3};
        points = varargin{4};
        options = varargin{5};
        if isempty(pdG)
            error('LPVANALYSIS:inputerror','Invalid plant')
        end
        if ~isa(varargin{2},'synConst.synConst')
            error('LPVANALYSIS:inputerror','The second arguments must be a constraints')
        end
        if strcmp(ftype,'lti')
            if (size(pdG,3) == 1)
                warning('For LTI plants, the argument lyapFcn is ignored');
                lyapfcn = 'cte';
            else
                error('LPVANALYSIS:inputerror',...
                    '3D plants are not supported')
            end
        end
        if strcmp(ftype,'lti') && (size(pdG,3) == 1)
            warning('LTI plants, the argument POINTS is ignored');
        elseif ischar(lyapfcn) && (strcmp(lyapfcn(1:3),'aff') || strcmp(lyapfcn(1:3),'pwa'))
            warning('For affine or PWA Lyapunov, the argument POINTS is ignored');
        end
        if ~isstruct(options)
            error('LPVANALYSIS:inputError',...
                'Invalid settings struct')
        end
        
    otherwise
        error('LPVANALYSIS:inputerror','Too many input arguments.')

end
% options.eigtol = 0;
options.penalty = 1e-9;
% option.solTol = 1e-4;

% %%%%%%%%%%%%%%%%%
% pdG
% constraints
% lyapfcn
% points
% %%%%%%%%%%%%%%%%%

% analyzed plant
if strcmp(ftype,'lpv')
    [~,~,ns,np,nv] = size(pdG);
    if isempty(points)
        lyapSet.LPVatPoints = ss(pdG);
        lyapSet.points = pdG.parset.points;
        
    elseif isnumeric(points) && (size(points,1) == np)
        lyapSet.LPVatPoints = ss(pdG,points);
        lyapSet.points = points;
        
    else
        error('LPVANALYSIS:inputError',...
            'POINTS must be a numeric matrix with %.0f rows',np)
    end
    lyapSet.parset = pdG.parset;

else
    ns = order(pdG(:,:,1));
    lyapSet.LPVatPoints = pdG;
    lyapSet.points = points;
end

% ------------------------------------------------------------------------
% List of constraints by type
constClasses         = arrayfun(@(x)class(x), constraints, 'UniformOutput', false);
constSet.Gain.idx    = find(ismember(constClasses,'synConst.Gain'));
constSet.Gain.idx    = constSet.Gain.idx(:).';
constSet.GainH2.idx  = find(ismember(constClasses,'synConst.GainH2'));
constSet.GainH2.idx  = constSet.GainH2.idx(:)';
constSet.Poles.idx   = find(ismember(constClasses,'synConst.Poles'));
constSet.Poles.idx   = constSet.Poles.idx(:)';
constSet.GainH2g.idx = find(ismember(constClasses,'synConst.GainH2g'));
constSet.GainH2g.idx = constSet.GainH2g.idx(:)';


% ------------------------------------------------------------------------
% Lyapunov variables 
%
lyapSet.inv = false; % don't use inverse Lyapunov function
if isempty(lyapfcn) || (ischar(lyapfcn) && strcmp(lyapfcn(1:3),'cte'))
    % constant Lyapunov function by default
    lyapSet.type = 'cte';
    lyapSet.parfcn  = 0;
    lyapSet.dparfcn = 0;
    if ischar(lyapfcn) && strcmp(lyapfcn(end),'i')
        lyapSet.inv = true;
    end

elseif (nv == 1)
    % when pdG is LTI, Lyapunov function is constant
    lyapSet.type = 'cte';
    lyapSet.parfcn  = 0;
    lyapSet.dparfcn = 0;
    if ischar(lyapfcn) && ~strcmp(lyapfcn,'cte')
        warnning('LPVANALYSIS','For classes ss, tf and zpk, only constant Lyapunov fcn')
    end
    if ischar(lyapfcn) && strcmp(lyapfcn(end),'i')
        lyapSet.inv = true;
    end

elseif ischar(lyapfcn) && strcmp(lyapfcn(1:3),'pwa')
    % PWA Lyapunov functions
    if (length(lyapfcn) >= 5) && strcmp(lyapfcn(4:5),'dX')
        % with dX~=0
        lyapSet.type = 'pwadX';
    else
        % with dX==0
        lyapSet.type = 'pwa';
    end
    lyapSet.parfcn  = 0;
    lyapSet.dparfcn = 0;
    if strcmp(lyapfcn(end),'i')
        lyapSet.inv = true;
    end

elseif ischar(lyapfcn) && strcmp(lyapfcn(1:3),'aff')
    % Affine Lyapunov functions
    if (length(lyapfcn) >= 5) && strcmp(lyapfcn(4:5),'dX')
        % with dX~=0
        lyapSet.type = 'affdX';
        lyapSet.parfcn  =@(p) p;
        lyapSet.dparfcn =@(p) 0;
    else
        % with dX==0
        lyapSet.type = 'aff';
        lyapSet.parfcn  =@(p) p;
        lyapSet.dparfcn = 0;
    end
    if strcmp(lyapfcn(end),'i')
        lyapSet.inv = true;
    end
    
elseif isstruct(lyapfcn)
    % general case
    aux = isnumeric(lyapfcn.parfcn) && (lyapfcn.parfcn == 0);
    aux = aux + (isnumeric(lyapfcn.parfcn) && (lyapfcn.dparfcn == 0));
    if (aux == 2)
        lyapSet.type = 'cte';
    else
        lyapSet.type = 'gral';
                                % computable controller
    end
    lyapSet.parfcn  = lyapfcn.parfcn;
    lyapSet.dparfcn = lyapfcn.dparfcn;
    if isfield(lyapfcn,'inv')
        lyapSet.inv = lyapfcn.inv;
    end
else
    error('LPVANALYSIS:inputError:LyapunovFunction',...
        'Invalid Lyapunov Function')
end
% optimization variables
switch lyapSet.type(1:3)
    case 'cte'
        % constant Lyapunov function
        lyapSet.X = sdpvar(ns);
        
    case {'pwa', 'aff'}
        % PWA/affine Lyapunov function
        lyapSet.X = sdpvar(ns,ns,nv);
        
    otherwise
        % parameter dependent Lyapunov function
        if isa(lyapSet.parfcn,'function_handle')
            np = size(pdG.parset);
            nf = length(lyapSet.parfcn(ones(np,1))) + 1;
            lyapSet.X = sdpvar(ns,ns,nf);
        else
            lyapSet.X = sdpvar(ns,ns);
        end
end

% ========================================================================
% LMI constraintss:
lmis = []; obj = 0;

% Gain - Hinf constrains
for ic = constSet.Gain.idx

    if isinf(constraints(ic).bound)
        constraints(ic).bound = sdpvar(1);
    end
    [lmiHinf,objHinf] = constHinf2LMI_a(pdG,lyapSet,constraints(ic),options);
    lmis = [lmis, lmiHinf];
    obj  = obj + objHinf;
    
end

% Gain - H2 constrains
for ic = constSet.GainH2.idx

    if isinf(constraints(ic).bound)
        constraints(ic).bound = sdpvar(1);
    end
    [lmiH2,objH2] = constH22LMI_a(pdG,lyapSet,constraints(ic),options);
    lmis = [lmis, lmiH2];
    obj = obj + objH2;
    
end

% Generalized Gain - H2 constrains
for ic = constSet.GainH2g.idx

    if isinf(constraints(ic).bound)
        constraints(ic).bound = sdpvar(1);
    end
    [lmiH2g,objH2g] = constH2g2LMI_a(pdG,lyapSet,constraints(ic),options);
    lmis = [lmis, lmiH2g];
    obj = obj + objH2g;
    
end

% Poles constrains
for ic = constSet.Poles.idx

    lmiPoles = constPoles2LMI_a(pdG,lyapSet,constraints(ic),options);
    lmis = [lmis, lmiPoles];
    
end

% common to all constraintss
[lmiPos,objPos] = constPos2LMI_a(pdG,lyapSet,options);
lmis = [lmis, lmiPos];
obj = obj + objPos;

% -------------------------------------------------------------------------
% Solving optimization problem
if strcmp(options.solver, 'sedumi')
    solverOpts = sdpsettings('solver', 'sedumi',...
                  'verbose', options.verb,...
                  'dualize', options.dualize,...
                  'sedumi.eps', options.solTol,...
                  'sedumi.cg.qprec', 1,...
                  'sedumi.stepdif', 2,...
                  'sedumi.maxiter', options.maxiter,...
                  'removeequalities', options.removeequalities);
    
elseif strcmp(options.solver, 'sdpt3')
    solverOpts = sdpsettings('solver', 'sdpt3',...
                  'sdpt3.steptol', option.solTol,...
                  'sdpt3.maxit', options.maxiter,...
                  'verbose', options.verb,...
                  'dualize', options.dualize,...
                  'removeequalities', options.removeequalities);

elseif strcmp(options.solver, 'mosek')
    solverOpts = sdpsettings('solver', 'mosek',...
                  'verbose', options.verb,...
                  'dualize', options.dualize,...
                  'removeequalities', options.removeequalities);
end
diagnostic = optimize(lmis, obj, solverOpts);

% -------------------------------------------------------------------------
% Results
if ~((diagnostic.problem == 0) || (diagnostic.problem == 4))
    % Infeasible
    
    fprintf('Infeasible analysis problem\n')
    fprintf('\tSolver says: %s\n',yalmiperror(diagnostic.problem));
                
    obj = Inf; 
    constSets = [constSet.Gain.idx constSet.GainH2.idx constSet.GainH2g.idx];
    for ii = constSets
        if isa(constraints(ii).bound,'sdpvar')
            constraints(ii).bound = inf;
        end
    end

else
    % Feasible
    
    obj  = value(obj);
    % summarizing results
    n = length(constraints); jj = 1;
    for ii = 1:n
        [s1,s2,s3] = char(constraints(ii));
        if ~isempty(s1)
            str1{jj} = s1; str2{jj} = s2; str3{jj} = s3;
            jj = jj + 1;
        end
    end
    idx_min = find(ismember(str1,'min'));
    idx_st  = find(ismember(str1,'s.t.'));
    
    fprintf('\n------------------------------------------------------------------\n')
    fprintf('Solved the following analysis problem:\n\n')
    if any(idx_min)
        fprintf('\tminimize:\n')
        for ii = idx_min(1:end-1)
            fprintf('\t\t%s + \t\t\t(%s)\n',str2{ii},str3{ii})
        end
        ii = idx_min(end);
        fprintf('\t\t%s  \t\t\t(%s)\n',str2{ii},str3{ii})
    end
    fprintf('\n')
    if any(idx_st)
        fprintf('\tsubject to:\n')
        for ii = idx_st
            if strcmp(str3{ii},'')
                fprintf('\t\t%s\n',str2{ii})
            else
                fprintf('\t\t%s \t(%s)\n',str2{ii},str3{ii})
            end
        end
    end
    fprintf('\nwhere:\n')
    if any(idx_min)
        for ii = idx_min(1:end)
            classConst = class(constraints(ii));
            switch classConst
                case {'synConst.Gain', 'synConst.GainH2', 'synConst.GainH2g'}
                    constraints(ii).bound = sqrt(value(constraints(ii).bound));
            end
            [~,str2,str3] = char(constraints(ii));
            fprintf('\t\t%s, \t(%s)\n',str2,str3);
        end
    end
    
    fprintf('\n')
    fprintf('Lyapunov functions\n');
    if lyapSet.inv
        strXcl = 'X^{-1}';
    else
        strXcl = 'X';
    end        
    
    if strcmp(lyapSet.type,'gral')
        fprintf('\t%s:  %s\n',strXcl,func2str(lyapSet.parfcn));
    elseif ischar(lyapSet.type) && strcmp(lyapSet.type(1:3),'aff')
        fprintf('\t%s:  %s\n',strXcl,'affine');
    elseif ischar(lyapSet.type) && strcmp(lyapSet.type(1:3),'pwa')
        fprintf('\t%s:  %s\n',strXcl,'PWA');
    else
        fprintf('\t%s:  %s\n',strXcl,'constant');
    end
    if strcmp(ftype,'lti') || isempty(pdG.parset.rate) || all(all(pdG.parset.rate)) == 0
        str_dpdt = ' (but dp/dt=0)';
    else
        str_dpdt = '';
    end        
    if strcmp(lyapSet.type,'gral')
        if isa(lyapSet.dparfcn,'function_handle')
            fprintf('\td%s: %s%s\n',strXcl,func2str(lyapSet.dparfcn),str_dpdt);
        else
            fprintf('\td%s: %s%s\n',strXcl,num2str(lyapSet.dparfcn),str_dpdt);
        end
    elseif ischar(lyapSet.type) && strcmp(lyapSet.type,'affdX')
        fprintf('\td%s: %s%s\n',strXcl,'[cte]*dp/dt',str_dpdt);
    elseif ischar(lyapSet.type) && strcmp(lyapSet.type,'pwadX')
        fprintf('\td%s: %s%s\n',strXcl,'[cte]*dp/dt',str_dpdt);
    else
        fprintf('\td%s: 0\n',strXcl);
    end
    
    % Gain - Hinf constrains
    for ic = constSet.Gain.idx
        fprintf('\nInfinity-norm at frozen points\n');
        
        nhinf = norm(lyapSet.LPVatPoints(constraints(ic).outputs,...
                     constraints(ic).inputs,:),inf);
                 
        for ii = 1:length(nhinf)
            fprintf('\tpoint(%.0f): norm_inf = %6.4f\n', ii, nhinf(ii));
        end

    end

    % Gain - H2 constrains
    for ic = [constSet.GainH2.idx constSet.GainH2g.idx]

        fprintf('\n2-norm at frozen points\n');
        
        n2 = norm(lyapSet.LPVatPoints(constraints(ic).outputs,...
                     constraints(ic).inputs,:),2);
                 
        for ii = 1:length(n2)
            fprintf('\tpoint(%.0f): norm_2 = %6.4f\n', ii, n2(ii));
        end

    end
    fprintf('------------------------------------------------------------------\n\n')

    
    
end

% function outputs
varargout{1} = constraints;
varargout{2} = obj;
Xopt = value(lyapSet.X);
if isa(lyapSet.type,'gral')
    X = pgss([],[],[],Xopt,pdG.parset,lyapSet.parfcn);
elseif ischar(lyapSet.type) && strcmp(lyapSet.type,'affdX')
    X = pass([],[],[],Xopt,pdG.parset);
elseif ischar(lyapSet.type) && strcmp(lyapSet.type,'pwadX')
    X = ppss([],[],[],Xopt,pdG.parset);
else
    X = Xopt;
end
varargout{3} = X;



% ========================================================================
% Sub-functions

function [sys,ftype] = checkSys(arg)

LTIclasses = {'ss','tf','zpk'};
LPVclasses = {'pass','ppss','pgss','pcss'};

if ismember(class(arg),LPVclasses)
    sys = arg;
    ftype = 'lpv';
    
elseif ismember(class(arg),LTIclasses)
    if isa(arg,'ss')
        sys = ss(arg);
    elseif (isa(arg,'zpk') || isa(arg,'tf'))
        sys = ss(minreal(ss(arg)));
    end
    ftype = 'lti';

else
    sys = [];
    ftype = 'none';

end

