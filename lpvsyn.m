function varargout = lpvsyn(pdG,varargin)

% LPVSYN: General controller synthesis for LPV models
%
% LPVSYN finds a controller pdK with input MEAS and output CTRL that solves
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
%   [pdK,constraint,obj] = LPVSYN(pdG,r)
%   [pdK,constraint,obj] = LPVSYN(pdG,meas,ctrl,[constraint],[lyapFcn],[ctrlFcn],[options])
%
% Inputs:
%   - pdG:      the augmented plant that can be pass, ppss, pgss, psys, ss, tf, zpk
%   - r:        1x2 vector with i/o size of the controller r=[ny nu]
%   - meas:     measure signals (names or indices)
%   - ctrl:     control signals (names or indices)
%   - constraint: synConst object specifying type and dist/perf signals
%   - lyapFcn:  Lyapunov function dependance (constant if not defined)
%   - ctrlFcn:  by default, the resulting controller will have the same dependance 
%               than the plant. This can be changed with ctrlFcn according
%               to:
%               * 0: robust controller (only for state feedback)
%               * -1: force dk=0 
%               * vector of indices corresponding to the terms of pdG
%               included in the controller dependence (the same dependance 
%               than the plant if not defined) 
%               * struct:
%                   * ctrlFcn.idxfcn: index vector of the terms of pdG included
%                   * ctrlFcn.XArYbnd: bound on term XArY
%                   * ctrlFcn.dk: negative value force dk=0
%                   * ctrlFcn.pdIn: set the aux variables Bk & Dk cte for
%                           case in which c2 & d21 are parameter dependent
%                   * ctrlFcn.pdOut: set the aux variables Ck & Dk cte for
%                           case in which b2 & d12 are parameter dependent
%   - options:  a struct with settings produced by LPVSETTINGS
%
% Outputs:
%   - pdK:      Controller. In general, pdK has the same parameter
%               dependency than pdG. Therefore, pdK can be a psys, pass, 
%               ppss or pcss depending on the pdG and LyapFcn.
%   - constraint: synConst object with resulting constraints
%   - obj:      objective function value

% fbianchi - 2020-07-07


% ========================================================================
% Processing input arguments
switch nargin
    case 1
        error('LPVSYN:inputerror','insufficient input arguments')
        
    case {2, 3}
        % only plant & i/o control map
        constraint = [];
        lyapfcn    = [];
        ctrlfcn    = [];
        options    = lpvsettings;
        
    case 4
        % plant, i/o control map & constraint
        constraint = varargin{3};
        lyapfcn    = [];
        ctrlfcn    = [];
        options    = lpvsettings;
        
    case 5
        % plant, i/o control map, constraint & lyapnov function
        constraint = varargin{3};
        lyapfcn    = varargin{4};
        ctrlfcn    = [];
        options    = lpvsettings;
        
    case 6
        % plant, i/o control map, constraint, lyapnov & control function
        constraint = varargin{3};
        lyapfcn    = varargin{4};
        ctrlfcn    = varargin{5};
        options    = lpvsettings;
        
    otherwise
        % plant, i/o control map, constraint, lyapnov, control function, &
        % options
        constraint = varargin{3};
        lyapfcn    = varargin{4};
        ctrlfcn    = varargin{5};
        options    = varargin{6};

end

% ------------------------------------------------------------------------
% Identifying synthesis according to the input parameters
%
% Synthesis structure:
%   - synSet.type = prFull/prA/basic -> type of ctrl variable elimination
%   - synSet.pdGclass                -> original class of pdG
%   - synSet.LPVatPoints             -> pdG evaluated at the points of pset
%
% Lyapunov structure:
%   - lyapSet.type:         'cte'/'pwa'/'aff(X,dX,Y,dY)'/'gral'
%   - lyapSet.parfcnX:      anonymous function
%   - lyapSet.dparfcnX:     anonymous function
%   - lyapSet.parfcnY:      anonymous function
%   - lyapSet.dparfcnY:     anonymous function
%   - lyapSet.Y,X,dY,dX:    lyapunov variables
%   
% Controller structure:
%   - ctrlSet.fback: state/output    -> type of feedback
%   - ctrlSet.sch: gs/robust         -> type of ctrl (robust only for state
%                                       feedback)
%   - ctrlSet.Ak,Bk,Ck,Dk: controller variables
%   - ctrlSet.extFcn: 0 = ctrller dependence equals to pdG, 
%                     1 = given the indices, 2 = given functions
%   - ctrlSet.parset: parameter set 
%   - ctrlSet.parfcn: parameter functions
%   - ctrlSet.interp: matrices interpolation -> aff/pwa/cte
%   - ctrlSet.nk: number of LTI model used in the description
%   - ctrlSet.idxfcn:  vector with the indices of the terms included in the
%                      controller
%   - ctrlSet.XArYbnd: bound on X*E and Y*F
%   - ctrlSet.XEbnd,XEbnd: decision variables 
%   These three variables are used when the controller parameter dependence 
%   is imposed
%   - ctrlSet.pdIn: if true then the  auxiliary variables dk and bk are 
%                   parameter dependent
%   - ctrlSet.pdOut: if true then the  auxiliary variables dk and ck are 
%                   parameter dependent
%

% ---------------------------> Plant pdG <--------------------------------
%
synSet.pdGclass = class(pdG);
if (exist('ispsys.m','file') == 2) && ispsys(pdG)
    % case lmitool object
    typ = psinfo(pdG);
    if strcmp(typ,'aff')
        pdG = pass(pdG);
        synSet.pdGclass = 'psys-aff';
    else
        pdG = ppss(pdG);
        synSet.pdGclass = 'psys-pol';
    end
    
elseif ~isa(pdG,'p_ss') 
    % case LTI models
    if isa(pdG,'ss') 
        pdG = ppss(pdG);
        
    elseif (isa(pdG,'zpk') || isa(pdG,'tf'))
        pdG = ppss(minreal(ss(pdG)));
        
    else
        error('LPVSYN:inputerror','Invalid plant pdG')
    end

end
%
% ---------------------> meas & ctrl signals <----------------------------
%
if (nargin == 2)
    % second argument is [ny nu]
    iomap = varargin{1};
    [pdG,ny,nu] = setIOctrl(pdG,iomap);

elseif (nargin >= 3)
    % 2nd & 3rd arguments are the list of meas & ctrl
    [pdG,ny,nu] = setIOctrl(pdG,varargin{1},varargin{2});
    
end
% check if the control channel is parameter indenpendent
[ioDiag,ioMsg] = ispd(pdG('meas','ctrl'));
% moved below Lyapunov functions because it isn't necessary if Lyapunov
% function is general
% if any(ioDiag(2:3))
%     error(ioMsg)
% end
synSet.pdMeasCtrl = ioDiag(2:3);
% plant evaluated at each point of the PARSET.POINTS with control channel
% information
synSet.LPVatPoints = ss(pdG);
% system dimensions
ns = order(pdG);        % order
np = npar(pdG);         % number of parameters
nv = nsys(pdG);         % number of LTI model in the description
%
% -------------------------> control scheme <-----------------------------
%
% gain scheduling by default
ctrlSet.sch = 'gs'; 
% checking if it's state-feedback
if (ny == 0)
    % state feedback
    ctrlSet.fback = 'state';
else
    % output feedback
    ctrlSet.fback = 'output';
end 
%
% ----------------------------> constraints <-----------------------------
%
if isempty(constraint)
    % (Default) L2-norm constraint with nw & nz deduced from ny & nu
    % 
    if (ny == 0)
        error('LPVSY:inputError:constraint',...
            'NY must be greater than zero')
    end
    [no,ni] = iosize(pdG);
    nz = no - ny;
    nw = ni - nu;
    constraint = synConst.Gain(1:nw,1:nz);
    
elseif ~isa(constraint,'synConst.synConst')
    error('LPVSY:inputError:constraint',...
        'Invalid constraint')
end
% List of constraint by type
% nConst              = length(constraint);
constClasses = arrayfun(@(x)class(x), constraint, 'UniformOutput', false);
constSet.Gain.idx = find(ismember(constClasses,'synConst.Gain'));
constSet.Gain.idx = constSet.Gain.idx(:).';
nConst(1) = length(constSet.Gain.idx);
constSet.GainH2.idx = find(ismember(constClasses,'synConst.GainH2'));
constSet.GainH2.idx = constSet.GainH2.idx(:)';
nConst(2) = length(constSet.GainH2.idx);
constSet.Poles.idx = find(ismember(constClasses,'synConst.Poles'));
constSet.Poles.idx = constSet.Poles.idx(:)';
nConst(3) = length(constSet.Poles.idx);
constSet.GainH2g.idx = find(ismember(constClasses,'synConst.GainH2g'));
constSet.GainH2g.idx = constSet.GainH2g.idx(:)';
nConst(4) = length(constSet.GainH2g.idx);

if (nConst(1) == 1) && all(nConst(2:end) == 0)
    % In case of having on 1 Hinf constraints => Full projection
    if (exist('klmi','file') == 2)
        synSet.type = 'prFull'; % to avoid using klmi
    else
        synSet.type = 'basic';
    end
    
    % ToDo: include projected version of H2 constraint
    % elseif (nConst(2) == 1) && all(nConst([1 3:end]) == 0)
    % In case of having on 2 H2 constraints => Full projection H2
    
    % synSet.type = 'prA';
else
    % otherwise, change of variable is used
    synSet.type = 'basic';
end

if isa(pdG,'pass') || isa(pdG,'pgss') || strcmp(ctrlSet.fback,'state') 
    % for class pass/pgss only basic preserves the controller dependency
    % so the type of elimination is over-written
    synSet.type = 'basic';
end
% WARNNING: this flag may change according to other synthesis options
    
% --------------------> Lyapunov variables <-----------------------------
%
lyapSet.parset = pdG.parset;
% to deal with inf parameter rates -> implies the Lyapunov function is not
% dependent of this parameters
if ~isempty(lyapSet.parset.rate)
    bnd = lyapSet.parset.rate;
    idx = isinf(bnd(:,1)) | isinf(bnd(:,2));
else
    idx = false(np,1);  % idx=true implies that fcn is independent of the 
                        % corresponding parameter
end

if isempty(lyapfcn) || (ischar(lyapfcn) && strcmp(lyapfcn,'cte'))
    % constant Lyapunov function by default
    lyapSet.type = 'cte';
    lyapSet.parfcnX = 0;
    lyapSet.dparfcnX = 0;
    lyapSet.parfcnY = 0;
    lyapSet.dparfcnY = 0;

elseif ismember(synSet.pdGclass,{'ss','tf','zpk'})
    % when pdG is LTI, Lyapunov function is constant
    lyapSet.type = 'cte';
    lyapSet.parfcnX = 0;
    lyapSet.dparfcnX = 0;
    lyapSet.parfcnY = 0;
    lyapSet.dparfcnY = 0;
    if ischar(lyapfcn) && ~strcmp(lyapfcn,'cte')
        warnning('LPVSYN','For classes ss, tf and zpk, only constant Lyapunov fcn')
    end

elseif ischar(lyapfcn) && strcmp(lyapfcn(1:3),'pwa')
    % PWA Lyapunov functions uses affine functions in each simplex
    if ~isa(pdG.parset, 'pset.Gral') || ~isa(pdG,'ppss')
        error('LPVSYN:inputError:LyapunovFunction',...
        'PWA Lyapunov functions requires a ppss plant with a pset.Gral parameter set');
    end
    lyapSet.type = lyapfcn;
    switch lyapSet.type
        case {'pwa'}
            nfY = size(synSet.LPVatPoints,3);
            nfX = size(synSet.LPVatPoints,3);
        case {'pwaX','pwadX'}
            nfY = 1;
            nfX = size(synSet.LPVatPoints,3);
        case {'pwaY','pwadY'}
            nfY = size(synSet.LPVatPoints,3);
            nfX = 1;
        otherwise
            error('LPVSYN:inputError:LyapunovFunction',...
                    'Invalid Lyapunov Function');
    end
    lyapSet.parfcnX  = 0;
    lyapSet.dparfcnX = 0;
    lyapSet.parfcnY  = 0;
    lyapSet.dparfcnY = 0;
    synSet.type = 'basic';
    

elseif ischar(lyapfcn) && strcmp(lyapfcn(1:3),'aff')
    % affine cases
    if ~isa(pdG,'pass')
        error('LPVSYN:inputError:LyapunovFunction',...
             'Affine Lyapunov functions are not available for class %s',class(pdG));
    end
    lyapSet.type = lyapfcn;
    auxE = eye(np); auxE(idx,:) = [];
    switch lyapSet.type
        case 'aff'
            lyapSet.parfcnX  =@(p) p(~idx);
            lyapSet.dparfcnX = 0;
            lyapSet.parfcnY  =@(p) p(~idx);
            lyapSet.dparfcnY = 0;
        case 'affX'
            lyapSet.parfcnX  =@(p) p(~idx);
            lyapSet.dparfcnX = 0;
            lyapSet.parfcnY  = 0;
            lyapSet.dparfcnY = 0;
        case 'affdX'
            lyapSet.parfcnX  =@(p) p(~idx);
            lyapSet.dparfcnX =@(p) auxE;
            lyapSet.parfcnY  = 0;
            lyapSet.dparfcnY = 0;
        case 'affY'
            lyapSet.parfcnX  = 0;
            lyapSet.dparfcnX = 0;
            lyapSet.parfcnY  =@(p) p(~idx);
            lyapSet.dparfcnY = 0;
        case 'affdY'
            lyapSet.parfcnX  = 0;
            lyapSet.dparfcnX = 0;
            lyapSet.parfcnY  =@(p) p(~idx);
            lyapSet.dparfcnY =@(p) auxE;
        otherwise
            error('LPVSYN:inputError:LyapunovFunction',...
                     'Invalid Lyapunov Function');
    end
    
elseif isstruct(lyapfcn)
    % general case
    bnd = isnumeric(lyapfcn.parfcnX) && (lyapfcn.parfcnX == 0);
    bnd = bnd + (isnumeric(lyapfcn.parfcnX) && (lyapfcn.dparfcnX == 0));
    bnd = bnd + (isnumeric(lyapfcn.parfcnY) && (lyapfcn.parfcnY == 0));
    bnd = bnd + (isnumeric(lyapfcn.dparfcnY) && (lyapfcn.dparfcnY == 0));
    if (bnd == 4)
        lyapSet.type = 'cte';
    else
        lyapSet.type = 'gral';
        synSet.type  = 'basic'; % overwrite in order to have an online 
                                % computable controller
    end
    lyapSet.parfcnX  = lyapfcn.parfcnX;
    lyapSet.dparfcnX = lyapfcn.dparfcnX;
    lyapSet.parfcnY  = lyapfcn.parfcnY;
    lyapSet.dparfcnY = lyapfcn.dparfcnY;
else
    error('LPVSYN:inputError:LyapunovFunction',...
        'Invalid Lyapunov Function')
end
% lyapSet.parset = pdG.parset;
% optimization variables
if strcmp(ctrlSet.fback,'output')
    % output feedback
    switch lyapSet.type(1:3)
        case 'cte'
            % constant Lyapunov function
            lyapSet.X = sdpvar(ns);
            lyapSet.Y = sdpvar(ns); 
        
        case 'pwa'
            % PWA Lyapunov function
            lyapSet.Y = sdpvar(ns,ns,nfY);
            lyapSet.X = sdpvar(ns,ns,nfX);
        
        otherwise
            % parameter dependent Lyapunov function
            if isa(lyapSet.parfcnX,'function_handle')
                nfx = length(lyapSet.parfcnX(ones(np,1))) + 1;
                lyapSet.X = sdpvar(ns,ns,nfx);
            else
                lyapSet.X = sdpvar(ns,ns);
            end
            if isa(lyapSet.parfcnY,'function_handle')
                nfy = length(lyapSet.parfcnY(ones(np,1))) + 1;
                lyapSet.Y = sdpvar(ns,ns,nfy);
            else
                lyapSet.Y = sdpvar(ns,ns);
            end
    end
else
    % state feedback
    switch lyapSet.type
        case 'cte'
            % constant Lyapunov function
            lyapSet.X = [];
            lyapSet.Y = sdpvar(ns); 
        
        case 'pwa'
            % PWA Lyapunov function
            lyapSet.X = [];
            lyapSet.Y = sdpvar(ns,ns,nv); 
        
        otherwise
            % parameter dependent Lyapunov function
            lyapSet.X = [];
            if isa(lyapSet.parfcnY,'function_handle')
                nfy = length(lyapSet.parfcnY(ones(np,1))) + 1;
                lyapSet.Y = sdpvar(ns,ns,nfy);
            else
                lyapSet.Y = sdpvar(ns,ns);
            end
    end
end

% ----------------------> Controller variables <--------------------------
%
ctrlSet.parset = pdG.parset;
ctrlSet.parfcn = NaN;
ctrlSet.extFcn = 0;
ctrlSet.nk = nsys(pdG);
ctrlSet.idxfcn = 1:ctrlSet.nk;
ctrlSet.XArYbnd = 0;
% ctrlSet.interp = 'aff';
ctrlSet.dk = true;
ctrlSet.pdIn = true;
ctrlSet.pdOut = true;
ctrlSet.class = synSet.pdGclass;
if isa(pdG,'pgss')
    ctrlSet.parfcn = pdG.parfcn;
    ctrlSet.interp = 'aff';
    
elseif isa(pdG,'pass')
    ctrlSet.interp = 'aff';
    
elseif isa(pdG,'ppss')
    ctrlSet.interp = 'pwa';
    
else
    error('LPVSYN:inputError:CtrlVariables','Invalid pdG')
    
end
if isempty(ctrlfcn)
    % the parameter dependency is copied from pdG
    ctrlSet.extFcn = 0;
    
else
    % the controller dependency is given
    synSet.type = 'basic';
    if isnumeric(ctrlfcn)
        if (ctrlfcn == 0)
            % robust state feedback
            if strcmp(ctrlSet.fback,'output')
                error('LPVSYN:inputError:CtrlVariables',...
                    'The robust controller only available for state feedback')
            end
            ctrlSet.parfcn = ctrlfcn;
            ctrlSet.sch = 'robust';
            ctrlSet.nk = 1;
            ctrlSet.idxfcn = 1;
            ctrlSet.interp = 'cte';
            
        elseif (ctrlfcn < 0)
            % dk = 0
            ctrlSet.dk = false;   
            
        elseif all(ctrlfcn > 0)
            if any(~ismember(1:ctrlSet.nk,ctrlfcn))
                % only indices in ctrlfcn -> Procedure 2
                if ~(isa(pdG,'pass') || isa(pdG,'pgss'))
                    error('LPVSYN:inputError:CtrlVariables',...
                        'The controller parameter dependence only can be imposed if pdG is pass or pgss')
                end
                ctrlSet.extFcn = 1;
                ctrlSet.idxfcn = ctrlfcn;
                if isa(pdG,'pgss')
                    ctrlSet.parfcn = pdG.parfcn;
                end
                ctrlSet.nk = length(ctrlfcn);
            end
        
        else
            error('LPVSYN:inputError:CtrlVariables','Invalid CtrlFcn')
        end
        
    elseif isstruct(ctrlfcn)
        
%         TODO: allow using different functions for the controller
%         if isfield(ctrlfcn,'parfcn') 
%             if isa(ctrlfcn.parfcn,'function_handle')
%                 ctrlSet.parfcn = ctrlfcn.parfcn;
%             else
%                 error('LPVSYN:inputError:CtrlVariables','Invalid CtrlFcn')
%             end
%             try
%                 auxpar = ones(np,1);
%                 auxfcn = ctrlfcn.parfcn(auxpar);
%                 ctrlSet.nk = 1 + size(auxfcn,1);
%                 ctrlSet.idxfcn = 1:ctrlSet.nk;
%             catch
%                 error('LPVSYN:inputError:CtrlVariables','Invalid CtrlFcn')
%             end
%             ctrlSet.class = 'pgss';
%             ctrlSet.interp = 'aff';
%         end
        
        if isfield(ctrlfcn,'idxfcn') && any(~ismember(1:ctrlSet.nk,ctrlfcn.idxfcn))
            % struct with indices and bound
            % used when the controller dependence is imposed
            
            if ~(isa(pdG,'pass') || isa(pdG,'pgss'))
                error('LPVSYN:inputError:CtrlVariables',...
                    'The controller parameter dependence only can be imposed if pdG is pass or pgss')
            end
            if isfield(ctrlfcn,'XArYbnd') && isnumeric(ctrlfcn.XArYbnd)
                % with bound on ||XArY|| -> Procedure 1
                ctrlSet.XArYbnd = ctrlfcn.XArYbnd;
            else
                % without bound on ||XArY|| -> Procedure 2
                ctrlSet.XArYbnd = 0;    
            end
            ctrlSet.extFcn = 1;
            ctrlSet.idxfcn = ctrlfcn.idxfcn;
            if isa(pdG,'pgss')
                ctrlSet.parfcn = pdG.parfcn;
            end
            ctrlSet.nk = length(ctrlfcn.idxfcn);
            
        end
        
        if isfield(ctrlfcn,'dk')
            if isnumeric(ctrlfcn.dk) && isscalar(ctrlfcn.dk)
                ctrlSet.dk = (ctrlfcn.dk > 0);
            else
                error('LPVSYN:inputError:CtrlVariables',...
                    'ctrlfcn.dk must be a numeric scalar')
            end
                
        end
        
        if isfield(ctrlfcn,'pdIn')
            if (isnumeric(ctrlfcn.pdIn) || islogical(ctrlfcn.pdIn)) && isscalar(ctrlfcn.pdIn)
                ctrlSet.pdIn = (ctrlfcn.pdIn > 0);
            else
                error('LPVSYN:inputError:CtrlVariables',...
                    'ctrlfcn.pdIn must be a numeric scalar')
            end
                
        end 
        
        if isfield(ctrlfcn,'pdOut')
            if (isnumeric(ctrlfcn.pdOut) || islogical(ctrlfcn.pdOut)) && isscalar(ctrlfcn.pdOut)
                ctrlSet.pdOut = (ctrlfcn.pdOut > 0);
            else
                error('LPVSYN:inputError:CtrlVariables',...
                    'ctrlfcn.pdOut must be a numeric scalar')
            end
            
        end
        
    else
        error('LPVSYN:inputError:CtrlVariables','Invalid CtrlFcn')
    
    end
end    

% controller variables
switch synSet.type
    % in case of synSet.sch = robust -> nk = 1
    case 'basic'
        if strcmp(ctrlSet.fback,'output')    % output feedback
            ctrlSet.A = sdpvar(ns,ns,ctrlSet.nk,'full'); % Ak
            % Bk
            if (ctrlSet.pdIn)
                % parameter dependent input channel
                nkB = ctrlSet.nk;
                ctrlSet.B = sdpvar(ns,ny,nkB,'full');
            else
                nkB = 1;
                auxB = sdpvar(ns,ny,'full');
                if strcmp(ctrlSet.interp,'pwa')
                    ctrlSet.B = repmat(auxB,1,1,ctrlSet.nk);
                else
                    ctrlSet.B = cat(3,auxB,zeros(ns,ny,ctrlSet.nk-nkB));
                end
            end
            % Ck
            if (ctrlSet.pdOut)
                % parameter dependent output channel
                nkC = ctrlSet.nk;
                ctrlSet.C = sdpvar(nu,ns,nkC,'full');
            else
                nkC = 1;
                auxC = sdpvar(nu,ns,1,'full');
                if strcmp(ctrlSet.interp,'pwa')
                    ctrlSet.C = repmat(auxC,1,1,ctrlSet.nk);
                else
                    ctrlSet.C = cat(3,auxC,zeros(nu,ns,ctrlSet.nk-nkC));
                end
            end
            if (ctrlSet.dk)
                nkD = min(nkB,nkC);
                if (nkD == ctrlSet.nk)
                    ctrlSet.D = sdpvar(nu,ny,ctrlSet.nk,'full');
                else
                    auxD = sdpvar(nu,ny,1,'full');
                    if strcmp(ctrlSet.interp,'pwa')
                        ctrlSet.D = repmat(auxD,1,1,ctrlSet.nk);
                    else
                        ctrlSet.D = cat(3,auxD,zeros(nu,ny,ctrlSet.nk-nkD));
                    end
                end
            else
                ctrlSet.D = zeros(nu,ny,ctrlSet.nk);
            end
            if (ctrlSet.extFcn == 1)
                ctrlSet.XEbnd = sdpvar(1);
                ctrlSet.YFbnd = sdpvar(1);
            end
        else
            % state feedback (Dk is the state feedback)
            ctrlSet.A = []; % Ak
            ctrlSet.B = []; % Bk
            ctrlSet.C = []; % Ck
            ctrlSet.D = sdpvar(nu,ns,ctrlSet.nk,'full'); % Dk
        end
        
    case 'prA'
        ctrlSet.A = []; % Ak
        ctrlSet.B = sdpvar(ns,ny,ctrlSet.nk,'full'); % Bk
        ctrlSet.C = sdpvar(nu,ns,ctrlSet.nk,'full'); % Ck
        if (ctrlSet.dk)
            ctrlSet.D = sdpvar(nu,ny,ctrlSet.nk,'full'); % Dk
        else
            ctrlSet.D = zeros(nu,ny,ctrlSet.nk);
        end
        
    otherwise
        % full projection => no need for controller variables
        ctrlSet.A = []; % Ak
        ctrlSet.B = []; % Bk
        ctrlSet.C = []; % Ck
        ctrlSet.D = []; % Dk
end


% dealing with different options for control and measure channel parameter
% dependence:
if ismember(synSet.pdGclass,{'pass','ppss'}) && ~strcmp(lyapSet.type,'gral')
    % for the moment: only gral Lyapunov case is admitted in case of plants
    % type pass and ppss
    if ((synSet.pdMeasCtrl(2) > 0) && (ctrlSet.pdIn > 0))
       error('LPVSYN:inputError',...
             '%s\n\tfilter the output channel or set ctrlfcn.pdIn = false',ioMsg);
    elseif ((synSet.pdMeasCtrl(1) > 0) && (ctrlSet.pdOut > 0))
       error('LPVSYN:inputError',...
             '%s\n\tfilter the input channel or set ctrlfcn.pdOut = false',ioMsg);
    end
    
elseif strcmp(synSet.pdGclass,'pgss') && strcmp(lyapSet.type,'cte')
    % in this case the controller will have product of the original
    % functions in pdG
    if ((synSet.pdMeasCtrl(1) == 1) && (ctrlSet.pdOut == 1)) || ...
            ((synSet.pdMeasCtrl(2) == 1) && (ctrlSet.pdIn == 1)) || ...
            all(synSet.pdMeasCtrl)
        
        if all([synSet.pdMeasCtrl ctrlSet.pdIn ctrlSet.pdOut]) && ...
                (ctrlSet.dk > 0)
            error('LPVSYN:inputError',...
                'The controller will not be pgss, set ctrlfcn.dk = 0');
        end
    end
    
end
% if the plant is pgss and the X o Y parameter dependent, there is no
% limitations: the controller will be pcss


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for debugging
if options.debug
    synSet
    lyapSet    
    ctrlSet
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ========================================================================
% LMI constraints:
lmis = []; obj = 0;

% Gain - Hinf constraints
for ic = constSet.Gain.idx

    if isinf(constraint(ic).bound)
        constraint(ic).bound = sdpvar(1);
    end
    [lmiHinf,objHinf] = constHinf2LMI(pdG,synSet,lyapSet,ctrlSet,constraint(ic),options);
    lmis = [lmis, lmiHinf];
    obj = obj + objHinf;
    
end

% Gain - H2 constraints
for ic = constSet.GainH2.idx

    if isinf(constraint(ic).bound)
        constraint(ic).bound = sdpvar(1);
    end
    [lmiH2,objH2,Q] = constH22LMI(pdG,synSet,lyapSet,ctrlSet,constraint(ic),options);
    lmis = [lmis, lmiH2];
    obj = obj + objH2;
    
end

% Generalized Gain - H2 constraints
for ic = constSet.GainH2g.idx

    if isinf(constraint(ic).bound)
        constraint(ic).bound = sdpvar(1);
    end
    [lmiH2g,objH2g] = constH2g2LMI(pdG,synSet,lyapSet,ctrlSet,constraint(ic),options);
    lmis = [lmis, lmiH2g];
    obj = obj + objH2g;
    
end

% Poles constraints
for ic = constSet.Poles.idx

    lmiPoles = constPoles2LMI(pdG,synSet,lyapSet,ctrlSet,constraint(ic),options);
    lmis = [lmis, lmiPoles];
    
end

% extra LMI to impose a parameter controller dependence
if (ctrlSet.extFcn == 1) 
    fXEbnd = 0; 
    for ii = 1:size(pdG.parset.points,2)
        
        % controller vars at par
        [~,~,~,~,E,F] = evalCtrlVars(ctrlSet,ii,pdG);
        
        % Lyapunov matrices at par
        [Y,X] = evalLyapFcn(lyapSet,ii);
        
        na = size(E,2); 
        if (na > 0)
            % to limit the sv of X*E*E'*X & Y*F*F'*Y
            lmis = [lmis, [eye(ns) X*E;E'*X ctrlSet.XEbnd*eye(na)] >= 0];
            lmis = [lmis, [eye(ns) Y*F;F'*Y ctrlSet.YFbnd*eye(na)] >= 0];
            fXEbnd = 1;
        end
    end
    if (ctrlSet.XArYbnd > 0)
        penXArY = ctrlSet.XArYbnd;
    else
        penXArY = 1e-5;
    end
    % only add this obj if there is at least one E with na > 0
    if (fXEbnd == 1)
        obj  = obj + penXArY*(ctrlSet.XEbnd + ctrlSet.YFbnd);
    else
        ctrlSet.XEbnd = 0;
        ctrlSet.YEbnd = 0;
    end
        
end

% common to all constraints
[lmiPos,objPos] = constPos2LMI(pdG,synSet,lyapSet,ctrlSet,options);
lmis = [lmis, lmiPos];
obj = obj + objPos;

% bound on controller matrices
if (options.Ctrlpenalty > 0) && strcmp(synSet.type,'basic')
    
    eCtrl = sdpvar(1);
    if strcmp(ctrlSet.fback,'output')
        lmis = [lmis, -eCtrl <= ctrlSet.A(:) <= eCtrl];
        lmis = [lmis, -eCtrl <= ctrlSet.B(:) <= eCtrl];
        lmis = [lmis, -eCtrl <= ctrlSet.C(:) <= eCtrl];
    end
    if (ctrlSet.dk)
        lmis = [lmis, -eCtrl <= ctrlSet.D(:) <= eCtrl];
    end
    lmis = [lmis, eCtrl >= 0];
    obj = obj + options.Ctrlpenalty*eCtrl;
end


% -------------------------------------------------------------------------
% Solving optimization problem
if strcmp(options.solver, 'sedumi')
    solverOpts = sdpsettings('solver', 'sedumi',...
                  'verbose', options.verb,...
                  'dualize', options.dualize,...
                  'removeequalities', options.removeequalities);
%                   'sedumi.eps', options.solTol,...
%                   'sedumi.cg.qprec', 1,... 
%                   'sedumi.cg.maxiter', 49, ...
%                   'sedumi.stepdif', 2,...
%                   'sedumi.maxiter', options.maxiter,...
    
elseif strcmp(options.solver, 'sdpt3')
    solverOpts = sdpsettings('solver', 'sdpt3',...
                  'sdpt3.steptol', options.solTol,...
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
    
    fprintf('Infeasible synthesis problem\n')
    fprintf('\tSolver says: %s\n',yalmiperror(diagnostic.problem));
                
    pdK = []; 
    obj = Inf; 
    synProb = [];
    constSets = [constSet.Gain.idx constSet.GainH2.idx constSet.GainH2g.idx];
    for ii = constSets
        if isa(constraint(ii).bound,'sdpvar')
            constraint(ii).bound = inf;
        end
    end

else
    % Feasible
    
    if options.subOpt
        % in this case the LMIs is solved again but with higher bound
        % constraints to avoid the numerical issues of almost optimal
        % solutions
        
        for ic = constSet.Gain.idx
            if isa(constraint(ic).bound,'sdpvar')
                bnd = value(constraint(ic).bound)*(1+options.subOptBnd)^2;
                lmis = replace(lmis, constraint(ic).bound, bnd);
                obj = replace(obj, constraint(ic).bound, 0);
                constraint(ic).bound = sqrt(bnd);
            end
        end
        for ic = constSet.GainH2.idx
            if isa(constraint(ic).bound,'sdpvar')
                bnd = value(constraint(ic).bound)*(1+options.subOptBnd)^2;
                lmis = replace(lmis, constraint(ic).bound, bnd);
                obj = replace(obj, constraint(ic).bound, 0);
                constraint(ic).bound = sqrt(bnd);
            end
        end
        for ic = constSet.GainH2g.idx
            if isa(constraint(ic).bound,'sdpvar')
                bnd = value(constraint(ic).bound)*(1+options.subOptBnd)^2;
                lmis = replace(lmis, constraint(ic).bound, bnd);
                obj = replace(obj, constraint(ic).bound, 0);
                constraint(ic).bound = sqrt(bnd);
            end
        end
        
        diagnostic = optimize(lmis, obj, solverOpts);
        
    end

    if (diagnostic.problem == 4)
        fprintf('\n\n')
        warning backtrace off
        warning('Solver says numerical problems, the solution might not be reliable');
        warning backtrace on
        
    elseif ~((diagnostic.problem == 0) || (diagnostic.problem == 4))
        % should be feasible but ..
    
        fprintf('Infeasible synthesis problem\n')
        fprintf('\tSolver says: %s\n',yalmiperror(diagnostic.problem));
        
        pdK = [];
        obj = Inf;
        synProb = [];
        constSets = [constSet.Gain.idx constSet.GainH2.idx constSet.GainH2g.idx];
        for ii = constSets
            if isa(constraint(ii).bound,'sdpvar')
                constraint(ii).bound = inf;
            end
        end
        return
    end
    
    
    % Lyapunov functions
    lyapSet.X = value(lyapSet.X);
    lyapSet.Y = value(lyapSet.Y);
    % controller variables
    if ~isempty(ctrlSet.A)
        ctrlSet.A = value(ctrlSet.A);
    end
    if ~isempty(ctrlSet.B)
        ctrlSet.B = value(ctrlSet.B);
    end
    if ~isempty(ctrlSet.C)
        ctrlSet.C = value(ctrlSet.C);
    end
    if ~isempty(ctrlSet.D)
        ctrlSet.D = value(ctrlSet.D);
    end
    % for the implementation of gridding
    synProb.synSet = synSet;
    synProb.lyapSet = lyapSet;
    synProb.ctrlSet = ctrlSet;
    synProb.const = constraint;
    synProb.pdG = pdG;
    
    obj  = value(obj);

    % summarizing results
    n = length(constraint); jj = 1;
    for ii = 1:n
        if ~isempty(constraint(ii))
            [s1,s2,s3] = char(constraint(ii));
%         if ~isempty(s1)
            str1{jj} = s1; str2{jj} = s2; str3{jj} = s3;
            jj = jj + 1;
        end
    end
    idx_min = find(ismember(str1,'min'));
    idx_st = find(ismember(str1,'s.t.'));
    
    fprintf('\n------------------------------------------------------------------\n')
    if options.subOpt
        subStr = sprintf(' (Suboptimal by %1.2f%%)',options.subOptBnd);
    else
        subStr = '';
    end
    fprintf('Solved the following synthesis problem%s:\n\n',subStr)
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
            classConst = class(constraint(ii));
            switch classConst
                case {'synConst.Gain', 'synConst.GainH2', 'synConst.GainH2g'}
                    constraint(ii).bound = sqrt(value(constraint(ii).bound));
            end
            [~,str2,str3] = char(constraint(ii));
            fprintf('\t\t%s, \t(%s)\n',str2,str3);
        end
    end
    
    % controller re-construction
    % the output depends on the synthesis (ppss, pass, or pcss)
    pdK = lpvctrl(pdG,synSet,lyapSet,ctrlSet,constraint);
    
    fprintf('\n')
    fprintf('Controller: class: %s, scheme: %s, feedback: %s\n',class(pdK),...
        ctrlSet.sch,ctrlSet.fback);

    fprintf('\n')
    fprintf('Lyapunov functions\n');
    if strcmp(ctrlSet.fback,'output')
        if isa(lyapSet.parfcnX,'function_handle')
            fprintf('\tX:  %s\n',func2str(lyapSet.parfcnX));
        elseif ischar(lyapfcn) && any(strcmp(lyapfcn,{'aff','affX','affdX'}))
            fprintf('\tX:  %s\n','affine');
        elseif ischar(lyapfcn) && any(strcmp(lyapfcn,{'pwa','pwaX','pwadX'}))
            fprintf('\tX:  %s\n','PWA');
        else
            fprintf('\tX:  %s\n','constant');
        end        
        if isa(lyapSet.dparfcnX,'function_handle')
            fprintf('\tdX: %s\n',func2str(lyapSet.dparfcnX));
        elseif ischar(lyapfcn) && strcmp(lyapfcn,'affdX')
            fprintf('\tdX:  %s\n','affine');
        elseif ischar(lyapfcn) && strcmp(lyapfcn,'pwadX')
            fprintf('\tdX:  %s\n','PWA');
        else
            fprintf('\tdX: %s\n',num2str(lyapSet.dparfcnX));
        end        
    end
    if isa(lyapSet.parfcnY,'function_handle')
        fprintf('\tY:  %s\n',func2str(lyapSet.parfcnY));
    elseif ischar(lyapfcn) && any(strcmp(lyapfcn,{'aff','affY','affdY'}))
        fprintf('\tY:  %s\n','affine');
    elseif ischar(lyapfcn) && any(strcmp(lyapfcn,{'pwa','pwaY','pwadY'}))
        fprintf('\tY:  %s\n','PWA');
    else
        fprintf('\tY:  %s\n','constant');
    end
    if isa(lyapSet.dparfcnY,'function_handle')
        fprintf('\tdY: %s\n',func2str(lyapSet.dparfcnY));
    elseif ischar(lyapfcn) && strcmp(lyapfcn,'affdY')
        fprintf('\tdY:  %s\n','affine');
    elseif ischar(lyapfcn) && strcmp(lyapfcn,'pwadY')
        fprintf('\tdY:  %s\n','PWA');
    else
        fprintf('\tdY: %s\n',num2str(lyapSet.dparfcnY));
    end
    fprintf('------------------------------------------------------------------\n')
    
    % Closed loop poles
    fprintf('Closed-loop eigenvalues:\n')
    if (exist('ispsys.m','file') == 2) && isnumeric(pdK) && ispsys(pdK)
        typ = psinfo(pdK);
        if strcmp(typ,'aff')
            K = ss(pass(pdK));
        else
            K = ss(ppss(pdK));
        end
    else
        K = ss(pdK);
    end
    G = ss(pdG);
    for ii = 1:size(G,3)
        if strcmp(ctrlSet.fback,'output')
            Gcl = lft(G(:,:,ii),K(:,:,ii));
            Acl = Gcl.A;
        else
            [A,B2] = ssdata(G(:,'ctrl',ii));
            if strcmp(ctrlSet.sch,'gs')
                Ksf = K(:,:,ii).D;
            else
                Ksf = K.D;
            end
            Acl = A + B2*Ksf;
        end
        [wn,ep,poles] = damp(Acl);
        ReMax = max(real(poles));
        if (ReMax > 0)
            fprintf('\tMinDecay = %6.4f\t MinDamping = %6.4f \tmaxFreq = %6.4f << unstable!\n',...
                -ReMax,min(ep),max(wn));
        else
            fprintf('\tMinDecay = %6.4f\t MinDamping = %6.4f \tmaxFreq = %6.4f\n',...
                -ReMax,min(ep),max(wn));
        end
    end

    
    fprintf('------------------------------------------------------------------\n')
    
    % factorization
    if strcmp(ctrlSet.fback,'output')
        fprintf('I - XY:\n');
        if isa(pdK,'pcss')
            points = pdG.parset.points;
            for ii = 1:size(points,2)
                if isnumeric(pdK.xFcn)
                    X = pdK.xFcn;
                else
                    bnd = ss(pdK.xFcn, points(:,ii));
                    X = bnd.D;
                end
                if isnumeric(pdK.yFcn)
                    Y = pdK.yFcn;
                else
                    bnd = ss(pdK.yFcn, points(:,ii));
                    Y = bnd.D;
                end
                ePos = eig([X eye(ns); eye(ns) Y]); 
                eXY = eig(eye(ns) - X*Y);
                fprintf('\tMax eig(I-XY) %+d,',max(eXY));
                fprintf('\tMin eig([X I;I Y]) %+d\n',min(ePos));
                
            end
        else
            X = lyapSet.X;
            Y = lyapSet.Y;
            ePos = eig([X eye(ns); eye(ns) Y]);
            eXY = eig(eye(ns) - X*Y);
            fprintf('\tMax eig(I-XY) %+d,',max(eXY));
            fprintf('\tMin eig([X I;I Y]) %+d\n',min(ePos));
            
        end
        fprintf('------------------------------------------------------------------\n\n')
    end
    
    % info for the case of imposing ctrller dependence
    if (ctrlSet.extFcn == 1)
        G = ss(pdG); nv = size(G,3);
        fprintf('-------------------------------------------------------\n');
        fprintf('bound XE = %6.4f\t bound XF = %6.4f\n',value(ctrlSet.XEbnd),value(ctrlSet.YFbnd));
        for ii = 1:nv
            % auxiliary variables
            [~,~,~,~,E,F] = evalCtrlVars(ctrlSet,ii,pdG);
            % Lyapunov matrices at par
            [Y,X] = evalLyapFcn(lyapSet,ii);
            na = size(E);
            if (na > 0)
                nXE = norm(X*E)^2;
                nYF = norm(Y*F)^2;
                fprintf('ii: %3.0f  ||X*E||^2 = %6.4f \t ||Y*F||^2 = %6.4f\n',...
                    ii,nXE,nYF);
            end
        end
    end
    
end

% function outputs
if strcmp(synSet.pdGclass,'psys-aff') || strcmp(synSet.pdGclass,'psys-pol')
    varargout{1} = constraint(1).bound;
    varargout{2} = pdK;
else
    varargout{1} = pdK;
    varargout{2} = constraint;
    varargout{3} = obj;
    if (nargout > 3)
        varargout{4} = synProb;
    end
end


