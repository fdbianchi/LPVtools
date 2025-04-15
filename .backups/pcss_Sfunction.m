function pcss_Sfunction(block)

% PCSS_SFUNCTION provides the on-line implementation of LPV controllers
% obtained from synthesis with parameter dependent Lyapunov functions
%
% Input parameter:
%   - sys:  sysect pcss
%   - x0:   0 or vector with the initial conditions

% fbianchi - 2020-07-02

setup(block);

%endfunction

% ===================================================================
% Function: setup 
%
function setup(block)

% Register the number of ports.
block.NumInputPorts  = 2;
block.NumOutputPorts = 1;

% Set up the port properties to be inherited or dynamic.
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% Register the parameters.
block.NumDialogPrms     = 2;

% system information
sys = block.DialogPrm(1).Data;
[ny,nu,ns,np] = size(sys);

% Input Port 1: system input
block.InputPort(1).Dimensions  = nu;
block.InputPort(1).DatatypeID  = 0;  % double
block.InputPort(1).Complexity  = 'Real';
if (all(all(sys.ctrller.D)) == 0)
    block.InputPort(1).DirectFeedthrough = false;
else
    block.InputPort(1).DirectFeedthrough = true;
end    

% Input Port 2: time-varying parameter
block.InputPort(2).Dimensions  = np;
block.InputPort(2).DatatypeID  = 0;  % double
block.InputPort(2).Complexity  = 'Real';
block.InputPort(2).DirectFeedthrough = true;

% Output Port 1: system output
block.OutputPort(1).Dimensions  = ny;
block.OutputPort(1).DatatypeID  = 0; % double
block.OutputPort(1).Complexity  = 'Real';

% Set up the continuous states.
block.NumContStates = ns;

% Register the sample times.
block.SampleTimes = [0 0];

% -----------------------------------------------------------------
% Options
% -----------------------------------------------------------------
block.SetAccelRunOnTLC(false);
block.SimStateCompliance = 'DefaultSimState';

% -----------------------------------------------------------------
% Register the methods called during update diagram/compilation.
% -----------------------------------------------------------------
block.RegBlockMethod('CheckParameters',      @CheckPrms);
block.RegBlockMethod('InitializeConditions', @InitializeConditions);
block.RegBlockMethod('Outputs',              @Outputs);
block.RegBlockMethod('Derivatives',          @Derivatives);
block.RegBlockMethod('PostPropagationSetup', @DoPostPropSetup);


% ===================================================================
% The local functions

% -------------------------------------------------------------------
% Checking parameters
% -------------------------------------------------------------------
function CheckPrms(block)
    
% checking system
sys = block.DialogPrm(1).Data;
if ~isa(sys,'pcss')
    error('System must be a pcss model')
end

% checking initial conditions
x0 = block.DialogPrm(2).Data;
ns = order(sys);
if isscalar(x0) && (x0 ~= 0)
    error('x0 must be a vector of %d elements or 0',ns)
elseif ~isvector(x0) && (length(x0) ~= ns)
    error('x0 must be a vector of %d elements or 0',ns)
end
  
function DoPostPropSetup(block)

  block.NumDworks = 2;
  
  [ns,nu] = size(block.DialogPrm(1).Data.ctrller.B(:,:,1));
   
  block.Dwork(1).Name            = 'ak';
  block.Dwork(1).Dimensions      = ns*ns;
  block.Dwork(1).DatatypeID      = 0;      % double
  block.Dwork(1).Complexity      = 'Real'; % real
  block.Dwork(1).UsedAsDiscState = false;

  block.Dwork(2).Name            = 'bk';
  block.Dwork(2).Dimensions      = ns*nu;
  block.Dwork(2).DatatypeID      = 0;      % double
  block.Dwork(2).Complexity      = 'Real'; % real
  block.Dwork(2).UsedAsDiscState = false;

  %endfunction

% -------------------------------------------------------------------
% Set initial conditions
% -------------------------------------------------------------------
function InitializeConditions(block)

    sys = block.DialogPrm(1).Data;
    ns = order(sys);

    x0 = block.DialogPrm(2).Data;
    if isscalar(x0) && (x0 == 0)
        block.ContStates.Data = zeros(ns,1);
    else
        block.ContStates.Data = x0;
    end

% -------------------------------------------------------------------
% Compute outputs
function Outputs(block)

% system data
sys = block.DialogPrm(1).Data;

% parameter value with limits
p    = block.InputPort(2).Data;  % parameter value
pmin = sys.ctrller.parset.range(:,1);
pmax = sys.ctrller.parset.range(:,2);
p    = min(max(p,pmin),pmax);

% model dimensions
ns = order(sys);

% Lyapunov function at par
flagX = isnumeric(sys.xFcn);
flagY = isnumeric(sys.yFcn);
if (flagX == 0) && (flagY == 0)
    Y = subs(sys.yFcn,p); Y = Y.D;
    X = subs(sys.xFcn,p); X = X.D;
    [U,S,V] = svd(eye(ns) - X*Y);
    S  = S^(0.5);
    N  = U*S;
    Mt = S*V';
    
elseif (flagX == 1) && (flagY == 0)
    X = sys.xFcn;
    Y = subs(sys.yFcn,p); Y = Y.D;
    N  = eye(ns);
    Mt = (eye(ns) - X*Y);
    
elseif (flagX == 0) && (flagY == 1)
    X = subs(sys.xFcn,p); X = X.D;
    Y = sys.yFcn;
    N  = (eye(ns) - X*Y);
    Mt = eye(ns);
    
end
               
% controller variables at p
Khat = subs(sys.ctrller,p);
Ak = Khat.A;
Bk = Khat.B;
Ck = Khat.C;
dk = Khat.D;

% augmented plant matrices
plant = subs(sys.plant,p);
A  = plant.A;
b2 = plant.B;
c2 = plant.C; 

% undoing variable change
ak = N\(Ak - X*(A-b2*dk*c2)*Y - Bk*c2*Y - X*b2*Ck)/Mt;
bk = N\(Bk - X*b2*dk);
ck = (Ck - dk*c2*Y)/Mt;

% save matrices to compute derivative
block.Dwork(1).Data = ak(:);
block.Dwork(2).Data = bk(:);

block.OutputPort(1).Data = ck*block.ContStates.Data + dk*block.InputPort(1).Data;


% -------------------------------------------------------------------
% Compute derivative
function Derivatives(block)

[ns,nu] = size(block.DialogPrm(1).Data.ctrller.B(:,:,1));

ak = reshape(block.Dwork(1).Data,ns,ns);
bk = reshape(block.Dwork(2).Data,ns,nu);

block.Derivatives.Data = ak*block.ContStates.Data + bk*block.InputPort(1).Data;

