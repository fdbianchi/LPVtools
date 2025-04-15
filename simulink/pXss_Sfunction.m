function pXss_Sfunction(block)

% PXSS_SFUNCTION implements affine, general or polytopic LPV models
%
% Input parameter:
%   - sys:  object pass or ppss
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
if (all(all(sys.D)) == 0)
    block.InputPort(1).DirectFeedthrough = false;
else
    block.InputPort(1).DirectFeedthrough = true;
end    

% Input Port 2: time-varying parameter
block.InputPort(2).Dimensions  = np;
block.InputPort(2).DatatypeID  = 0;  % double
block.InputPort(2).Complexity  = 'Real';
block.InputPort(2).DirectFeedthrough = true;

% Output Port 3: system output
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
block.RegBlockMethod('Outputs',              @Outputs);
if block.NumContStates > 0
    block.RegBlockMethod('Derivatives',          @Derivatives);
    block.RegBlockMethod('InitializeConditions', @InitializeConditions);
end


% ===================================================================
% The local functions

% -------------------------------------------------------------------
% Checking parameters
% -------------------------------------------------------------------
function CheckPrms(block)
    
% checking system
sys = block.DialogPrm(1).Data;
if ~isa(sys,'p_ss')
    error('System must be pass or ppss object')
end

% checking initial conditions
x0 = block.DialogPrm(2).Data;
ns = order(sys);
if isscalar(x0) && (x0 ~= 0)
    error('x0 must be a vector of %d elements or 0',ns)
elseif ~isvector(x0) && (length(x0) ~= ns)
    error('x0 must be a vector of %d elements or 0',ns)
end
  
% -------------------------------------------------------------------
% Set initial conditions
% -------------------------------------------------------------------
function InitializeConditions(block)

if block.NumContStates > 0
    x0 = block.DialogPrm(2).Data;
    if isscalar(x0) && (x0 == 0)
        sys = block.DialogPrm(1).Data;
        ns = order(sys);
        block.ContStates.Data = zeros(ns,1);
    else
        block.ContStates.Data = x0;
    end
end

% -------------------------------------------------------------------
% Compute outputs
function Outputs(block)

% system data
sys = block.DialogPrm(1).Data;

% parameter value with limits
p    = block.InputPort(2).Data;  % parameter value
pmin = sys.parset.range(:,1);
pmax = sys.parset.range(:,2);
p    = min(max(p,pmin),pmax);

% matrix coefficients
if isa(sys,'pass')
    coeff = [1; p];
    idx   = 1:length(coeff);
elseif isa(sys,'pgss')
    coeff = [1; sys.parfcn(p)];
    idx   = 1:length(coeff);
elseif isa(sys,'ppss')
    [coeff, idx] = cvxdec(sys.parset,p);
end

% system matrices
nv = length(coeff); 
[no,ni] = size(sys.D(:,:,1)); ns = size(sys.A(:,:,1),1);
C = reshape(reshape(sys.C(:,:,idx),no*ns,nv,1)*coeff,no,ns,1);
D = reshape(reshape(sys.D(:,:,idx),ni*no,nv,1)*coeff,no,ni,1);
        
if block.NumContStates > 0
    block.OutputPort(1).Data = C*block.ContStates.Data + D*block.InputPort(1).Data;
else
    block.OutputPort(1).Data = D*block.InputPort(1).Data;
end    

% -------------------------------------------------------------------
% Compute derivative
function Derivatives(block)

% system data
sys = block.DialogPrm(1).Data;

% parameter value with limits
p    = block.InputPort(2).Data;  % parameter value
pmin = sys.parset.range(:,1);
pmax = sys.parset.range(:,2);
p    = min(max(p,pmin),pmax);

% matrix coefficients
if isa(sys,'pass')
    coeff = [1; p];
    idx   = 1:length(coeff);
elseif isa(sys,'pgss')
    coeff = [1; sys.parfcn(p)];
    idx   = 1:length(coeff);
elseif isa(sys,'ppss')
    [coeff, idx] = cvxdec(sys.parset,p);
end

% system matrices
nv = length(coeff);
[no,ni] = size(sys.D(:,:,1)); ns = size(sys.A(:,:,1),1);
A = reshape(reshape(sys.A(:,:,idx),ns*ns,nv,1)*coeff,ns,ns,1);
B = reshape(reshape(sys.B(:,:,idx),ns*ni,nv,1)*coeff,ns,ni,1);

block.Derivatives.Data = A*block.ContStates.Data + B*block.InputPort(1).Data;

