function pcssMatCoeff_Sfunction(block)

% PCSS_MATCOEFF_SFUNCTION undoes the controller change of variables given
% the auxiliary matrices, the Lyapunov matrices and the plant matrices
%
% Input parameter:
%   - ns: number of states
%   - ny: number of inputs
%   - nu: number of outputs
%   - flagX: 1 is X is numeric and 0 otherwise
%   - flagY: 1 is Y is numeric and 0 otherwise

% fbianchi - 2020-07-02

setup(block);

%endfunction

% ===================================================================
% Function: setup 
%
function setup(block)

% Register the number of ports.
block.NumInputPorts  = 4;
block.NumOutputPorts = 1;

% Set up the port properties to be inherited or dynamic.
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% Register the parameters.
block.NumDialogPrms = 5;

% system information
ns = block.DialogPrm(1).Data;
ny = block.DialogPrm(2).Data;
nu = block.DialogPrm(3).Data;

% Input Port 1: system input
block.InputPort(1).Dimensions = [ns+ny ns+nu];
block.InputPort(1).DatatypeID = 0;  % double
block.InputPort(1).Complexity = 'Real';
block.InputPort(1).DirectFeedthrough = true;

% Input Port 2: system input
block.InputPort(2).Dimensions = [ns ns];
block.InputPort(2).DatatypeID = 0;  % double
block.InputPort(2).Complexity = 'Real';
block.InputPort(2).DirectFeedthrough = true;

% Input Port 3: system input
block.InputPort(3).Dimensions = [ns ns];
block.InputPort(3).DatatypeID = 0;  % double
block.InputPort(3).Complexity = 'Real';
block.InputPort(3).DirectFeedthrough = true;

% Input Port 4: time-varying parameter
block.InputPort(4).Dimensions = [ns+nu ns+ny];
block.InputPort(4).DatatypeID = 0;  % double
block.InputPort(4).Complexity = 'Real';
block.InputPort(4).DirectFeedthrough = true;

% Output Port 1: system output
block.OutputPort(1).Dimensions = [ns+ny ns+nu];
block.OutputPort(1).DatatypeID = 0; % double
block.OutputPort(1).Complexity = 'Real';

% Register the sample times.
block.SampleTimes = [-1 0];

% -----------------------------------------------------------------
% Options
% -----------------------------------------------------------------
block.SetAccelRunOnTLC(false);
block.SimStateCompliance = 'DefaultSimState';

% -----------------------------------------------------------------
% Register the methods called during update diagram/compilation.
% -----------------------------------------------------------------
% block.RegBlockMethod('CheckParameters',      @CheckPrms);
% block.RegBlockMethod('InitializeConditions', @InitializeConditions);
block.RegBlockMethod('Outputs',              @Outputs);
% block.RegBlockMethod('Derivatives',          @Derivatives);
% block.RegBlockMethod('PostPropagationSetup', @DoPostPropSetup);


% ===================================================================
% The local functions

% -------------------------------------------------------------------
% Compute outputs
function Outputs(block)

% system dimensions
ns = block.DialogPrm(1).Data;
ny = block.DialogPrm(2).Data;
nu = block.DialogPrm(3).Data;

flagX = block.DialogPrm(4).Data;
flagY = block.DialogPrm(5).Data;

% auxiliary controller matrices
K = block.InputPort(1).Data;

% Lyapunov matrices
X = block.InputPort(2).Data;
Y = block.InputPort(3).Data;

% plant matrices
P = block.InputPort(4).Data;
A = P(1:ns,1:ns);
b2 = P(1:ns,1+ns:end);
c2 = P(1+ns:end,1:ns);

% Lyapunov function at par
if (flagX == 0) && (flagY == 0)
    % X & Y parameter dependent
    Ins = eye(ns);
    Iny = eye(ny);
    Inu = eye(nu);
    
    % The ambiguity in the factorization produces discontinuities in the control
    % using QR => 
    % [Q,R] = qr(eye(ns) - X*Y);
    % Ni = Q';
    % Mti = inv(R);
     
    % using SVD
    % [U,S,V] = svd(eye(ns) - X*Y);
    % S  = S^(-0.5);
    % Ni = S*U';
    % Mti = V*S;
  
    Ni = Ins/X;
    Mti = (Ni - Y)\Ins;

    % undoing variable change
    Kf = blkdiag(Ni,eye(ny))*...
        ([Ins -X*b2; zeros(ny,ns) Iny]*K*[Ins zeros(ns,nu); -c2*Y Inu] - ...
        blkdiag(X*A*Y,zeros(ny,nu)))*blkdiag(Mti,Inu);
   
elseif (flagX == 1) && (flagY == 0)
    % X constant & Y parameter dependent
    Ins = eye(ns);
    Iny = eye(ny);
    Inu = eye(nu);
    % N  = eye(ns);
    Mti = Ins/(Ins - X*Y);

    % undoing variable change
    Kf = ([Ins -X*b2; zeros(ny,ns) Iny]*K*[Ins zeros(ns,nu); -c2*Y Inu] - ...
          blkdiag(X*A*Y,zeros(ny,nu)))*blkdiag(Mti,Inu);
    
elseif (flagX == 0) && (flagY == 1)
    % X parameter dependent & X constant
    Ins = eye(ns);
    Iny = eye(ny);
    Inu = eye(nu);
    Ni  = (Ins - X*Y)\Ins;
    % Mt = eye(ns);

    % undoing variable change
    Kf = blkdiag(Ni,Iny)*...
        ([Ins -X*b2; zeros(ny,ns) Iny]*K*[Ins zeros(ns,nu); -c2*Y Inu] - ...
        blkdiag(X*A*Y,zeros(ny,nu)));
    
end

block.OutputPort(1).Data = Kf;


