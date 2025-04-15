function coeff_Sfunction(block)

% COEFF_SFUNCTION computes the coefficients to produce the parameter 
%   system matrices for pgss, pass and ppss
%
% Arguments:
%   - nu: number of time-varying parameters
%   - no: number of resulting coefficients
%   - fcn: function to compute the coeffienct for pgss
%   - type: type of interpolation pass, ppss or pgss
%
% fbianchi - 2023-12-12

setup(block);

%endfunction

% ===================================================================
% Function: setup 
%
function setup(block)

% Register the number of ports.
block.NumInputPorts  = 1;
block.NumOutputPorts = 1;

% Set up the port properties to be inherited or dynamic.
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% Register the parameters.
block.NumDialogPrms = 4;

% system information
type = block.DialogPrm(4).Data;
switch type
    case {'pass', 'pgss'}
        ny = [block.DialogPrm(2).Data 1];
    case 'ppss'
        ny = [block.DialogPrm(2).Data 1];
    case 'mat'
        ny = block.DialogPrm(2).Data;
end

% Input Port 1: system input
% block.InputPort(1).Dimensions  = nu;
block.InputPort(1).DatatypeID = 0;  % double
block.InputPort(1).Complexity = 'Real';
% block.InputPort(1).DirectFeedthrough = true;
 
% Output Port 1: system output
block.OutputPort(1).Dimensions = ny;
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
block.RegBlockMethod('Outputs', @Outputs);
%
block.RegBlockMethod('SetInputPortDimensions', @SetInpPortDims);


% ===================================================================
% The local functions

function SetInpPortDims(block, idx, di)
  
block.InputPort(idx).Dimensions = di;
  
% -------------------------------------------------------------------
% Compute outputs
function Outputs(block)

% system data
ny = block.DialogPrm(2).Data;
fcn = block.DialogPrm(3).Data;
type = block.DialogPrm(4).Data;

% parameter value with limits
p = block.InputPort(1).Data;  % parameter value

% matrix coefficients
switch type
    case 'pass'
        coeff = [1; p];
    case 'pgss'
        coeff = [1; fcn(p)];
    case 'ppss'
        if (ny == 1)
            coeff = 1;
        else
            coeff = cvxdec(fcn,p);
        end
    case 'mat'
        coeff = eye(ny(1));
end

block.OutputPort(1).Data = coeff;

