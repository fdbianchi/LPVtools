function [no,ni] = iosize(obj)

% IOSIZE(pdG) returns the number of outputs and inputs of an LPV model pdG
%
% Use:
%   [ny,nu] = IOSIZE(pdG)
%   s = IOSIZE(pdG)
%
% where
%   - ny: number of outputs
%   - nu: number of inputs
%   - s: [ny,nu] 
%
% See also pass, ppss, pgss, pcss, iosize

% fbianchi - 2020-06-28

[no,ni] = size(obj.ctrller.D(:,:,1));

if (nargout < 2)
    no = [no,ni];
end


