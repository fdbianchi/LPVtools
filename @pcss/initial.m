function varargout = initial(varargin)

% INITIAL plot for LPV models at frozen parameter values
%
% Use:
%   see INITIAL
%
% INITIAL plots the unforced response of the LTI obtained from evaluating
% the lpv models at the grid points stored in the parameter set.

% fbianchi - 2020-07-03

% the pcss object is transformed in multidimensional ss obj
idxPCSS = find(cellfun(@(x)isa(x,'pcss'),varargin));
for ii = idxPCSS
    varargin{ii} = ss(varargin{ii});
end

% the result is then passed to the initial function
if (nargout == 0)
    initial(varargin{:});
else
    varargout = initial(varargin{:});
end

