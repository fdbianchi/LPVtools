function varargout = impulse(varargin)

% IMPULSE plot for the LPV model at frozen parameter values
%
% Use:
%   see IMPULSE
%
% IMPULSE plots the impulse reponse of the LTI obtained from evaluating
% the lpv models evaluated at vertices the parameter set or at the grid 
% points stored in the parameter set

% fbianchi - 2020-07-02

% the pcss object is transformed in multidimensional ss obj
idxPCSS = find(cellfun(@(x)isa(x,'pcss'),varargin));
for ii = idxPCSS
    varargin{ii} = ss(varargin{ii});
end

% the result is then passed to the impulse function
if (nargout == 0)
    impulse(varargin{:});
else
    varargout = impulse(varargin{:});
end

end
