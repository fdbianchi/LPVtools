function varargout = bodemag(varargin)

% BODEMAG plot for LPV models at frozen parameter values
%
% Use:
%   see BODEMAG
%
% BODEMAG plots the magnitude of the frequency response of the LTI models
% obtained from evaluating the LPV model at the points stored in 
% the parameter set
%
% See also pass, ppss, pgss, pcss, bodemag

% fbianchi - 2020-07-03

% the pcss object is transformed into multidimensional ss obj
idxPCSS = find(cellfun(@(x)isa(x,'pcss'),varargin));
for ii = idxPCSS
    varargin{ii} = ss(varargin{ii});
end

% the result is then passed to the bodemag function
if (nargout == 0)
    bodemag(varargin{:});
else
    varargout = bodemag(varargin{:});
end

end
