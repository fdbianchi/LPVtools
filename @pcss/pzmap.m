function varargout = pzmap(varargin)

% PZMAP plot for LPV models at frozen parameter values
%
% Use:
%   see PZMAP
%
% PZMAP plots pole-zero map of the LTI models obtained from evaluating 
% the LPV model at the points stored in the parameter set
%
% See also pass, ppss, pgss, pcss, pzmap

% fbianchi - 2021-03-31

% the pcss object is transformed in multidimensional ss obj
idxPCSS = find(cellfun(@(x)isa(x,'pcss'),varargin));
for ii = idxPCSS
    varargin{ii} = ss(varargin{ii});
end

% the result is then passed to the pzmap function
if (nargout == 0)
    pzplot(varargin{:});
else
    varargout = pzplot(varargin{:});
end

