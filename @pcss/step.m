function varargout = step(varargin)

% STEP response plot for LPV models at frozen parameter values  
%
% Use:
%   see STEP
%
% STEP plots the step reponse for LPV models at the points stored in
% the parameter set
%
% See also pass, ppss, pgss, pcss, step

% fbianchi - 2020-07-03

% the pcss object is transformed in multidimensional ss obj
idxPCSS = find(cellfun(@(x)isa(x,'pcss'),varargin));
for ii = idxPCSS
    varargin{ii} = ss(varargin{ii});
end

% the result is then passed to the step function
if (nargout == 0)
    step(varargin{:});
else
    varargout = step(varargin{:});
end

