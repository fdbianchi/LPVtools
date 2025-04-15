function varargout = step(varargin)

% STEP response plot for LPV models at frozen parameter values  
%
% Use:
%   see STEP
%
% STEP plots the step reponse for the LPV model at the points stored in
% the parameter set.
%
% See also pass, ppss, pgss, pcss, step

% fbianchi - 2021-03-31

% the p_ss object is transformed in multidimensional ss obj
idxP_SS = find(cellfun(@(x)isa(x,'p_ss'),varargin));
for ii = idxP_SS
    varargin{ii} = ss(varargin{ii});
end

% the result is then passed to the step function
if (nargout == 0)
    step(varargin{:});
else
    varargout = step(varargin{:});
end

