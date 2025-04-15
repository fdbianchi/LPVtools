function varargout = bode(varargin)

% BODE plot for LPV models at frozen parameter values
%
% Use:
%   see BODE
%
% BODE plots the magnitude and phase of the frequency response of the LTI 
% models obtained from evaluating the LPV models at the points stored in 
% in the parameter set.
%
% See also pass, ppss, pgss, pcss, bode

% fbianchi - 2021-03-31


% the p_ss object is transformed in multidimensional ss obj
idxP_SS = find(cellfun(@(x)isa(x,'p_ss'),varargin));
for ii = idxP_SS
    varargin{ii} = ss(varargin{ii});
end

% the result is then passed to the bode function
if (nargout == 0)
    bode(varargin{:});
else
    varargout = bode(varargin{:});
end

