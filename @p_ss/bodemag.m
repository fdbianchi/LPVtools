function varargout = bodemag(varargin)

% BODEMAG plot for LPV models at frozen parameter values
%
% Use:
%   see BODEMAG
%
% BODEMAG plots the magnitude of the frequency response of the LTI models
% obtained from evaluating the LPV models at the points stored in 
% in the parameter set.
%
% See also pass, ppss, pgss, pcss, bodemag

% fbianchi - 2021-03-31

% the p_ss object is transformed in multidimensional ss obj
idxP_SS = find(cellfun(@(x)isa(x,'p_ss'),varargin));
for ii = idxP_SS
    varargin{ii} = ss(varargin{ii});
end

% the result is then passed to the bodemag function
if (nargout == 0)
    bodemag(varargin{:});
else
    varargout = bodemag(varargin{:});
end

end
