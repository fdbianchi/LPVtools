function varargout = sigma(varargin)

% SIGMA singular value plot for LPV models at frozen parameter values 
%
% Use:
%   see SIGMA
%
% SIGMA plots the singular values of the LTI models obtained from 
% evaluating the LPV model at the points stored in the parameter set.
% 
% See also pass, ppss, pgss, pcss, sigma

% fbianchi - 2021-03-31

% the p_ss object is transformed in multidimensional ss obj
idxP_SS = find(cellfun(@(x)isa(x,'p_ss'),varargin));
for ii = idxP_SS
    varargin{ii} = ss(varargin{ii});
end

% the result is then passed to the sigma function
if (nargout == 0)
    sigma(varargin{:});
else
    varargout = sigma(varargin{:});
end

end
