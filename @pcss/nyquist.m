function varargout = nyquist(varargin)

% NYQUIST plot for LPV models at frozen parameter values 
%
% Use:
%   see NYQUIST
%
% NYQUIST creates the Nyquist plot of the LTI models obtained from 
% evaluating the LPV models at the points stored in the parameter set
%
% See also pass, ppss, pgss, pcss, nyquist

% fbianchi - 2021-07-02

% the pcss object is transformed in multidimensional ss obj
idxPCSS = find(cellfun(@(x)isa(x,'pcss'),varargin));
for ii = idxPCSS
    varargin{ii} = ss(varargin{ii});
end

% the result is then passed to the nyquist function
if (nargout == 0)
    nyquist(varargin{:});
else
    varargout = nyquist(varargin{:});
end