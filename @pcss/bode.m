function varargout = bode(varargin)

% BODE plot for LPV models at frozen parameter values
%
% Use:
%   see BODE
%
% BODE plots the magnitude and phase of the frequency response of the LTI 
% models obtained from evaluating the LPV model at the points stored in 
% the parameter set
%
% See also pass, ppss, pgss, pcss, bode

% fbianchi - 2020-07-02

% the pcss object is transformed in multidimensional ss obj
idxPCSS = find(cellfun(@(x)isa(x,'pcss'),varargin));
for ii = idxPCSS
    varargin{ii} = ss(varargin{ii});
end

% the result is then passed to the bode function
if (nargout == 0)
    bode(varargin{:});
else
    varargout = bode(varargin{:});
end

end

