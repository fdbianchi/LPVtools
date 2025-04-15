function np = npar(obj)

% NPAR(pdG) returns the number of parameters in the LPV model pdG.
%
% See also pass, ppss, pgss, pcss

% fbianchi - 2021-03-31

np = size(obj.parset);


