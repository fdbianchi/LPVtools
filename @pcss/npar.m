function np = npar(obj)

% NPAR(pdG) returns the number of parameters in the LPV model pdG.

% fbianchi - 2020-02-20

np = size(obj.ctrller.parset);

