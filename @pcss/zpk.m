function sys = zpk(pdG,pval)

% ZPK converts a pcss object into a 3-D zpk object
%
% Use:
%   G = ZPK(pdG)
%   G = ZPK(pdG,pval)
%
% where 
%   - pdG:  LPV model
%   - pval: a matrix of (np x nv) with the parameter values in which the
%           LPV model is evaluated. 
%           (If pval is empty, the model is evaluated at the values stored 
%           in the parameter set)
%   - G:    zero-pole-gain model corresponding to the frozen parameter 
%           values (zpk object)
%
% ZPK(pdG) evaluates an LPV model pdG at frozen parameter values 
%   and produces a zero-pole-gain model G
%
% See also pass, ppss, pgss, pcss, zpk

% fbianchi - 2020-07-03

if (nargin > 1)
    % the system is evaluated at pval
    sys = zpk(ss(pdG,pval));

else
    % the system is evaluated at the points in the parameter set
    sys = zpk(ss(pdG));

end

