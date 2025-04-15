function G = ss(obj,pval)

% SS(pdG) evaluates an LPV model pdG at frozen parameter values 
%   and produces a state-space model G
%
% Use:
%   G = SS(pdG)
%   G = SS(pdG,pval)
%
% where 
%   - pdG:  LPV model
%   - pval: a matrix of (np x nv) with the parameter values at which the
%           LPV model is evaluated. 
%           (If pval is empty, the model is evaluated at the values stored 
%           in the parameter set)
%   - G:    state-space model corresponding to the frozen parameter values
%           (ss object)
%
% See also pass, pgss, ppss, pcss, ss, tf, zpk, subs

% fbianchi - 2021-03-31


if (nargin < 2)
    % all points in parset
    pval = obj.parset.points;
end
    
% the p_ss object is transformed into multidimensional ss obj
[G,msg] = subs(obj,pval);
if isempty(G)
    error(msg)
end


