function pv = pvec(obj)

% PVEC converts a pset.Gral object into a pvec (lmitool) object
%
% Use: 
%   pv = pvec(set)
%
% Inputs: 
%   set: pset object
%
% Output:
%   pv: pvec object
%
% See also PSET, PVEC

% fbianchi - 2021-03-31

pv = pvec('pol',obj.points);
