function pv = pvec(obj)

% PVEC converts a pset.Box object into a pvec (lmitool) object
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
% See also PSET.BOX, PVEC

% fbianchi - 2021-03-29

if isempty(obj.rate)
    pv = pvec('box',obj.range);
else
    pv = pvec('box',obj.range,obj.rate);
end

