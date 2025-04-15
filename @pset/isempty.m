function bool = isempty(obj)

% ISEMPTY checks whether the pset is empty or not
%
% See also pset.Box, pset.Grid, pset.Hull, pset.Gral

% fbianchi - 2021-03-31

bool = isempty(obj.points);


