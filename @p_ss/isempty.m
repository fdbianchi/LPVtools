function bool = isempty(obj)

% ISEMPTY(pdG) returns true if the LPV model pdG has no inputs or no outputs 
%   and false otherwise  
%
% See also pass, ppss, pgss, pcss

% fbianchi - 2021-03-31

bool = isempty(obj.D(:,:,1));
