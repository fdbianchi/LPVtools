function bool = isempty(obj)

% ISEMPTY(lpvsys) returns true if the model has no inputs or no outputs and
% false otherwise  

% fbianchi - 2020-06-28

bool = isempty(obj.ctrller);
