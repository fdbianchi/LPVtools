function n = nsys(obj)

% NSYS(pdG) returns the number of LTI models used in the representation 
% corresponding to the LPV model pdG.
%
% See also pass, ppss, pgss, pcss

% fbianchi - 2021-03-31

n = size(obj.D,3);

        
