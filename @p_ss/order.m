function ns = order(obj)

% ORDER(pdG) returns the order of the LPV model pdG
%
% See also pass, ppss, pgss, pcss

% fbianchi - 2021-03-31

ns = size(obj.A(:,:,1),1);

       
