function Go = minus(G1,G2)

% MINUS subtracts LPV or LTI models
%
% Use:
%   G3 = G2 - G1
%   G3 = MINUS(G1,G2)
%
% where:
%   G1 or G2 are pass, ppss, pgss, ss, tf, or zpk objects.
%
% See also pass, ppss, pgss, pcss, minus

% fbianchi - 2021-03-31

Go = parallel(G1,-G2);
