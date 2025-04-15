function Go = plus(G1,G2)

% PLUS connects LPV or LTI models in parallel
%
% Use:
%   G3 = G2 + G1
%   G3 = PLUS(G2,G1)
%
% where:
%   G1 or G2 are pass, ppss, pgss, ss, tf, or zpk objects.
%
% See also pass, ppss, pgss, pcss, plus

% fbianchi - 2021-03-31

Go = parallel(G1,G2);
