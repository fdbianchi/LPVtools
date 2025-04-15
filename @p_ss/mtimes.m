function Go = mtimes(G2,G1)

% MTIMES connects LPV or LTI models in series
%
% Use:
%   G3 = G2*G1
%   G3 = MTIMES(G2,G1)
%
% where:
%   G1 or G2 are pass, ppss, pgss, ss, tf, or zpk objects.
%
% See also pass, ppss, pgss, pcss, mtimes

% fbianchi - 2021-03-31

Go = series(G1,G2);
