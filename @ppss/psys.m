function pdGo = psys(obj)

% PSYS converts a ppss object into a psys object (lmitool)
%
% Use:
%   pdGo = PSYS(pdGi)
%
% where:
%   - pdGi: ppss object
%   - pdGo: psys object
%
% See also ppss, psys

% fbianchi - 2021-03-31

% copy system matrices
pdGo = [];
for ii = 1:nsys(obj)
    pdGo = [pdGo, ltisys(obj.A(:,:,ii),obj.B(:,:,ii),obj.C(:,:,ii),obj.D(:,:,ii))];
end
pdGo = psys(pdGo);

% copy the parameter set
pv   = pvec(obj.parset);
pdGo = addpv(pdGo,pv);


end