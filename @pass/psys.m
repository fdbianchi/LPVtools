function pdGo = psys(obj)

% PSYS converts a pass object into a psys object (lmitool)
%
% Use:
%   pdGo = PSYS(pdGi)
%
% where:
%   - pdGi: pass object
%   - pdGo: psys object
%
% See also pass, psys

% fbianchi - 2021-03-31


% copy system matrices
pdGo = ltisys(obj.A(:,:,1),obj.B(:,:,1),obj.C(:,:,1),obj.D(:,:,1));
for ii = 2:nsys(obj)
    pdGo = [pdGo, ltisys(obj.A(:,:,ii),obj.B(:,:,ii),obj.C(:,:,ii),...
        obj.D(:,:,ii),0)];
end

% copy the parameter set
if isa(obj.parset,'pset.Hull') || isa(obj.parset,'pset.Grid') || isa(obj.parset,'pset.Gral')
    pv = pvec(obj.parset.points);
else
    pv   = pvec(obj.parset);
end

% new representation
pdGo = psys(pv,pdGo);


