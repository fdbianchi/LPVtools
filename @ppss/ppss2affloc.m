function pdGloc = ppss2affloc(pdG, idx)

% PPSS2AFFLOC produces an affine local model from a PWA model
%
% Use:
%   pdGloc = ppss2affloc(pdG, idx)
%
% where:
%   pdG: ppss object with the PWA model
%   idx: index corresponding to the simplex
%   pdGloc: pass model valid for the idx simplex
%
% pdGloc is an affine lpv model valid in the simplex IDX

% fbianchi - 2023-08-10

if ~isa(pdG, 'ppss')
    error('PPSS2AFFLOC:inputError','Invalid pdG class')
end

if ~(isnumeric(idx) && isscalar(idx))
    error('PPSS2AFFLOC:inputError','idx must be a numeric scalar')
end
    

% tolerance to remove very small values
ntol = 1e-6;

% plant dimensions
[no,ni,ns] = size(pdG);
% LTI models
Gs = ss(pdG);

% Parameter set info
Pset = pdG.parset;
simplices = Pset.simplices;
if (idx > size(simplices,1))
    error('PPSS2AFFLOC:inputError','idx exceeds the number of simplices')
else
    simplex = Pset.simplices(idx, :);
    nv = length(simplex);
end    

% center point in the simplex
p0 = Pset.points(:,simplex)*ones(nv,1)/nv;

% transformation from ppss to pass
Delta = [ones(1,nv);
         Pset.points(:,simplex) - p0*ones(1,nv)];

% matrix A     
aux = reshape(Gs.A(:,:,simplex),ns*ns,nv);
aux = aux/Delta;
aux(abs(aux) < max(ntol,ntol*max(max(aux)))) = 0;
Aaff = reshape(aux,ns,ns,nv);

% matrix B
aux = reshape(Gs.B(:,:,simplex),ns*ni,nv);
aux = aux/Delta;
aux(abs(aux) < max(ntol,ntol*max(max(aux)))) = 0;
Baff = reshape(aux,ns,ni,nv);

% matrix C
aux = reshape(Gs.C(:,:,simplex),no*ns,nv);
aux = aux/Delta;
aux(abs(aux) < max(ntol,ntol*max(max(aux)))) = 0;
Caff = reshape(aux,no,ns,nv);

% matrix D
aux = reshape(Gs.D(:,:,simplex),no*ni,nv);
aux = aux/Delta;
aux(abs(aux) < max(ntol,ntol*max(max(aux)))) = 0;
Daff = reshape(aux,no,ni,nv);

% new local parameter set
PsetLoc = pset.Gral(Pset.points(:,simplex) - repmat(p0,1,nv));

% new local affine model
pdGloc = pass(Aaff, Baff, Caff, Daff,...
            PsetLoc,...
            'StateName', pdG.StateName,...
            'InputName', pdG.InputName,...
            'OutputName', pdG.OutputName,...
            'InputGroup', pdG.InputGroup,...
            'OutputGroup', pdG.OutputGroup);





