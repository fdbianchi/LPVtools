function pdG = denormalize(pdGn, Pset)

% DENORMALIZE transforms the model pdGn with parameters taking value 
% between -1 and 1 into one with a new parameter set Pset
%
% Use:
%   pdG = DENORMALIZE(pdGn, Pset)
%
% See also normalize

% fbianchi - 2024-01-16

if ~isa(Pset,'pset')
    error('PASS:DENORMALIZE:inputError','PSET must be a pset object')
end

% model dimensions
[~,~,~,np,nv] = size(pdGn);

% new parameter set
p0 = mean(Pset.points,2);
pd = max(Pset.points,[],2) - p0;

% matrix modification
A(:,:,1) = pdGn.A(:,:,1);
B(:,:,1) = pdGn.B(:,:,1);
C(:,:,1) = pdGn.C(:,:,1);
D(:,:,1) = pdGn.D(:,:,1);
for ii = 2:nv
    fact0 = p0(ii-1)/pd(ii-1);
    A(:,:,1) = A(:,:,1) - fact0*pdGn.A(:,:,ii);
    B(:,:,1) = B(:,:,1) - fact0*pdGn.B(:,:,ii);
    C(:,:,1) = C(:,:,1) - fact0*pdGn.C(:,:,ii);
    D(:,:,1) = D(:,:,1) - fact0*pdGn.D(:,:,ii);
    
    factd = 1/pd(ii-1);
    A(:,:,ii) = factd*pdGn.A(:,:,ii);
    B(:,:,ii) = factd*pdGn.B(:,:,ii);
    C(:,:,ii) = factd*pdGn.C(:,:,ii);
    D(:,:,ii) = factd*pdGn.D(:,:,ii);

end


pdG = pass(A, B, C, D, Pset,...
            'StateName',pdGn.StateName,...
            'InputName',pdGn.InputName,...
            'OutputName',pdGn.OutputName,...
            'InputGroup',pdGn.InputGroup,...
            'OutputGroup',pdGn.OutputGroup);
