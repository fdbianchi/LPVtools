function pdGn = normalize(pdG)

% NORMALIZE transforms the model pdG into a model with parameters taking
% value between -1 and 1
%
% Use:
%   pdGn = NORMALIZE(pdG)
%
% See also denormalize

% fbianchi - 2024-01-16


% model dimensions
[~,~,~,np,nv] = size(pdG);

% normalization of the parameter set
parset = pdG.parset;
p0 = mean(parset.points,2);
pd = max(parset.points,[],2) - p0;
if isa(parset, 'pset.Box')
    newParset = pset.Box(repmat([-1 1],np,1),...
                         parset.rate,...
                         parset.ParameterNames);

elseif isa(parset, 'pset.Grid')
    for ii = 1:np
        subset{ii} = (unique(pset.points(ii,:)) - p0(ii))/pd(ii);
    end
    newParset = pset.Grid(subset,...
                         parset.rate,...
                         parset.ParameterNames);

elseif isa(parset, 'pset.Hull')
    newPoints = (parset.points - repmat(p0,1,nv))/diag(pd);
    newParset = pset.Hull(newPoints,...
                         parset.rate,...
                         parset.ParameterNames);
    
elseif isa(parset, 'pset.Gral')
    newPoints = (parset.points - repmat(p0,1,nv))/diag(pd);
    newParset = pset.Gral(newPoints,...
                         parset.rate,...
                         parset.ParameterNames);

end

% matrix modification
A(:,:,1) = pdG.A(:,:,1);
B(:,:,1) = pdG.B(:,:,1);
C(:,:,1) = pdG.C(:,:,1);
D(:,:,1) = pdG.D(:,:,1);
for ii = 2:nv
    A(:,:,1) = A(:,:,1) + p0(ii-1)*pdG.A(:,:,ii);
    B(:,:,1) = B(:,:,1) + p0(ii-1)*pdG.B(:,:,ii);
    C(:,:,1) = C(:,:,1) + p0(ii-1)*pdG.C(:,:,ii);
    D(:,:,1) = D(:,:,1) + p0(ii-1)*pdG.D(:,:,ii);
    
    A(:,:,ii) = pd(ii-1)*pdG.A(:,:,ii);
    B(:,:,ii) = pd(ii-1)*pdG.B(:,:,ii);
    C(:,:,ii) = pd(ii-1)*pdG.C(:,:,ii);
    D(:,:,ii) = pd(ii-1)*pdG.D(:,:,ii);

end


pdGn = pass(A, B, C, D, newParset,...
            'StateName',pdG.StateName,...
            'InputName',pdG.InputName,...
            'OutputName',pdG.OutputName,...
            'InputGroup',pdG.InputGroup,...
            'OutputGroup',pdG.OutputGroup);
