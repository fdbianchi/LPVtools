function pdGn = normalize(pdG)

% NORMALIZE transforms the model pdG into a model with functions between -1 and 1
%
% Use:
%   pdGn = NORMALIZE(pdG)
%
% See also denormalize

% fbianchi - 2024-01-16


% model dimensions
[~,~,~,~,nv] = size(pdG);

% normalization of the parameter functions
range = pdG.parset.range;
f0 = zeros(nv - 1,1);
fd = zeros(nv - 1,1);
faux =@(x,i) x(i);
for ii = 1:nv-1
    fmin = fminbnd(@(x) faux(pdG.parfcn(x),1), range(ii,1), range(ii,2));
    fmax = fminbnd(@(x) faux(-pdG.parfcn(x),1), range(ii,1), range(ii,2));
    
    f0(ii) = (fmax + fmin)/2;
    fd(ii) = (fmax - fmin)/2;
    
end
newfcn =@(p) (pdG.parfcn(p) - f0)./fd;

% matrix modification
A(:,:,1) = pdG.A(:,:,1);
B(:,:,1) = pdG.B(:,:,1);
C(:,:,1) = pdG.C(:,:,1);
D(:,:,1) = pdG.D(:,:,1);
for ii = 2:nv
    A(:,:,1) = A(:,:,1) + f0(ii-1)*pdG.A(:,:,ii);
    B(:,:,1) = B(:,:,1) + f0(ii-1)*pdG.B(:,:,ii);
    C(:,:,1) = C(:,:,1) + f0(ii-1)*pdG.C(:,:,ii);
    D(:,:,1) = D(:,:,1) + f0(ii-1)*pdG.D(:,:,ii);
    
    A(:,:,ii) = fd(ii-1)*pdG.A(:,:,ii);
    B(:,:,ii) = fd(ii-1)*pdG.B(:,:,ii);
    C(:,:,ii) = fd(ii-1)*pdG.C(:,:,ii);
    D(:,:,ii) = fd(ii-1)*pdG.D(:,:,ii);

end


pdGn = pgss(A, B, C, D, pdG.parset, newfcn,...
            'StateName',pdG.StateName,...
            'InputName',pdG.InputName,...
            'OutputName',pdG.OutputName,...
            'InputGroup',pdG.InputGroup,...
            'OutputGroup',pdG.OutputGroup);
