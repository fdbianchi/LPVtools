function pdG = denormalize(pdGn, fcn)

% DENORMALIZE transforms the model pdGn with parameter funtions ranging  
% between -1 and 1 into one with a new functions FCN
%
% Use:
%   pdG = DENORMALIZE(pdGn, fcn)
%
% See also normalize

% fbianchi - 2024-01-16

if ~isa(fcn,'function_handle')
    error('PGSS:NORMALIZE:inputError','FCn must be an anonymous function')
end

% model dimensions
[~,~,~,~,nv] = size(pdGn);

% new parameter functions
range = pdGn.parset.range;
f0 = zeros(nv - 1,1);
fd = zeros(nv - 1,1);
faux =@(x,i) x(i);
for ii = 1:nv-1
    fmin = fminbnd(@(x) faux(fcn(x),1), range(ii,1), range(ii,2));
    fmax = fminbnd(@(x) faux(-fcn(x),1), range(ii,1), range(ii,2));
    
    f0(ii) = (fmax + fmin)/2;
    fd(ii) = (fmax - fmin)/2;
    
end

% matrix modification
A(:,:,1) = pdGn.A(:,:,1);
B(:,:,1) = pdGn.B(:,:,1);
C(:,:,1) = pdGn.C(:,:,1);
D(:,:,1) = pdGn.D(:,:,1);
for ii = 2:nv
    fact0 = f0(ii-1)/fd(ii-1);
    A(:,:,1) = A(:,:,1) - fact0*pdGn.A(:,:,ii);
    B(:,:,1) = B(:,:,1) - fact0*pdGn.B(:,:,ii);
    C(:,:,1) = C(:,:,1) - fact0*pdGn.C(:,:,ii);
    D(:,:,1) = D(:,:,1) - fact0*pdGn.D(:,:,ii);
    
    factd = 1/fd(ii-1);
    A(:,:,ii) = factd*pdGn.A(:,:,ii);
    B(:,:,ii) = factd*pdGn.B(:,:,ii);
    C(:,:,ii) = factd*pdGn.C(:,:,ii);
    D(:,:,ii) = factd*pdGn.D(:,:,ii);

end


pdG = pgss(A, B, C, D, pdGn.parset, fcn,...
            'StateName',pdGn.StateName,...
            'InputName',pdGn.InputName,...
            'OutputName',pdGn.OutputName,...
            'InputGroup',pdGn.InputGroup,...
            'OutputGroup',pdGn.OutputGroup);
