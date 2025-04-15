function usys = uss(obj)

% USS(pdG) converts a pass object into a uss object
%
% USS(pdG) converts an affine LPV model into an uncertain model
%   in which the parameters are considered as uncertain real parameters
%
% Use:
%   usys = USS(pdG)
%
% See also pass, uss

% fbianchi - 2021-03-31
    
% independent terms
A = obj.A(:,:,1); B = obj.B(:,:,1);
C = obj.C(:,:,1); D = obj.D(:,:,1);

% Nominal and range of p
np    = size(obj.A,3) - 1;
pval  = obj.parset.points;
pnom  = mean(pval,2);
range = obj.parset.range;

for ii = 1:np
    % each parameter is converted into ureal object
    par = ureal(obj.parset.ParameterNames{ii},pnom(ii,:),...
                'range',range(ii,:));
    
    % building the parameter matrices
    if any(any(obj.A(:,:,ii+1)))
        A = A + obj.A(:,:,ii+1)*par; 
    end
    if any(any(obj.B(:,:,ii+1)))
        B = B + obj.B(:,:,ii+1)*par; 
    end
    if any(any(obj.C(:,:,ii+1)))
        C = C + obj.C(:,:,ii+1)*par; 
    end
    if any(any(obj.D(:,:,ii+1)))
        D = D + obj.D(:,:,ii+1)*par; 
    end
end

% creating the uss system
usys = ss(A,B,C,D,...
        'StateName',obj.StateName,...
        'InputName',obj.InputName,...
        'OutputName',obj.OutputName,...
        'InputGroup',obj.InputGroup,...
        'OutputGroup',obj.OutputGroup);

