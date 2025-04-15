function usys = uss(obj)

% USS(pdG) converts a pgss object into a uss object
%
% USS(pdG) converts an general LPV model into an uncertain model
%   in which the fcn_i(p) are considered as uncertain real parameters
%
% Use:
%   usys = USS(pdG)
%
% See also pgss, uss

% fbianchi - 2021-04-05
    
% independent terms
A = obj.A(:,:,1); B = obj.B(:,:,1);
C = obj.C(:,:,1); D = obj.D(:,:,1);

% Nominal and range of fcn(p)
[~,nv] = size(obj.parset);
np = size(obj.A,3) - 1;
pval = zeros(np,nv);
for ii = 1:nv
    pval(:,ii) = obj.parfcn(obj.parset.points(:,ii));
end
pnom = mean(pval,2);
pmin = min(pval,[],2);
pmax = max(pval,[],2);

for ii = 1:np
    % each parameter is converted into ureal object
    par = ureal(sprintf('p%d',ii),pnom(ii,:),...
                'range',[pmin(ii,:),pmax(ii,:)]);
    
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
          'StateName', obj.StateName,...
          'InputName', obj.InputName,...
          'OutputName', obj.OutputName,...
          'InputGroup', obj.InputGroup,...
          'OutputGroup', obj.OutputGroup);


