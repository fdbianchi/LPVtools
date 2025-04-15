function bool = isequal(psetA,psetB)

% ISEQUAL returns true is both sets have the same points, the same 
% parameter rate limits and the same parameter names, and false otherwise.
%
% See also PSET.Box, PSET.Grid, PSET.Hull, PSET.Gral

% fbianchi - 2021-05-03

tol = 10*eps;

[npA,nvA] = size(psetA);
[npB,nvB] = size(psetB);

if (npA ~= npB) || (nvA ~= nvB)
    % checking dimensions
    bool = false;
    
else
    diffPoints = max(max(abs(psetA.points - psetB.points)));
    boolRange = diffPoints <= tol;
    
    rateTest = [isempty(psetA.rate) isempty(psetB.rate)];
    if sum(rateTest) == 2
        boolRate = true;
        
    elseif sum(rateTest) == 1
        boolRate = false;
        
    else
        diffRate = max(max(abs(psetA.rate - psetB.rate)));
        boolRate = diffRate <= tol;
        
    end
    
    boolNames = all(ismember(psetA.ParameterNames,psetB.ParameterNames));

    bool = boolRange && boolRate && boolNames;

end
