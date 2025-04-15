function str = char(obj)

% CHAR method for pcss class

% fbianchi - 2020-07-06

if isempty(obj)
    str = sprintf('\tEmpty pcss class.');
    
else
    % system dimensions
    [ny,nu,ns,np] = size(obj);
    
    str1 = sprintf('On-line implemented LPV controller of order %.0f\n',ns);
    str2 = sprintf('\twith %.0f input(s), %.0f output(s) and %.0f parameter(s)\n',nu,ny,np);
    str  = [str1 str2];

end
