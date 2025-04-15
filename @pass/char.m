function str = char(obj)

% CHAR method for pass class

% fbianchi - 2021-03-31

if isempty(obj)
    str = sprintf('\tEmpty affine LPV model.');
    
else
    % system dimensions
    [no,ni] = iosize(obj);
    [np,~] = size(obj.parset);
    ns = order(obj);
    
    str = sprintf(['Affine LPV model of order %.0f and %.0f parameter(s)\n',...
                   '\twith %.0f input(s) and %.0f output(s).\n'],...
                   ns,np,ni,no);
    
end
