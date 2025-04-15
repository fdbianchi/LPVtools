function str = char(obj)

% CHAR method for pgss class

% fbianchi - 2021-03-31

if isempty(obj)
    str = sprintf('\tEmpty general LPV model.');
    
else
    % system dimensions
    [no,ni] = iosize(obj);
    [np,nv] = size(obj.parset);
    ns = order(obj);
    nf = nsys(obj) - 1;

    str = sprintf(['General LPV model of order %.0f and %.0f parameter(s)',...
        '(%.0f function(s), %.0f point(s))\n\t with %.0f ',...
        'input(s) and %.0f output(s).\n'],...
        ns,np,nf,nv,ni,no);
    
end
