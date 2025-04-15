function str = char(obj)

% CHAR method for ppss class

% fbianchi - 2021-03-31

if isempty(obj)
    str = sprintf('\tEmpty polytopic LPV model.');
    
else
    % system dimensions
    [ny,nu,ns,~,nv] = size(obj);
    
    str = sprintf(['Polytopic LPV model of order %.0f and %.0f vertices\n',...
                   '\t with %.0f input(s) and %.0f output(s)\n'],...
                   ns,nv,nu,ny);

end
