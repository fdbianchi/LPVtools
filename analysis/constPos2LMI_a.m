function [lmi,obj] = constPos2LMI_a(pdG,lyapSet,options)

% *** for internal use ***
%
% CONSTPOS2LMI_A returns a set of LMI for the positiveness of Xcl
%
% Use: 
%   [lmi,obj] = constPos2LMI_A(pdG,lyapSet,options)

% fbianchi - 2021-07-05

% number of systems in the description
nv = size(lyapSet.LPVatPoints,3);

% system order
ns = order(pdG);

% optimization settings
eigtol  = options.eigtol;
penalty = options.penalty;

lmi = [];
obj = 0;

switch lyapSet.type
    
    case 'cte'
        lmiName = 'Xcl>0: Xcl=cte';
        X = evalLyapFcn_a(lyapSet);
        lmi  = [(X >= eigtol*eye(ns)):lmiName];
        obj  = penalty*norm(X);
        
    otherwise
        
        for ii = 1:nv
            X = evalLyapFcn_a(lyapSet,ii);
            lmiName = sprintf('Xcl>0: Xcl(p) @(%d)',ii);
            lmi  = [lmi, (X >= eigtol*eye(ns)):lmiName];
            obj  = obj + penalty*norm(X)/nv;
        end
        
        
end
        
