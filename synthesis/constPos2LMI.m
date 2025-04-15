 function [lmi,obj,lmisChecked] = constPos2LMI(pdG,synSet,lyapSet,ctrlSet,options)

% *** for internal use ***
%
% CONSTPOS2LMI returns a set of LMI for the positiveness of Xcl
%
% Use: 
%   [lmi,obj] = constPos2LMI(pdG,synSet,lyapSet,eigtol)

% fbianchi - 2020-04-24

% number of systems in the description
nv = size(pdG.parset.points,2);
% system order
ns = order(pdG);

% optimization settings
eigtol  = options.eigtol;
penalty = options.penalty;
beta    = options.XYpenalty;

lmi = [];
lmisChecked = [];
obj = 0;
switch ctrlSet.fback         % type of feedback
    
    case 'state'            % state feedback
        
        switch lyapSet.type
            case 'cte'
                lmiName = 'Xcl>0: Xcl=cte - SF';
                Y = evalLyapFcn(lyapSet);
                lmi  = [(Y >= eigtol*eye(ns)):lmiName];
                obj  = penalty*norm(Y);

            case 'pwa'
                % ToDo: check Lyapunov LMI for case PWA

            otherwise

                for ii = 1:nv
                    Y = evalLyapFcn(lyapSet,ii);
                    lmiName = sprintf('Xcl>0: Xcl(p) @(%d) - SF',ii);
                    lmi  = [lmi, (Y >= eigtol*eye(ns)):lmiName];
                    obj  = obj + penalty*norm(Y);
                end

        end
        
    case 'output'           % output feedback
    
        switch lyapSet.type
            case 'cte'
                % constant function
                [Y, X] = evalLyapFcn(lyapSet);
                auxMat = [X beta*eye(ns);beta*eye(ns) Y];
                lmiName = 'Xcl>0 : Xcl=cte';
                if isnumeric(auxMat)
                    aux.lmiValue = all(eig(auxMat) >= 0);
                    aux.lmiName  = lmiName;
                    lmisChecked = [lmisChecked, aux];
                else
                    lmi = [lmi, (auxMat >= eigtol*eye(2*ns)):lmiName];
                    obj = penalty*(norm(X) + norm(Y));
                end

            case 'pwa'
                % PWA function
                for ii = 1:nv
                    [X, Y] = evalLyapFcn(lyapSet,ii);
                    auxMat = [X beta*eye(ns);beta*eye(ns) Y];
                    lmiName = sprintf('Xcl>0: Xcl PWA @(%d)',ii);
                    if isnumeric(auxMat)
                        aux.lmiValue = all(eig(auxMat) >= 0);
                        aux.lmiName  = lmiName;
                        lmisChecked = [lmisChecked, aux];
                    else
                        lmi = [lmi, (auxMat >= eigtol*eye(2*ns)):lmiName];
                        obj = obj + penalty*(norm(X) + norm(Y))/nv;
                    end                   
                end

            otherwise
                % affine or general parameter dependence
                for ii = 1:nv
                    [Y, X] = evalLyapFcn(lyapSet,ii);
                    auxMat = [X beta*eye(ns);beta*eye(ns) Y];
                    lmiName = sprintf('Xcl>0: Aff/Xcl(p) @(%d)',ii);
                    if isnumeric(auxMat)
                        aux.lmiValue = all(eig(auxMat) >= 0);
                        aux.lmiName  = lmiName;
                        lmisChecked = [lmisChecked, aux];
                    else
                        lmi = [lmi, (auxMat >= eigtol*eye(2*ns)):lmiName];
                        obj = obj + penalty*(norm(X) + norm(Y))/nv;
                    end
                end

        end
        
    otherwise
        error('LPVSYN:constPos2LMI:LMI',...
            'Type of feedback invalid')
    
end


