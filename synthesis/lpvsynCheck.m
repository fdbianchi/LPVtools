function [bool,constraint,obj] = lpvsynCheck(synProb,pv,options)

% LPVSYNCHECK checks if a synthesis problem is feasible in a denser grid of
% parameter. This is used to implement a gridding synthesis scheme.
%
% Use:
%
% Initial controller design with a grid of n points
% [pdK,constOut,obj,synSet] = lpvsyn(pdG,y,u,const,lyapSet);
% 
% then checking in a denser grid
% pvDenser = pset.Grid(range,2*n);
% [bool,constraint,obj] = lpvsynCheck(synSet,pvDenser);
%
% if bool=1 then the controller satisfies the constraints in the denser
% grid and pdK is a satisfactory controller.
%
% See also lpvsyn

% fbianchi - 2021-06-14

if (nargin < 3)
    options = lpvsettings;
end

% synthesis problem information
synSet  = synProb.synSet;

lyapSet = synProb.lyapSet;
if strcmp(lyapSet.type(1:3),'pwa')
    error('LPVSYNCHECK:inputerror','LPVSYNCHECK not available for PWA Lyapunov functions')
elseif strcmp(lyapSet.type(1:3),'aff')
    error('LPVSYNCHECK:inputerror','LPVSYNCHECK not available for affine Lyapunov functions')
end
lyapSet.parset = pv;

ctrlSet = synProb.ctrlSet;
ctrlSet.parset = pv;

constraint = synProb.const;

pdG = synProb.pdG;
synSet.LPVatPoints = ss(pdG,pv.points);


% -------------------------------------------------------------------------
% List of constraint by type
constClasses = arrayfun(@(x)class(x),constraint,'UniformOutput',false);
constSet.Gain.idx = find(ismember(constClasses,'synConst.Gain'));
constSet.GainH2.idx = find(ismember(constClasses,'synConst.GainH2'));
constSet.GainH2g.idx = find(ismember(constClasses,'synConst.GainH2g'));
constSet.Poles.idx = find(ismember(constClasses,'synConst.Poles'));

% -------------------------------------------------------------------------
% LMI constraints:
lmis = []; obj = 0;
lmiChecked = [];

% Gain - Hinf constrains
for ic = constSet.Gain.idx

    [lmiHinf,objHinf,lmiHinfChecked] = constHinf2LMI(pdG,synSet,lyapSet,ctrlSet,constraint(ic),options);
    lmis = [lmis, lmiHinf];
    lmiChecked = [lmiChecked, lmiHinfChecked];
    obj = obj + objHinf;
    
end

% Gain - H2 constrains
for ic = constSet.GainH2.idx

    [lmiH2,objH2] = constH22LMI(synSet,pdG,lyapSet,ctrlSet,constraint(ic),options);
    lmis = [lmis, lmiH2];
    obj = obj + objH2;
    
end

% Gain - H2g constrains
for ic = constSet.GainH2g.idx

    [lmiH2g,objH2g] = constH2g2LMI(synSet,pdG,lyapSet,ctrlSet,constraint(ic),options);
    lmis = [lmis, lmiH2g];
    obj = obj + objH2g;
    
end

% Poles constrains
for ic = constSet.Poles.idx

    [lmiPoles, lmiPolesChecked] = constPoles2LMI(pdG,synSet,lyapSet,ctrlSet,constraint(ic),options);
    lmiChecked = [lmiChecked, lmiPolesChecked];
    lmis = [lmis, lmiPoles];
    
end

% common to all constraints
[lmiPos,objPos,lmiPosChecked] = constPos2LMI(pdG,synSet,lyapSet,ctrlSet,options);
lmiChecked = [lmiChecked, lmiPosChecked];
lmis = [lmis, lmiPos];
obj = obj + objPos;

% -------------------------------------------------------------------------
% Solver options: 
if strcmp(options.solver, 'sedumi')
    solverOpts = sdpsettings('solver', 'sedumi',...
                  'sedumi.eps', options.solTol,...
                  'sedumi.cg.qprec', 1,... 
                  'sedumi.cg.maxiter', 49, ...
                  'sedumi.stepdif', 2,...
                  'sedumi.maxiter', options.maxiter,...
                  'verbose', options.verb,...
                  'dualize', options.dualize,...
                  'removeequalities', options.removeequalities);
    
elseif strcmp(options.solver, 'sdpt3')
    solverOpts = sdpsettings('solver', 'sdpt3',...
                  'sdpt3.steptol', options.solTol,...
                  'sdpt3.maxit', options.maxiter,...
                  'verbose', options.verb,...
                  'dualize', options.dualize,...
                  'removeequalities', options.removeequalities);
    
elseif strcmp(options.solver, 'mosek')
    solverOpts = sdpsettings('solver', 'mosek',...
                  'verbose', options.verb,...
                  'dualize', options.dualize,...
                  'removeequalities', options.removeequalities);
end


% -------------------------------------------------------------------------
% Solving optimization problem
diagnostic = solvesdp(lmis,obj,solverOpts);

if ~((diagnostic.problem == 0) || (diagnostic.problem == 4))
    % Infeasible
    disp(yalmiperror(diagnostic.problem));
    obj = Inf; 
    for ic = constSet.Gain.idx
        if isa(constraint(ic).bound,'sdpvar')
            constraint(ic).bound = inf;
        end
    end
    bool = false;
                
else
    % feasible
    obj = value(obj);
    for ic = constSet.Gain.idx
        if isa(constraint(ic).bound,'sdpvar')
            constraint(ic).bound = sqrt(value(constraint(ic).bound));
        end
    end
    % Gain - H2 constrains
    for ic = constSet.GainH2.idx
        if isa(constraint(ic).bound,'sdpvar')
            constraint(ic).bound = sqrt(value(constraint(ic).bound));
        end
    end
    % Gain - H2g constrains
    for ic = constSet.GainH2g.idx
        if isa(constraint(ic).bound,'sdpvar')
            constraint(ic).bound = sqrt(value(constraint(ic).bound));
        end
    end
    bool = all(cell2mat({lmiChecked(:).lmiValue}));
end



