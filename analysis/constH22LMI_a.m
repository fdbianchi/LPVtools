function [lmis,obj,Q] = constH22LMI_a(LPVatPoints,lyapSet,constraint,options)

% *** for internal use ***
%
% CONSTHINF2LMI_A returns a set of LMIs for H2 performance constraint
%
% Use: 
%   [lmis,obj] = CONSTH22LMI_A(LPVatPoints,lyapSet,constraint,optionss)

% fbianchi - 2021-07-05

% objective initial value
obj = 0;

% optimization settings
eigtol  = options.eigtol;

% number of systems in the description
nv = size(LPVatPoints,3);

% disturbance and performance input-output map
ioChannel.perf = constraint.outputs;
nz = length(ioChannel.perf);
ioChannel.dist = constraint.inputs;
nw = length(ioChannel.dist);

% additional variable
Q  = sdpvar(nz);

% mu is an optimization variable or not depending on soft or hard
% constraints
if isa(constraint.bound,'sdpvar')
    mu2 = constraint.bound;
else
    mu2 = constraint.bound^2;
end

% LMIs
lmis = [];

for ii=1:nv
    
    % system matrices at par
    [A,B,C,D] = ssdata(LPVatPoints(ioChannel.perf,ioChannel.dist,ii));
    
    if any(any(D))
        error('LPVANALYSIS:constH22LMI_A:LMI',...
            'The constraint D=0 is infeasible')
    end
    
    % Lyapunov matrices at par
    [X,dX] = evalLyapFcn_a(lyapSet,ii);
    
    for jj = 1:length(dX)
        
        % H2 constraints
        matB = blkvar;
        matB(1,1) = dX{jj} + X*A + (X*A)';
        matB(1,2) = X*B;
        matB(2,2) = -eye(nw);
        matB = sdpvar(matB);
        
        lmiName = sprintf('H2: LMIB @(%d)',ii);
        lmis = [lmis, (matB <= -eigtol*eye(size(matB))):lmiName];
        
    end
    
    lmiQ = blkvar;
    lmiQ(1,1) = X;
    lmiQ(1,2) = C';
    lmiQ(2,2) = Q;
    lmiQ = sdpvar(lmiQ);
    lmiName = sprintf('H2: LMIQ @(%d)',ii);
    lmis = [lmis, (lmiQ >= eigtol*eye(size(lmiQ))):lmiName];
    
end

lmiName = 'H2: Trace(Q) < mu';
lmis = [lmis, (trace(Q) <= mu2):lmiName];

% objective function
if isa(constraint.bound,'sdpvar')
    % soft constraint
    obj = obj + mean(constraint.bound)*constraint.factor;
  
end

