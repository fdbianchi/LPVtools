
% LPVtools: Example of current loop control of a PMSM: design of a LPV PID
% 
%   Plant:          general LPV (ppss)
%                   1 parameter,
%                   2 inputs, 2 outputs
%   Constraint:     Hinf/Poles
%   Lyapunov Fcn:   Constant/Parameter dependent
%   Ctrller Fcn:    No specified

% fbianchi - 2021-04-05


% cleaning
clearvars
close all
clc

% ========================================================================
% Modelling

% model parameters
Rc = 0.3;  Lc = 4.6e-3;
Rg = 0.22; Lg = 2.9e-3; phi = 0.2591; 
R  = Rc+Rg; L = Lc+Lg; Tsw = 1/12e3;

% parameter envelope
pv  = pset.Box([200 1600]*pi/10);

% lpv model
A(:,:,1) = [-R/L 0; 0 -R/L]; 
A(:,:,2) = [0 -1; 1 0]; 
B = [-1/L 0; 0 -1/L];
C = eye(2); D = zeros(2);

pdG = pass(A,B,C,D,pv);
pdG.y = 'y';
pdG.u = 'u';

% ========================================================================
% Control design
%

pdG = ppss(pdG);

% vertex 1
[a,b,c] = ssdata(pdG(1));
A   = [a, zeros(2);-c zeros(2)];
B1  = [zeros(2);eye(2)];
B2  = [b; zeros(2)];
C1  = [zeros(2), eye(2); zeros(2,4)];
% C2  = [zeros(2), eye(2); eye(2) zeros(2)];
C2  = eye(4);
D11 = zeros(4,2);
D12 = [zeros(2); eye(2)];
D21 = zeros(4,2);
D22 = zeros(4,2);
Sa(:,:,1) = ss(A,[B1 B2],[C1;C2],[D11, D12; D21,D22]);
% vertex 2
[a,b,c] = ssdata(pdG(2));
A = [a, zeros(2);-c zeros(2)];
Sa(:,:,2) = ss(A,[B1 B2],[C1;C2],[D11, D12; D21,D22]);

pdGau = ppss(Sa,pv);
pdGau.u = {'r(1)','r(2)','u(1)','u(2)'};
pdGau.y = {'us(1)','us(2)','u(1)','u(2)','y(1)','y(2)','y(3)','y(4)'};
        
% weighting functions
We = ss(0.1);
Wu = ss(0.1);
Wout = append(We,We,Wu,Wu,eye(4));
pdGw = Wout*pdGau;
pdGw.u = {'r(1)','r(2)','u(1)','u(2)'};
pdGw.y = {'us(1)','us(2)','u(1)','u(2)','y(1)','y(2)','y(3)','y(4)'};


% Controller design
y = [];
u = {'u(1)','u(2)'};
z = {'us(1)','us(2)','u(1)','u(2)'};
w = {'r(1)','r(2)'};
r = [4 2];

const1   = synConst.Gain(1:2,1:4);

% State-feedback gain scheduling 
[pdKgs,const_gs] = lpvsyn(pdGw,0,3:4,const1);
g1 = const_gs.bound;
% State-feedback robust control 
[pdKrb,const_rb] = lpvsyn(pdGw,0,3:4,const1,[],0);
g2 = const_rb.bound;
% State-feedback gain scheduling + poles constraint
const2(1) = synConst.Gain(1:2,1:4);
d = 0.7;
const2(2) = synConst.Poles('MinDamping',d,'MaxFreq',500);
[pdKgs2,const_gs] = lpvsyn(pdGw,0,3:4,const2);
g1b = const_gs(1).bound;
% Output-feedback gain scheduling 
[pdKof,const_of] = lpvsyn(pdGw,5:8,3:4,const2);
g3 = const_of(1).bound;


% parameter dependant Lyapunov functions
%   lyapSet.parfcnX  = 0;       % lyapSet.dparfcnX = 0;
%   lyapSet.parfcnY  = @(p) p;  % lyapSet.dparfcnY = 0;
lyapSet = createLyapFcn(0,0,@(p) p,0);
[pdK0,const_gs] = lpvsyn(pdGw,y,u,const2,lyapSet);
g1c = const_gs(1).bound;

% analysis
pdGclgs = lft(pdGw(1:4,:),pdKgs);
[constA,objA] = lpvanalysis(pdGclgs,const1);
g1_A = constA(1).bound;

pdGclrb = lft(pdGw(1:4,:),pdKrb);
[constA,objA] = lpvanalysis(pdGclrb,const1);
g2_A = constA(1).bound;

pdGclof = lft(pdGw,pdKof);
[constA,objA] = lpvanalysis(pdGclof,const1);
g3_A = constA(1).bound;

pdGclgs2 = lft(pdGw(1:4,:),pdKgs2);
[constA,objA] = lpvanalysis(pdGclgs2,const2);
g1b_A = constA(1).bound;

pdGcl0 = lft(pdGw(1:4,:),pdK0);
lyapSet = createLyapFcn('cl',@(p) p);
lyapSet.inv = true;
pdGcl0 = ppss(pdGcl0,pv);
[constA,objA] = lpvanalysis(pdGcl0,const2(1),lyapSet,pv.points);
g1c_A = constA(1).bound;

% results
fprintf('\n-------------------------------------------------------------------\n')
fprintf('Comparison for Hinf synthesis\n')
fprintf('Gain-scheduling (state-feedback):          %5.4f (analysis %5.4f)\n',g1,g1_A)
fprintf('Robust (state-feedback):                   %5.4f (analysis %5.4f)\n',g2,g2_A)
fprintf('Gain-scheduling (output-feedback):         %5.4f (analysis %5.4f)\n',g3,g3_A)
fprintf('Gain-scheduling (state-feedback+poles):    %5.4f (analysis %5.4f)\n',g1b,g1b_A)
fprintf('Gain-scheduling (state-fdback+poles+Y(p)): %5.4f (analysis %5.4f)\n',g1c,g1c_A)
fprintf('-------------------------------------------------------------------\n')
        
% Checking performance at vertices
Gaw = ss(pdGw);
Kgs = ss(pdKgs);
Krb = ss(pdKrb);
Kof = ss(pdKof);

Gcl_sf = lft(Gaw, Kgs);
Gcl_rb = lft(Gaw, Krb);
Gcl_of = lft(Gaw, Kof);

w = logspace(0, 6, 1000);
for ii = 1:2
    aux = sigma(Gcl_sf(:,:,ii), w);
    sg_sf(ii,:) = aux(1,:);
    aux = sigma(Gcl_rb(:,:,ii), w);
    sg_rb(ii,:) = aux(1,:);
    aux = sigma(Gcl_of(:,:,ii), w);
    sg_of(ii,:) = aux(1,:);
end

cl = lines(3);
hl(1:2) = semilogx(w, sg_sf, 'Color', cl(1,:));
hold on
hl(3:4) = semilogx(w, sg_rb, 'Color', cl(2,:));
hl(5:6) = semilogx(w, sg_of, 'Color', cl(3,:));
xlabel('Freq. (rad/sec)')
ylabel('Sigma')

legend(hl(1:2:6),'State Feedback', 'Robust', 'Output Feedback')

