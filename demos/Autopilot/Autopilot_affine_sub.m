
% LPVtools: Autopilot Example (see "LMI control toolbox manual", pp 7-10)
%
%   Plant:          affine LPV (psys/pass)
%   Constraint:     Hinf/Poles
%   Lyapunov Fcn:   Constant/Parameter dependent
%   Ctrller Fcn:    No specified
%
% using suboptimal option to compute a suboptimal controller with better
% numerical representation

% fbianchi - 2024-03-11


% cleanig 
clc
clearvars
close all

% ========================================================================
% Modelling

% Operating range
Zmin =  0.5;   Zmax =   4;
Mmin =  0.0;   Mmax = 106;  
% system matrices
a0 = [ 0 1; 0 0]; a1 = [-1 0; 0 0]; a2 = [ 0 0;-1 0];
b0 = [0;1];       b1 = [0;0];       b2 = [0;0];
c0 = [-1 0;0 1];  c1 = [ 0 0; 0 0]; c2 = [ 0 0; 0 0];
d0 = [0;0];       d1 = [0;0];       d2 = [0;0];

% ------------------------------------------------------------------------
% Using lmitool

% parameter set
pv = pvec('box',[Zmin Zmax; Mmin Mmax]);
Vx = polydec(pv);
% Affine model:
s0 = ltisys(a0,b0,c0,d0);
s1 = ltisys(a1,b1,c1,d1,0); % Z_al 
s2 = ltisys(a2,b2,c2,d2,0); % M_al 
pdGlmi = psys(pv,[s0 s1 s2]);

% ------------------------------------------------------------------------
% Using LPVtool

% system matrices
A = cat(3,a0,a1,a2);
B = cat(3,b0,b1,b2);
C = cat(3,c0,c1,c2);
D = cat(3,d0,d1,d2);
% parameter set
Pset = pset.Box([Zmin Zmax; Mmin Mmax],[-10 10],{'Z','M'});

% affine lpv model
pdG = pass(A,B,C,D,Pset,'InputName','u','OutputName','y');


% ========================================================================
% Control design

% Weight on S
nf1 = 2.0101;  df1 = [1.0000e+00   2.0101e-01];
% Weight on KS
nf2 = [9.6785e+00   2.9035e-02   0            0];
df2 = [1.0000e+00   1.2064e+04   1.1360e+07   1.0661e+10];

% ------------------------------------------------------------------------
% Using lmitool

W1lmi = ltisys('tf',nf1,df1);
W2lmi = ltisys('tf',nf2,df2);

% augmented plant
inputs  = 'r';
outputs = 'e=r-G(1);K';
K_in    = 'K:e;G(2)';
G_in    = 'G:K';
[pdGauLmi,r] = sconnect(inputs,outputs,K_in,G_in,pdGlmi);

% augmented plant + weight
pdGawLmi = smult(pdGauLmi,sdiag(W1lmi,W2lmi,eye(2)));

% H-infinity constraint
[gLmi,pdkLmi] = hinfgs(pdGawLmi,r);
pdKLmi = addpv(pdkLmi,pvec('pol',Vx));


pdGclLmi = slft(pdGawLmi,pdkLmi);
[gLmi_a,P] = quadperf(pdGclLmi);

% ------------------------------------------------------------------------
% Using LPVtool (more options)

W1 = mat2lti(W1lmi);
W2 = mat2lti(W2lmi);

% augmented plant
er  = sumblk('e = r - y(1)');
pdGau = connect(pdG,er,{'r','u'},{'e','u','e','y(2)'});
% augmented plant + weight
Wout  = balreal(append(W1,W2,1,1));   
% Wout  = append(W1,W2,1,1);   
Wout.y = {'et','ut','e','y(2)'};
pdGaw = Wout*pdGau;
% pdGaw = ppss(pdGawLmi);
pdGaw.y = {'et','ut','e','y(2)'};
pdGaw.u = {'r','u'};


% computing suboptimal controllers
opt = lpvsettings('solver','mosek','subOpt',true,'XYpenalty',1.6,'Ctrlpenalty',1e-5);

% H-infinity constraint
const(1) = synConst.Gain('r',{'et','ut'});
const(2) = synConst.Poles('MaxFreq',1e5);

% --------------------------
% constant Lyapunov function
[pdK1,constOut] = lpvsyn(pdGaw,r);
ghinf_1 = constOut(1).bound;
[pdK2,constOut] = lpvsyn(pdGaw,{'e','y(2)'},'u',const,[],[],opt);
ghinf_2 = constOut(1).bound;

pdGcl2 = lft(pdGaw,pdK2);
constA = lpvanalysis(pdGcl2,const,[],[],opt);
ghinf_2_A = constA(1).bound;

% --------------------------
% X & Y parameter dependent, dX=0, dY=0
[pdK3a,constOut] = lpvsyn(pdGaw,{'e','y(2)'},'u',const,'aff',[],opt);
ghinf_3a = constOut(1).bound;
% 
% when using @-function, we need to used griding
Pset_gr = pset.Grid(pdGaw.parset.range,2,pdGaw.parset.rate);
pdGaw_gr = pdGaw;
pdGaw_gr.parset = Pset_gr;
lyapSet = createLyapFcn(@(p) p,0,@(p) p,0);
[pdK3,constOut]  = lpvsyn(pdGaw_gr,{'e','y(2)'},'u',const,lyapSet,[],opt);
ghinf_3 = constOut(1).bound;

% --------------------------
% X parameter dependent & Y cte, dX~=0, dY=0
[pdK4a,constOut] = lpvsyn(pdGaw,{'e','y(2)'},'u',const,'affdX',[],opt);
ghinf_4a = constOut(1).bound;
%
% when using @-function, we need to used griding
lyapSet = createLyapFcn(@(p) p,@(p) eye(2));
[pdK4,constOut]  = lpvsyn(pdGaw_gr,{'e','y(2)'},'u',const,lyapSet,[],opt);
ghinf_4 = constOut(1).bound;

% --------------------------
% X cte & Y parameter dependent, dX=0, dY~=0
[pdK5a,constOut] = lpvsyn(pdGaw,{'e','y(2)'},'u',const,'affdY',[],opt);
ghinf_5a = constOut(1).bound;
% 
% when using @-function, we need to used griding
lyapSet = createLyapFcn(0,0,@(p) p,@(p) eye(2));
[pdK5,constOut]  = lpvsyn(pdGaw_gr,{'e','y(2)'},'u',const,lyapSet,[],opt);
ghinf_5 = constOut(1).bound;


% ========================================================================
% Evaluation/Validation

fprintf('\n')
fprintf('===================================================================\n')
fprintf('Comparison for Hinf synthesis\n')
fprintf('hinfgs:                      %5.4f\n',gLmi)
fprintf('lpvsyn(pdG,r):               %5.4f\n',ghinf_1)
fprintf('lpvsyn(pdG,u,y,const):       %5.4f\n',ghinf_2)
fprintf('lpvsyn(pdG,u,y,const):       %5.4f (analysis)\n',ghinf_2_A)
fprintf('------------------------------------------------------------------\n')
fprintf('lpvsyn(pdG,u,y,const,lyap):  %5.4f (X=X(p), dX=0,  Y=Y(p), dY=0)\n',ghinf_3)
fprintf('lpvsyn(pdG,u,y,const,lyap):  %5.4f (aff)\n',ghinf_3a)
fprintf('------------------------------------------------------------------\n')
fprintf('lpvsyn(pdG,u,y,const,lyap):  %5.4f (X=X(p), dX~=0, Y=Y0,   dY=0)\n',ghinf_4)
fprintf('lpvsyn(pdG,u,y,const,lyap):  %5.4f (affdX)\n',ghinf_4a)
fprintf('------------------------------------------------------------------\n')
fprintf('lpvsyn(pdG,u,y,const,lyap):  %5.4f (X=X0,   dX=0,  Y=Y(p), dY~=0)\n',ghinf_5)
fprintf('lpvsyn(pdG,u,y,const,lyap):  %5.4f (affdY)\n',ghinf_5a)
fprintf('==================================================================\n')

