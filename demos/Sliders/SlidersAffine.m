
% LPVtools: Coupled Sliders on Rotating Rods (see Becker's Thesis, pp 89)
%
%   Plant:          affine LPV (psys/pass)
%   Constraint:     Hinf/Poles
%   Lyapunov Fcn:   Constant/Parameter dependent
%   Ctrller Fcn:    No specified

% fbianchi - 2021-04-01


% cleaning
clearvars
close all
clc

% ========================================================================
% Modelling

% System data
M1 = 1;     % mass first slider (kg)
M2 = 0.5;   % mass second slider (kg)
b  = 1;     % damping coefficient in the slots (kg/sec)
k  = 200;   % spring constant (N/m)
%
% Varying parameters
%   - p1: squared angular velocity of the rod 1
%   - p2: squared angular velocity of the rod 2
P1_min = 0;   P1_max =  9;
P2_min = 0;   P2_max = 25;
Prange = [P1_min P1_max;P2_min P2_max];
Prate  = [-5 5;-5 5];

% ------------------------------------------------------------------------
% Using lmitool

pv = pvec('box',Prange);

% LPV plant
A0 = [0,       0,       1,      0;
      0,       0,       0,      1;
     -k/M1,   -k/M1,   -b/M1,   0;
     -k/M2,   -k/M2,    0,     -b/M2];
A1 = [0,       0,       0,      0;
      0,       0,       0,      0;
      1,       0,       0,      0;
      0,       0,       0,      0];
A2 = [0,       0,       0,      0;
      0,       0,       0,      0;
      0,       0,       0,      0;
      0,       1,       0,      0];
B = [0 0 0;0 0 0;1/M1 0.1/M1 0;0 0 0.1/M2];
C = [0 1 0 0];
D = [0 0 0];
s0 = ltisys(A0,B,C,D);
s1 = ltisys(A1,B*0,C*0,D*0,0);
s2 = ltisys(A2,B*0,C*0,D*0,0);

pdGlmi = psys(pv,[s0 s1 s2]);

% ------------------------------------------------------------------------
% Using LPVtool

Pset = pset.Box([P1_min P1_max;P2_min P2_max],[-5 5]);

A(:,:,1) = [0,       0,       1,      0;
            0,       0,       0,      1;
           -k/M1,   -k/M1,   -b/M1,   0;
           -k/M2,   -k/M2,    0,     -b/M2];
A(:,:,2) = [0,       0,       0,      0;
            0,       0,       0,      0;
            1,       0,       0,      0;
            0,       0,       0,      0];
A(:,:,3) = [0,       0,       0,      0;
            0,       0,       0,      0;
            0,       0,       0,      0;
            0,       1,       0,      0];
B = [0 0 0;0 0 0;1/M1 0.1/M1 0;0 0 0.1/M2];
C = [0 1 0 0];
D = [0 0 0];

pdG = pass(A,B,C,D,Pset);
pdG.u = {'ua','d1','d2'};
pdG.y = 'y';

% ========================================================================
% Control design

% ------------------------------------------------------------------------
% Using lmitool

% weigthing functions
We = ltisys('tf',[0.3 1.2],[1 0.04]);
Wn = ltisys('tf',[1 0.4],[0.01 400]);
Wu = ltisys('tf',[1 0.1],[0.01 125]);
Wa = 0.00001;
Act= ltisys('tf',1,[0.01 1]);
% augmented plant
inputs  = 'r;v;d(2)';
outputs = 'pdG(1)-r;K;Act';
K_in    = 'K:r;pdG(1)+v';
pdG_in  = 'pdG:Act;d';
Act_in  = 'Act:K';
[pdGauLmi,r] = sconnect(inputs,outputs,K_in,pdG_in,pdGlmi,Act_in,Act);
% augmented plant + weights
Wo = sdiag(We,Wu,Wa,1,1);
Wi = sdiag(1,Wn,1,1,1);
pdGauwLmi = smult(Wi,pdGauLmi,Wo);
% H-infinity constraint
[glmi,pdKlmi] = hinfgs(pdGauwLmi,r);
pdKlmi = addpv(pdKlmi,pv);

pdGclLmi = slft(pdGauwLmi,pdKlmi);
[glmi_A,P] = quadperf(pdGclLmi);

% ------------------------------------------------------------------------
% Using LPVtool

% weigthing functions
We = tf([0.3 1.2],[1 0.04]); We.u = 'e';   We.y = 'et';
Wu = tf([1 0.1],[0.01 125]); Wu.u = 'u';   Wu.y = 'ut';
Wa = tf(0.00001,1);          Wa.u = 'ua';  Wa.y = 'uat';
Wn = tf([1 0.4],[0.01 400]); Wn.u = 'n';   Wn.y = 'nt';
% augmented plant
Act= tf(1,[0.01 1]);         Act.u = 'u';  Act.y = 'ua';
s1 = sumblk('e = y - r');
s2 = sumblk('yn = n + y');
pdGau = connect(pdG,s1,s2,Act,{'r','n','d1','d2','u'},{'e','u','ua','r','yn'});
% to check
for ii=1:4
    G_o(:,:,ii) = mat2lti(psinfo(pdGauLmi,'sys',ii));
end
figure
bodemag(G_o,'r:',pdGau,'k-');

% augmented plant + weights
Wo = append(We,Wu,Wa,1,1); Wo.y = {'et','ut','at','r','yn'};
Wi = append(1,Wn,1,1,1);   Wi.u = {'r','n','d1','d2','u'};
pdGaw = Wo*pdGau*Wi;

% Controller synthesis: H-infinity
y = {'r','yn'};
u = {'u'};
z = {'et','ut','at'};
w = {'r','n','d1','d2'};
const(1) = synConst.Gain(w,z);
% const(2) = synConst.Poles('MaxFreq',1e5);
const(2) = synConst.Poles('MinDamping',0.09,'MaxFreq',1e5);

% -----------------------------------------------------
% constant Lyapunov function

[pdK1,constOut] = lpvsyn(pdGaw,r);
ghinf_1 = constOut.bound;
[pdK2,constOut] = lpvsyn(pdGaw,y,u,const);
ghinf_2 = constOut(1).bound;

% checking CL performance
pdGcl2 = lft(pdGaw,pdK2);
constA = lpvanalysis(pdGcl2,const);
ghinf_2_A = constA(1).bound;


% -----------------------------------------------------
% Using parameter-dependent Lyapunov functions

nsyn = 3;   % # points in the grid for synthesis
nchk = 5;   % # points in the grid for checking

% 1) X & Y parameter dependent, dX=0, dY=0
%    - using multiconvexity
[pdK3a,constOut] = lpvsyn(pdGaw,y,u,const,'aff');
ghinf_3a = constOut(1).bound;
%    - using gridding
PsetGrid = pset.Grid(Prange,nsyn,Prate);
pdGawGrid = pdGaw;
pdGawGrid.parset = PsetGrid;
lyapSet = createLyapFcn(@(p) p,0,@(p) p,0);
[pdK3,constOut,~,synSet] = lpvsyn(pdGawGrid,y,u,const,lyapSet);
ghinf_3 = constOut(1).bound;
% checking CL performance in a denser grid
pvDenser = pset.Grid(Prange,nchk,Prate);
[~,constChk] = lpvsynCheck(synSet,pvDenser);
ghinf_3_C = constChk(1).bound;

% 2) X parameter dependent & Y cte, dX~=0, dY=0
%    - using multiconvexity
[pdK4a,constOut] = lpvsyn(pdGaw,y,u,const,'affdX');
ghinf_4a = constOut(1).bound;
%    - using gridding
lyapSet = createLyapFcn(@(p) p,@(p) eye(2));
[pdK4,constOut,~,synSet] = lpvsyn(pdGawGrid,y,u,const,lyapSet);
ghinf_4 = constOut(1).bound;
% checking CL performance in a denser grid
[~,constChk] = lpvsynCheck(synSet,pvDenser);
ghinf_4_C = constChk(1).bound;

% 3) X cte & Y parameter dependent, dX=0, dY~=0
%    - using multiconvexity
[pdK5a,constOut] = lpvsyn(pdGaw,y,u,const,'affdY');
ghinf_5a = constOut(1).bound;
%    - using gridding
lyapSet = createLyapFcn(0,0,@(p) p,@(p) eye(2));
[pdK5,constOut,~,synSet] = lpvsyn(pdGawGrid,y,u,const,lyapSet);
ghinf_5 = constOut(1).bound;
% checking CL performance in a denser grid
[~,constChk] = lpvsynCheck(synSet,pvDenser);
ghinf_5_C = constChk(1).bound;

       
% =========================================================================
% Evaluation/Validation
fprintf('\n')
fprintf('==================================================================\n')
fprintf('Comparison for Hinf synthesis\n')
fprintf('hinfgs:                        %5.4f (analysis %5.4f)\n',glmi, glmi_A)
fprintf('lpvsyn(pdG,r):                 %5.4f\n',ghinf_1)
fprintf('lpvsyn(pdG,u,y,const):         %5.4f (analysis %5.4f)\n',ghinf_2,ghinf_2_A)
fprintf('                               (MinDamping = %4.2f, MaxFreq = %d)\n',const(2).MinDamping, const(2).MaxFreq)
fprintf('------------------------------------------------------------------\n')
fprintf('lpvsyn(pdG,u,y,const,lyap):\n')
fprintf('  X=X(p), dX=0, Y=Y(p), dY=0:  %5.4f (synthesis, grid of %1.0f points)\n',ghinf_3,nsyn)
fprintf('                               %5.4f (checking, grid of %1.0f points)\n',ghinf_3_C,nchk)
fprintf('  aff:                         %5.4f\n',ghinf_3a)
fprintf('------------------------------------------------------------------\n')
fprintf('lpvsyn(pdG,u,y,const,lyap):\n')
fprintf('  X=X(p), dX~=0, Y=Y0, dY=0:   %5.4f (synthesis, grid of %1.0f points)\n',ghinf_4,nsyn)
fprintf('                               %5.4f (checking, grid of %1.0f points)\n',ghinf_4_C,nchk)
fprintf('  affdX:                       %5.4f\n',ghinf_4a)
fprintf('------------------------------------------------------------------\n')
fprintf('lpvsyn(pdG,u,y,const,lyap):\n')
fprintf('  X=X0, dX=0, Y=Y(p), dY~=0:   %5.4f (synthesis, grid of %1.0f points)\n',ghinf_5,nsyn)
fprintf('                               %5.4f (checking, grid of %1.0f points)\n',ghinf_5_C,nchk)
fprintf('  affdX:                       %5.4f\n',ghinf_5a)
fprintf('==================================================================\n')


% ======================================================================

% pole placement checking
idx = 1;                     % system index

Gw    = ss(pdGaw);           % augmented plant
GwGrid = ss(pdGawGrid);

K1 = ss(pdK1);
Gcl1  = lft(Gw,K1);
poles_1 = eig(Gcl1(:,:,1)); % new lpvsyn

K2 = ss(pdK2);
Gcl2  = lft(Gw,K2);
poles_2 = eig(Gcl2(:,:,1)); % new lpvsyn (constraint)

K3 = ss(pdK3);
Gcl3  = lft(GwGrid,K3);
poles_3 = eig(Gcl3(:,:,1)); % new lpvsyn with pole placement


figure
xlim = min([min(real(poles_1)),...
            min(real(poles_2)),min(real(poles_3))]); 
ylim = max([max(imag(poles_1)),...
            max(imag(poles_2)),max(imag(poles_3))]); 

% pole placement area
rad = sqrt(xlim^2+ylim^2);
d = const(2).MinDamping;
beta = atan(sqrt(1-d^2)/d);
ang = linspace(-beta,beta,50);
lim = max(-xlim,ylim);
% lim = 1e4;
ppareax = -[0 rad*cos(ang) 0]; 
ppareay = [0 rad*sin(ang) 0]; 
set(gca,'NextPlot','add',...
        'Color',0.8*[1 1 1],...
        'Xlim',[-lim 0],...
        'Ylim',lim*[-1 1])
patch(ppareax,ppareay,[1 1 1])

plot([xlim 0],[0 0],'k-');

plot(poles_1,'ks');
plot(poles_2,'bs');
plot(poles_2,'ms');

