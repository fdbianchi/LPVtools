
% LPVtools: Autopilot Example (see "LMI control toolbox manual", pp 7-10)
%
%   Plant:          affine LPV (psys/pass)
%   Constraint:     Hinf/Poles
%   Lyapunov Fcn:   PWA
%   Ctrller Fcn:    No specified
%
% PWA version

% fbianchi - 2021-04-01


% cleanig 
clc
clearvars
close all

% ------------------------------------------------------------------------
% Modelling

% Operating range
Zmin =  0.5;   Zmax =   4;
Mmin =  0.0;   Mmax = 106;  

% system matrices
a0 = [ 0 1; 0 0]; a1 = [-1 0; 0 0]; a2 = [ 0 0;-1 0];
b0 = [0;1];       b1 = [0;0];       b2 = [0;0];
c0 = [-1 0;0 1];  c1 = [ 0 0; 0 0]; c2 = [ 0 0; 0 0];
d0 = [0;0];       d1 = [0;0];       d2 = [0;0];

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
W1 = tf(nf1,df1);
W2 = tf(nf2,df2);
Wout  = balreal(append(W1,W2,1,1));   
Wout.y = {'et','ut','e','y(2)'};

% augmented plant
er  = sumblk('e = r - y(1)');
pdGau = connect(pdG,er,{'r','u'},{'e','u','e','y(2)'});

% augmented plant + weight
pdGaw = Wout*pdGau;     % affine model
pdGawp = ppss(pdGaw);   % PWA model
pdGawp.parset = pset.Gral(Pset.points,[-10 10],{'Z','M'});

% H-infinity constraint
const(1) = synConst.Gain('r',{'et','ut'});
const(2) = synConst.Poles('MaxFreq',1e5);


opt = lpvsettings('solver','sedumi');

% constant Lyapunov function
% affine model
[pdK1,constOut] = lpvsyn(pdGaw,[2 1]);
ghinf_1 = constOut(1).bound;

[pdK2,constOut] = lpvsyn(pdGaw,{'e','y(2)'},'u',const,[],[],opt);
ghinf_2 = constOut(1).bound;
% pwa model
[pdK2p,constOut] = lpvsyn(pdGawp,{'e','y(2)'},'u',const,[],[],opt);
ghinf_2p = constOut(1).bound;

% X & Y parameter dependent, dX=0, dY=0
lyapSet = createLyapFcn(@(p) p,0,@(p) p,0);
[pdK3,constOut]  = lpvsyn(pdGaw,{'e','y(2)'},'u',const,lyapSet,[],opt);
ghinf_3 = constOut(1).bound;

% affine model
[pdK3a,constOut] = lpvsyn(pdGaw,{'e','y(2)'},'u',const,'aff',[],opt);
ghinf_3a = constOut(1).bound;

% pwa model
[pdK3p,constOut] = lpvsyn(pdGawp,{'e','y(2)'},'u',const,'pwa',[],opt);
ghinf_3p = constOut(1).bound;


fprintf('\n')
fprintf('===================================================================\n')
fprintf('Comparison for Hinf synthesis\n')
fprintf('lpvsyn(pdG,r):                    %5.4f\n',ghinf_1)
fprintf('lpvsyn(pdG,u,y,const):            %5.4f\n',ghinf_2)
fprintf('lpvsyn(pdG,u,y,const):            %5.4f (converting pdG into PWA) \n',ghinf_2p)
fprintf('lpvsyn(pdG,u,y,const,f@()):       %5.4f\n',ghinf_3)
fprintf('lpvsyn(pdG,u,y,const,''aff''):      %5.4f\n',ghinf_3a)
fprintf('lpvsyn(pdG,u,y,const,''pwa''):      %5.4f\n',ghinf_3p)
fprintf('------------------------------------------------------------------\n')


%% -----------------------------------------------------------------------
% Simulations

simFile = 'Autopilot_mdl.slx';
load_system(simFile)

% afin controller without poles constraints
pdK = pdK1;

set_param([simFile(1:end-4) '/LPV Controller'],'Commented','off')
set_param([simFile(1:end-4) '/LPV Controller X(p)'],'Commented','on')
set_param([simFile(1:end-4) '/sw'],'sw','1')


sim(simFile);
t = simData.time;
a = simData.signals(1).values(:,1);
q = simData.signals(2).values(:,1);
p1 = simData.signals(3).values(:,1);
p2 = simData.signals(3).values(:,2);

ha1 = subplot(4,1,1);
plot(t,a); hold on
ylabel('\alpha'); xlabel('time (s)')
ha2 = subplot(4,1,2);
plot(t,q); hold on
ylabel('q'); xlabel('time (s)')
ha3 = subplot(4,1,3);
plot(t,p1); hold on
ylabel('p1'); xlabel('time (s)')
ha4 = subplot(4,1,4);
plot(t,p2); hold on
ylabel('p2'); xlabel('time (s)')

ylim(ha1,[0 1.5])
ylim(ha2,[-15 5])

% afin controller with poles constraints
pdK = pdK2;

sim(simFile);
t = simData.time;
a = simData.signals(1).values(:,1);
q = simData.signals(2).values(:,1);

plot(ha1,t,a);
plot(ha2,t,q);

% PWA controller without poles constraints
pdK = pdK2p;

sim(simFile);
t = simData.time;
a = simData.signals(1).values(:,1);
q = simData.signals(2).values(:,1);

plot(ha1,t,a);
plot(ha2,t,q);

% f@() controller without poles constraints
pdK = pdK3;

set_param([simFile(1:end-4) '/LPV Controller'],'Commented','on')
set_param([simFile(1:end-4) '/LPV Controller X(p)'],'Commented','off')
set_param([simFile(1:end-4) '/sw'],'sw','0')

sim(simFile);
t = simData.time;
a = simData.signals(1).values(:,1);
q = simData.signals(2).values(:,1);

plot(ha1,t,a);
plot(ha2,t,q);

% affine controller without poles constraints
pdK = pdK3a;

sim(simFile);
t = simData.time;
a = simData.signals(1).values(:,1);
q = simData.signals(2).values(:,1);

plot(ha1,t,a);
plot(ha2,t,q);

% pwa controller without poles constraints
pdK = pdK3p;

sim(simFile);
t = simData.time;
a = simData.signals(1).values(:,1);
q = simData.signals(2).values(:,1);

plot(ha1,t,a);

legend(ha1,'X cte','X cte + PP','PWA cte + PP',...
       'X f@() + PP', 'X affine + PP', 'X PWA + PP')

plot(ha2,t,q);

save_system(simFile)
close_system(simFile)
