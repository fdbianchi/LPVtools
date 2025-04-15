
% LPVtools: Nonlinear Autopilot Example (see "Modeling and Hinf control for 
%           switched linear parameter-varying missile autopilot", Lim & How, 2003)
%
%   Plant:          PWA LPV (ppss)
%   Constraint:     Hinf/Poles
%   Lyapunov Fcn:   PWA
%   Ctrller Fcn:    No specified

% fbianchi - 2023-08-09


% cleanig 
clc
clearvars
close all


% ========================================================================
% Nonlinear model

f1 =@(M,a) (0.0001*a^2 - 0.0112*abs(a) - 0.2010*(2 - M/3))*M*cos(a*pi/180);
f2 =@(M,a) (0.0152*a^2 - 1.3765*abs(a) + 3.6001*(-7 + 8*M/3))*M^2;
f3 =@(M,a) -0.0403*M*cos(a*pi/180);
f4 =@(M) -14.542*M^2;
f5 =@(M,a) (0.0001*a^2 - 0.0063*abs(a) - 0.1130*(2 - M/3))*M^2;
f6 =@(M) -0.0226*M^2;

Af =@(M,a) [f1(M,a), 1; f2(M,a), 0];
Bf =@(M,a) [f3(M,a); f4(M)];
Cf =@(M,a) [f5(M,a), 0; 0, 1];
Df =@(M,a) [f6(M); 0];

% ========================================================================
% Grid of points and parameter set
N = 3;
Points = pgrid([2, 4; -20, 20], [N; 2*N-1]);
rate = [-1.5 1.5;-100 100];
rate = [-1.5 1.5;-inf inf];
Pset = pset.Gral(Points, rate, {'M', 'a'});

% ========================================================================
% PWA LPV model
A = zeros(2,2,N);
B = zeros(2,1,N);
C = zeros(2,2,N);
D = zeros(2,1,N);
for ii = 1:size(Points,2)
    
    p = Points(:,ii);
    A(:,:,ii) = Af(p(1), p(2));
    B(:,:,ii) = Bf(p(1), p(2));
    C(:,:,ii) = Cf(p(1), p(2));
    D(:,:,ii) = Df(p(1), p(2));

end
pdG = ppss(A, B, C, D, Pset);
pdG.u = 'up';   pdG.y = {'n','q'};


% ========================================================================
% LPV design

% weights

Wref = tf(144*[-0.05 1], [1 2*0.8*12 144]);
Wref.u = 'nc';  Wref.y = 'ref';

We = tf(17.321, [1 0.0577]);
We.u = 'e';     We.y = 'er';

Wdc = tf([1 0.25], 25*[0.005 1]);
Wdc.u = 'dc';    Wdc.y = 'edc';

Wd = tf(0.01);
Wd.u = 'd';    Wd.y = 'dr';

Wn1 = tf(0.001);
Wn1.u = 'dn1';  Wn1.y = 'dn1r';

Wn2 = tf(0.001);
Wn2.u = 'dn2';  Wn2.y = 'dn2r';

Act = tf(1,[1/100 1]);
Act.u = 'dc';   Act.y = 'ua';

F = tf(1,[1/100 1]);
F.u = 'n'; F.y = 'nf';

% augmented plant + weight
sb1 = sumblk('err = nc - nf - dn1r');
sb2 = sumblk('y = q + dn2r');
sb3 = sumblk('up = ua + dr');
sb4 = sumblk('e = ref - n');
pdGaw = connect(pdG,sb1,sb2,sb3,sb4,Wref,We,Wdc,Wd,Wn1,Wn2,Act,F,...
                {'nc','d','dn1','dn2','dc'},...
                {'er','edc','err','y'});


% H-infinity constraint
const(1) = synConst.Gain({'nc','d','dn1','dn2'},{'er','edc'});
const(2) = synConst.Poles('MaxFreq',1e5);

opt = lpvsettings('solver','mosek');


% constant Lyapunov function
[pdKc,constOut] = lpvsyn(pdGaw,3:4,5,const,[],[],opt);
gCte = constOut(1).bound;

% different options for PWA Lyapunov function
[pdKpwa,constOut] = lpvsyn(pdGaw,3:4,5,const,'pwa',[],opt);
gpwa = constOut(1).bound;

[pdKpwax,constOut] = lpvsyn(pdGaw,3:4,5,const,'pwaX',[],opt);
gpwax = constOut(1).bound;

[pdKpway,constOut] = lpvsyn(pdGaw,3:4,5,const,'pwaY',[],opt);
gpway = constOut(1).bound;

[pdKpwadx,constOut] = lpvsyn(pdGaw,3:4,5,const,'pwadX',[],opt);
gpwadx = constOut(1).bound;

[pdKpwady,constOut] = lpvsyn(pdGaw,3:4,5,const,'pwadY',[],opt);
gpwady = constOut(1).bound;


% using general Lyapunov functions
lyapSet = createLyapFcn(0, 0, @(p) p(1), @(p) [1 0]);
% lyapSet = createLyapFcn(@(p) p,0);
[pdKg,constOut] = lpvsyn(pdGaw,3:4,5,const,lyapSet,[],opt);
gg = constOut(1).bound;


fprintf('\n')
fprintf('===================================================================\n')
fprintf('Comparison for Hinf synthesis\n')
fprintf('lpvsyn(pdG,u,y,const):            %5.4f\n',gCte)
fprintf('lpvsyn(pdG,u,y,const,''pwa''):      %5.4f\n',gpwa)
fprintf('lpvsyn(pdG,u,y,const,''pwaX''):     %5.4f\n',gpwax)
fprintf('lpvsyn(pdG,u,y,const,''pwaY''):     %5.4f\n',gpway)
fprintf('lpvsyn(pdG,u,y,const,''pwadX''):    %5.4f\n',gpwadx)
fprintf('lpvsyn(pdG,u,y,const,''pwadY''):    %5.4f\n',gpwady)
fprintf('lpvsyn(pdG,u,y,const,f@()):       %5.4f\n',gg)
fprintf('------------------------------------------------------------------\n')

% =======================================================================
% Checking at frozen parameter values

gridSet = pgrid([2 4;5 20],3);
Gset = ss(pdGaw, gridSet);

for ii = 1:size(gridSet,2)

    % plant
    p = gridSet(:,ii);
    A = Af(p(1), p(2));
    B = Bf(p(1), p(2));
    C = Cf(p(1), p(2));
    D = Df(p(1), p(2));  
    G = ss(A,B,C,D,'InputName','up','OutputName',{'n','q'});

    % H-infinite design
    Gaw = connect(G,sb1,sb2,sb3,sb4,Wref,We,Wdc,Wd,Wn1,Wn2,Act,F,...
                {'nc','d','dn1','dn2','dc'},...
                {'er','edc','err','y'});
    [Kaux, CL, gAux] = hinfsyn(Gaw, 2, 1, 'method', 'lmi');
    
    Khinf(:,:,ii) = Kaux;
    gHinf(ii) = gAux;
    
    Klpv = ss(pdKc,p);
    Gcl  = lft(Gaw,Klpv);
    gLPV(ii) = norm(Gcl,inf);
    
    Kpwa = ss(pdKpwa,p);
    Gcl  = lft(Gaw,Kpwa);
    gPWA(ii) = norm(Gcl,inf);

    Kpwax = ss(pdKpwax,p);
    Gcl  = lft(Gaw,Kpwax);
    gPWAx(ii) = norm(Gcl,inf);

    Kpway = ss(pdKpway,p);
    Gcl  = lft(Gaw,Kpway);
    gPWAy(ii) = norm(Gcl,inf);

    Kpwadx = ss(pdKpwadx,p);
    Gcl  = lft(Gaw,Kpwadx);
    gPWAdx(ii) = norm(Gcl,inf);

    Kpwady = ss(pdKpwady,p);
    Gcl  = lft(Gaw,Kpwady);
    gPWAdy(ii) = norm(Gcl,inf);
    
    Kg = ss(pdKg,p);
    Gcl  = lft(Gaw,Kg);
    gg(ii) = norm(Gcl,inf);
    
end

v = 1:size(gridSet,2);
plot(v,gHinf,'-x',v,gLPV,'-s',v,gPWA,'-^')
hold on
plot(v,gPWAx,'-o',v,gPWAdx,'-o')
plot(v,gPWAy,'-^',v,gPWAdy,'-^')
plot(v,gg,'.-')
xlabel('Plant indices')
ylabel('gamma')
title('Closed-loop performance at frozen parameter values')
legend('Hinf','Cte','pwa','pwaX','pwadX','pwaY','pwadY','fcn')


%% -----------------------------------------------------------------------
% Simulations

simFile = 'Autopilot_pwa_mdl.slx';
load_system(simFile);

% afin controller without poles constraints
pdK = pdKc;

set_param([simFile(1:end-4) '/LPV Controller'],'Commented','off')
set_param([simFile(1:end-4) '/LPV Controller X(p)'],'Commented','on')
set_param([simFile(1:end-4) '/sw'],'sw','1')


sim(simFile);
t = simData.time;
ar = simData.signals(1).values(:,1);
ad = simData.signals(1).values(:,2);
a = simData.signals(1).values(:,3);
u = simData.signals(2).values(:,1);
p1 = simData.signals(3).values(:,1);
p2 = simData.signals(4).values(:,1);

% color lines
clines = lines(7);
figure('Position', [680   200   800   780])

ha1 = subplot(4,1,1);
plot(t,ar,'Color',clines(1,:),'LineWidth',1); hold on
plot(t,ad,'Color',clines(2,:),'LineWidth',1);
plot(t,a,'Color',clines(3,:),'LineWidth',1);
ylabel('\alpha'); xlabel('time (s)')
ha2 = subplot(4,1,2);
plot(t,u,'Color',clines(3,:),'LineWidth',1); hold on
ylabel('u'); xlabel('time (s)')
ha3 = subplot(4,1,3);
plot(t,p1,'Color',clines(1,:),'LineWidth',1); hold on
ylabel('p1'); xlabel('time (s)')
ha4 = subplot(4,1,4);
plot(t,p2,'Color',clines(1,:),'LineWidth',1); hold on
ylabel('p2'); xlabel('time (s)')

ylim(ha1,[-20 20])
ylim(ha2,[-40 40])


% PWA controller without poles constraints
pdK = pdKpwa;

set_param([simFile(1:end-4) '/LPV Controller'],'Commented','on')
set_param([simFile(1:end-4) '/LPV Controller X(p)'],'Commented','off')
set_param([simFile(1:end-4) '/sw'],'sw','0')

sim(simFile);
t = simData.time;
a = simData.signals(1).values(:,3);
u = simData.signals(2).values(:,1);

plot(ha1,t,a,'Color',clines(4,:),'LineWidth',1);
plot(ha2,t,u,'Color',clines(4,:),'LineWidth',1);

% PWA_Y controller without poles constraints
pdK = pdKpway;

sim(simFile);
t = simData.time;
a = simData.signals(1).values(:,3);
u = simData.signals(2).values(:,1);

plot(ha1,t,a,'Color',clines(5,:),'LineWidth',1);
plot(ha2,t,u,'Color',clines(5,:),'LineWidth',1);

% PWA_dY controller without poles constraints
pdK = pdKpwady;

sim(simFile);
t = simData.time;
a = simData.signals(1).values(:,3);
u = simData.signals(2).values(:,1);

plot(ha1,t,a,'Color',clines(6,:),'LineWidth',1);
plot(ha2,t,u,'Color',clines(6,:),'LineWidth',1);


% % f@() controller without poles constraints
pdK = pdKg;

sim(simFile);
t = simData.time;
a = simData.signals(1).values(:,3);
u = simData.signals(2).values(:,1);

plot(ha1,t,a,'Color',clines(7,:),'LineWidth',1);
plot(ha2,t,u,'Color',clines(7,:),'LineWidth',1);


legend(ha1,'Ref', 'Target', 'K(Xcl=cte)','K(PWA)','K(PWAY)', 'K(PWAdY)', 'K(X(f@()))',...
    'Orientation', 'Horizontal', 'Box', 'off')


save_system(simFile);

