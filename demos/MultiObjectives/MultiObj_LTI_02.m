
% LPVtools: Example of multi-objective synthesis (LTI)
%
% Scherer, C. W., Gahinet, P., & Chilali, M. (1997). Multiobjective output-
% feedback control via LMI optimization. IEEE Transactions on Automatic 
% Control, 42(7), 896–911. doi:10.1109/9.599969 (Example 7)
%
%   Plant:          LTI (ppss)
%   Constraint:     Hinf/H2/Poles
%   Lyapunov Fcn:   Constant
%   Ctrller Fcn:    No specified

% fbianchi - 2021-06-15

% cleanig 
clc
clearvars
close all

% ========================================================================
% 1) Modelling:

% Plant
A = [0        0       1        0;...
     0        0       0        1;...
    -0.1010  -0.1681 -0.04564 -0.01075;... 
     0.06082 -2.1407 -0.05578 -0.1273];
B = [0        0       0; 
     0        0       0; 
     0.1179   0.1441  0.1476; 
     0.1441   1.7057 -0.7557];
C = [eye(2) zeros(2)];
D = zeros(2,3);
G = ss(A,B,C,D);
G.u = {'F','M','Fu'};
G.y = {'Y','phi'};

figure
subplot(2,2,1)
sigma(G(1,1))
title('Open loop F -> Y')
subplot(2,2,2)
sigma(G(2,1))
title('Open loop F -> Phi')
subplot(2,2,3)
sigma(G(1,2))
title('Open loop M -> Y')
subplot(2,2,4)
sigma(G(2,2))
title('Open loop M -> Phi')

% actuator
Act = tf(1,[0.7 1]);
Act.u = 'u';
Act.y = 'Fu';

Wy = tf([1 0.1],[1 0]);
Wy.u = 'Y';
Wy.y = 'Yt';

Wphi = tf(0.1,1);
Wphi.u = 'phi';
Wphi.y = 'phit';

Wu = tf(0.1,1);
Wu.u = 'u';
Wu.y = 'ut';

% augmented plant
Gaug = connect(G,Act,Wy,Wphi,Wu,{'F','M','u'},{'Yt','phit','ut','Yt','phi'});


% ========================================================================
% Design 1: Hinfinity

[Klti,~,g_o,~] = hinfsyn(Gaug,2,1,'method','lmi');

const(1) = synConst.Gain([1 2],[1 2 3]);
opts = lpvsettings('solver','sedumi');
[Klpv,constOut] = lpvsyn(Gaug,[4 5],3,const,[],[],opts);
g_n = constOut(1).bound;

Gcl = lft(Gaug,ss(Klpv));
[constA,objA] = lpvanalysis(Gcl,const);
g_nA = constA(1).bound;

Gaug_check = connect(G,Act,Wy,Wphi,{'F','M','u'},{'Y','phit','u','Yt','phi'});
Gcl_lti = lft(Gaug_check,Klti);
n_Hinf = norm(Gcl_lti,inf); 
Gcl_lpv = lft(Gaug_check,ss(Klpv));
n_LPV = norm(Gcl_lpv,inf); 

fprintf('\n--------------------------------------------------------\n')
fprintf('Comparison for Hinf synthesis\n')
fprintf('h2hinfsyn:  %6.4f (norm = %6.4f)\n',g_o,n_Hinf)
fprintf('lpvsyn:     %6.4f (norm = %6.4f, analysis %6.4f)\n',g_n,n_LPV,g_nA)
fprintf('--------------------------------------------------------\n')

figure
subplot(2,2,1)
sigma(Gcl_lti(1,1),Gcl_lpv(1,1),{0.01 100})
title('CL model, design #1, F -> Y')
subplot(2,2,2)
sigma(Gcl_lti(2,1),Gcl_lpv(2,1),{0.01 100})
title('CL model, design #1, F -> Phi')
subplot(2,2,3)
sigma(Gcl_lti(1,2),Gcl_lpv(1,2),{0.01 100})
title('CL model, design #1, M -> Y')
subplot(2,2,4)
sigma(Gcl_lti(2,2),Gcl_lpv(2,2),{0.01 100})
title('CL model, design #1, M -> Phi')
legend('hinfsyn','lvpsyn')


% ========================================================================
% Design 2: Hinfinity + H2

bndhinf = 0.5;
Gaug = connect(G,Act,Wy,Wphi,{'F','M','u'},{'Yt','phit','u','Yt','phi'});

[Klti,~,NORMZ,~] = h2hinfsyn(Gaug,2,1,1,[0 1],'hinfmax',bndhinf);
g_o = NORMZ(2);

const(1) = synConst.GainH2([1 2],3);
const(2) = synConst.Gain([1 2],[1 2],'bound',bndhinf);
[Klpv,constOut] = lpvsyn(Gaug,[4 5],3,const,[],[],opts);
g_n = constOut(1).bound;

Gcl = lft(Gaug,ss(Klpv));
[constA,objA] = lpvanalysis(Gcl,const);
g_nA = constA(1).bound;

Gcl_lti = lft(Gaug_check,Klti);
n2_lti = norm(Gcl_lti(3,1:2),2); 
nh_lti = norm(Gcl_lti(1:2,1:2),inf); 
Gcl_lpv = lft(Gaug_check,Klpv);
n2_lpv = norm(Gcl_lpv(3,1:2),2); 
nh_lpv = norm(Gcl_lpv(1:2,1:2),inf); 

fprintf('\n--------------------------------------------------------\n')
fprintf('Comparison for min H2 s.t. Hinf < %3.2f\n',bndhinf)
fprintf('h2hinfsyn:  %6.4f (norm 2 = %6.4f, norm_inf = %6.4f)\n',...
    g_o,n2_lti,nh_lti)
fprintf('lpvsyn:     %6.4f (norm 2 = %6.4f, norm_inf = %6.4f)\n',...
    g_n,n2_lpv,nh_lpv)
fprintf('                   (analysis %6.4f)\n',g_nA)
fprintf('--------------------------------------------------------\n')


figure
subplot(2,2,1)
sigma(Gcl_lti(1,1),Gcl_lpv(1,1),{0.01 100})
title('CL model, design #2, F -> Y')
subplot(2,2,2)
sigma(Gcl_lti(2,1),Gcl_lpv(2,1),{0.01 100})
title('CL model, design #2, F -> Phi')
subplot(2,2,3)
sigma(Gcl_lti(1,2),Gcl_lpv(1,2),{0.01 100})
title('CL model, design #2, M -> Y')
subplot(2,2,4)
sigma(Gcl_lti(2,2),Gcl_lpv(2,2),{0.01 100})
title('CL model, design #2, M -> Phi')
legend('hinfsyn','lvpsyn')


% ========================================================================
% Design 3: Hinfinity + H2 + pole placement

% lmi region
theta = 6*pi/7;
L = zeros(2); L(1,1) = L(1,1)+2i;
M = [sin(theta/2) -cos(theta/2); cos(theta/2) sin(theta/2)];
reg = [L M];
bndhinf = 0.5;

[Klti,CL,NORMZ,INFO] = h2hinfsyn(Gaug,2,1,1,[0 1],'hinfmax',bndhinf,'region',reg);
g_o = NORMZ(2);

const(1) = synConst.GainH2([1 2],3);
const(2) = synConst.Gain([1 2],[1 2],'bound',bndhinf);
const(3) = synConst.Poles('MinDamping',cos(theta/2));
[Klpv,constOut] = lpvsyn(Gaug,[4 5],3,const,[],[],opts);
g_n = constOut(1).bound;

Gcl = lft(Gaug,ss(Klpv));
[constA,objA] = lpvanalysis(Gcl,const);
g_nA = constA(1).bound;

Gcl_lti = lft(Gaug_check,Klti);
n2_lti = norm(Gcl_lti(3,1:2),2); 
nh_lti = norm(Gcl_lti(1:2,1:2),inf); 
Gcl_lpv = lft(Gaug_check,Klpv);
n2_lpv = norm(Gcl_lpv(3,1:2),2); 
nh_lpv = norm(Gcl_lpv(1:2,1:2),inf); 

fprintf('\n--------------------------------------------------------\n')
fprintf('Comparison for min H2 s.t. Hinf < %3.2f + PP\n',bndhinf)
fprintf('h2hinfsyn:  %6.4f (norm 2 = %6.4f, norm_inf = %6.4f)\n',...
    g_o,n2_lti,nh_lti)
fprintf('lpvsyn:     %6.4f (norm 2 = %6.4f, norm_inf = %6.4f)\n',...
    g_n,n2_lpv,nh_lpv)
fprintf('                   (analysis %6.4f)\n',g_nA)
fprintf('--------------------------------------------------------\n')

figure
subplot(2,2,1)
sigma(Gcl_lti(1,1),Gcl_lpv(1,1),{0.01 100})
title('CL model, design #3, F -> Y')
subplot(2,2,2)
sigma(Gcl_lti(2,1),Gcl_lpv(2,1),{0.01 100})
title('CL model, design #3, F -> Phi')
subplot(2,2,3)
sigma(Gcl_lti(1,2),Gcl_lpv(1,2),{0.01 100})
title('CL model, design #3, M -> Y')
subplot(2,2,4)
sigma(Gcl_lti(2,2),Gcl_lpv(2,2),{0.01 100})
title('CL model, design #3, M -> Phi')
legend('hinfsyn','lvpsyn')

% poles
poles_lti = eig(Gcl_lti);
poles_lpv = eig(Gcl_lpv);

figure
xlim = min([min(real(poles_lti)),min(real(poles_lpv))]); 
ylim = max([max(imag(poles_lti)),max(imag(poles_lpv))]); 

% pole placement area
rad = sqrt(xlim^2 + ylim^2);
d = const(3).MinDamping;
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

plot(poles_lti,'ks');
plot(poles_lpv,'rs','MarkerSize',11);

