
% =========================================================================
% Model of a a satellite consisting of two rigid bodies joined by a 
% flexible link (the “boom”). The boom is modeled as a spring with torque 
% constant k and viscous damping f taking values in [0.09, 0.4] and
% [0.0038, 0.04], respectively
%   
% From: LMI toolbox manual, page 4-13
%
% fbianchi - 2020-08-24
%
% =========================================================================

clearvars
close all

% system parameters
J1 = 1; J2 = 0.1;

% LPV model(affine)
E0 = diag([1 1 J1 J2]);
A0 = [zeros(2) eye(2); zeros(2,4)]; 
Ak = [zeros(2,4); [-1 1;1 -1] zeros(2)];
Af = [zeros(2,4); zeros(2) [-1 1;1 -1]]; 
B  = [0 0;0 0;1 1;0 0];                     % b = [b1 b2] 
C  = [0 1 0 0; 1 0 0 0; 0 1 0 0; 0 0 0 0];  % Hinf: c1 
D  = [0 0;0 0;0 0;0 1];                     % H2: c2

% range of parameter values 
pv = pvec('box',[0.09 0.4 ; 0.0038 0.04]);

% parameter-dependent plant
pdGaug = psys(pv,[ltisys(A0,B,C,D,E0),... 
                  ltisys(Ak,0*B,0*C,0*D,0),... 
                  ltisys(Af,0*B,0*C,0*D,0)]);
pdGaug = aff2pol(pdGaug);         

% -------------------------------------------------------------------------
% Controller design: LMI toolbox          

% Pole placement region
minDecay = 0.1; theta = 0.75*pi;
minDamping = 1/sqrt(1 + tan(theta/2)^2);
regRe(1,1) = 2*minDecay;
regRe(1,4) = 1;
regRe(2:3,5:6) = [sin(theta/2) -cos(theta/2); cos(theta/2) sin(theta/2)];
regIm = zeros(3,6); regIm(1,1) = 1; regIm(2,2) = 2;
region = regRe + regIm*1i;

% robust control
[gli0,gl20,Klmi0,Gcl_lmi0] = msfsyn(pdGaug,[3 1],[0 0 1 0],region);

% pareto analysis
giSet = [0.01, 0.1, 0.2, 0.5];
N = length(giSet);
KlmiSet = zeros(1,4,N);
for ii = 1:length(giSet)
    [gi_i,g2_i,K_i,Pcl_i] = msfsyn(pdGaug,[3 1],[giSet(ii) 0 0 1],region);
    gl2(ii) = g2_i;
    KlmiSet(:,:,ii) = K_i;
    Gcl_lmi{ii} = ppss(Pcl_i,pv);
end

% -------------------------------------------------------------------------
% Controller design: LPV toolbox

% new LPV model from psys object
pdGw = ppss(pdGaug);

% constraints
const(1) = synConst.Gain(1,1);
const(2) = synConst.Poles('MinDamping',minDamping,'MinDecay',minDecay);
const(3) = synConst.GainH2(1,2:4);

% State-feedback gain scheduling 
opts = lpvsettings('solver','mosek','eigtol',1e-6);
[Klpv,constR] = lpvsyn(pdGw,0,2,const(1:2),[],0,opts);
gpi0 = constR(1).bound;

% pareto analysis
KlpvSet = zeros(1,4,N);
% pdGw_ext = ppss(pdGw.A,B,[C;eye(4)],[D;zeros(4,2)],pdGw.parset);
for ii = 1:length(giSet)
    const(1) = synConst.Gain(1,1,'bound',giSet(ii));
    [Klpv_i,constR] = lpvsyn(pdGw,0,2,const([1 2 3]),[],0,opts);
    gp2(ii) = constR(3).bound;
    KlpvSet(:,:,ii) = Klpv_i.D;
    
    Gcl_lpv{ii} = lft(pdGw,Klpv_i);
end

fprintf('\n')
fprintf('-------------------------\n')
fprintf('Robust: %2.3f (lmitool)\n',gli0);
fprintf('Robust: %2.3f (lpvtool)\n',gpi0);
fprintf('-------------------------\n\n')


figure
plot(giSet,gl2,'s-')
hold on
plot(giSet,gp2,'s-')
ylabel('H2 performance')
xlabel('H-infinity performance')


%% -------------------------------------------------------------------------
% Pole placement verifications

id_sys = 2;

eig_lmi = eig(Gcl_lmi{id_sys});
eig_lpv = eig(Gcl_lpv{id_sys});

xlim = min(min(real([squeeze(eig_lmi),squeeze(eig_lpv)])));
ylim = min(min(imag([squeeze(eig_lmi),squeeze(eig_lpv)])));

d = const(2).MinDamping;
beta = atan(sqrt(1-d^2)/d);
lim = max(-xlim,ylim);
figure
set(gca,'NextPlot','add',...
        'Xlim',[-lim*cos(beta) 0],...
        'Ylim',lim*sin(beta)*[-1 1])
ppareax = [0 -lim*cos(beta) 0 0 -lim*cos(beta) 0]; 
ppareay = [0 lim*sin(beta) lim*sin(beta) -lim*sin(beta) -lim*sin(beta) 0]; 
patch(ppareax,ppareay,0.8*[1 1 1],'LineStyle','none')
ppareax = [0 -1 -1 0]*const(2).MinDecay; 
ppareay = [1 1 -1 -1]*lim*sin(beta); 
patch(ppareax,ppareay,0.8*[1 1 1],'LineStyle','none')
for ii = 1:4
    h = plot(real(eig_lmi(:,:,ii)),imag(eig_lmi(:,:,ii)),'rx',...
             real(eig_lpv(:,:,ii)),imag(eig_lpv(:,:,ii)),'bx');
end
set(gca, 'Layer', 'top')
legend(h,'LMItool','LPVtool')
title('Eigenvalues of the closed-loop plant, for g=0.1')

 
figure
impulse(Gcl_lmi{id_sys}(1,1),'r',Gcl_lpv{id_sys}(1,1),'b')
legend('LMItool','LPVtool')
title('Impulse response w -> x2, for g=0.1')

