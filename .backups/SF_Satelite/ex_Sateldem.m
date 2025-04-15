
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

% system parameters
J1 = 1; J2 = 0.1;

% LPV model(affine)
E0 = diag([1 1 J1 J2]);
A0 = [zeros(2) eye(2); zeros(2,4)]; 
Ak = [zeros(2,4); [-1 1;1 -1] zeros(2)];
Af = [zeros(2,4); zeros(2) [-1 1;1 -1]]; 
B  = [0 0;0 0;1 1;0 0];                     % b = [b1 b2] 
C  = [0 1 0 0; 1 0 0 0; 0 1 0 0; 0 0 0 0];  % H2: c1 
D  = [0 0;0 0;0 0;0 1];                     % Hinf: c2

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
[gi0,g20,pdKlmil] = msfsyn(pdGaug,[3 1],[0 0 1 0],region);

% -------------------------------------------------------------------------
% Controller design: LPV toolbox

% new LPV model from psys object
pdGw = ppss(pdGaug);

% constraints
const(1) = synConst.Gain(1,4);
const(2) = synConst.Poles('MinDamping',minDamping,'MinDecay',minDecay);

% State-feedback gain scheduling 
[pdKgs,constR] = lpvsyn(pdGw,0,2,const);
gi1 = constR(1).bound;

fprintf('\n')
fprintf('-------------------------\n')
fprintf('Robust: %2.3f\n',gi0);
fprintf('GS:     %2.3f\n',gi1);
fprintf('-------------------------\n\n')

return

% -------------------------------------------------------------------------
% Pole placement verifications
eig_r = [];
eig_gs = [];
eig_gspp = [];
[typ,nv] = psinfo(pdKgs);
for ii=1:4
   [A,B,C,D] = ltiss(psinfo(pdGaug,'sys',ii));
   eig_r = [eig_r; eig(A+B(:,2)*pdKr)];
   Gcl_r(:,:,ii) = ss(A+B(:,2)*pdKr,B(:,1),C+D(:,2)*pdKr,D(:,1));
   
   [xx,xx,xx,Kgs] = ltiss(psinfo(pdKgs,'sys',min(ii,nv)));
   eig_gs = [eig_gspp; eig(A+B(:,2)*Kgs)];
   Gcl_gs(:,:,ii) = ss(A+B(:,2)*Kgs,B(:,1),C+D(:,2)*Kgs,D(:,1));
   
end
plot(real(eig_r),imag(eig_r),'rx',...
     real(eig_gs),imag(eig_gs),'bx');
set(gca,'xlim',[-3 0],'ylim',[-3 3]);
legend('Robust','GS')
hold on
% cone 
x=[-9:1:0]; plot(x,x/tan(pi/8),'m'); plot(x,-x/tan(pi/8),'m'); 
% half-plane
y=[-10:1:10]; plot(-0.1*ones(size(y)),y,'m'); hold off
 
figure
impulse(Gcl_r(1,:),'r',Gcl_gs(1,:),'b')
legend('Robust','GS')
