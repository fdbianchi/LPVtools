

clear all
clc
close all

% ======================================================================
% pdG1: LPV & pdG2: LTI
ns1 = 3; nz1 = 1; nw1 = 1; nv1 = 3;
ns2 = 2; nz2 = 0; nw2 = 0;
nu  = 2; ny  = 1;

%%  
% channel u parameter dependent

A = rand(ns1,ns1,nv1);
B = rand(ns1,nz1+nu,nv1);
C = rand(nw1+ny,ns1,nv1);
C(nz1+1:end,:,:) = 0;
D = rand(nz1+ny,nw1+nu,nv1);
D(nz1+1:end,:,:) = 0;
pv = pset.Box([0 2]);
fcn =@(p) [p;cos(p)];
pdG1 = pgss(A,B,C,D,pv,fcn);

G2 = rss(ns2,nz2+nu,nw2+ny);
% G2.D = zeros(nz2+nu,nw2+ny);

pdGo = lft(pdG1,G2);

Go = ss(pdGo);
Gs = lft(ss(pdG1),G2);

figure
bodemag(Go,'b-',Gs,'r^-')
legend('LPV lft','LTi lft')
title('Case u channel parameter dependent')


%% 
% channel y parameter dependent

A = rand(ns1,ns1,nv1);
B = rand(ns1,nz1+nu,nv1);
B(:,nw1+1:end,:) = 0;
C = rand(nw1+ny,ns1,nv1);
D = rand(nz1+ny,nw1+nu,nv1);
D(:,nw1+1:end,:) = 0;

pdG1 = pgss(A,B,C,D,pv,fcn);

G2 = rss(ns2,nz2+nu,nw2+ny);
% G2.D = zeros(nz2+nu,nw2+ny);

pdGo = lft(pdG1,G2);

Go = ss(pdGo);
Gs = lft(ss(pdG1),G2);

figure
bodemag(Go,'b-',Gs,'r^-')
legend('LPV lft','LTi lft')
title('Case y channel parameter dependent')


%%
% ======================================================================
% pdG1: LPV & pdG2: LPV

ns1 = 3; nz1 = 1; nw1 = 1; nv1 = 3;
ns2 = 2; nz2 = 0; nw2 = 0; nv2 = 3;
nu  = 2; ny  = 1;

A = rand(ns1,ns1,nv1);
B = rand(ns1,nz1+nu,nv1);
B(:,nw1+1:end,:) = 0;
C = rand(nw1+ny,ns1,nv1);
C(nz1+1:end,:,:) = 0;
D = rand(nz1+ny,nw1+nu,nv1);
D(nz1+1:end,:,:) = 0;
D(:,nw1+1:end,:) = 0;
pdG1 = pgss(A,B,C,D,pv,fcn);

Ak = rand(ns2,ns2,nv2);
Bk = rand(ns2,ny,nv2);
Ck = rand(nu,ns2,nv2);
Dk = rand(nu,ny,nv2);
pdG2 = pgss(Ak,Bk,Ck,Dk,pv,fcn);

pdGo = lft(pdG1,pdG2);

Go = ss(pdGo);
Gs = lft(ss(pdG1),ss(pdG2));

figure
bodemag(Go,'b-',Gs,'r^-')
legend('LPV lft','LTi lft')
title('Both plants LPV')

