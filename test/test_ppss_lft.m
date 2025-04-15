

clear all
clc
close all

% ======================================================================
% pdG1: LPV & pdG2: LTI
ns1 = 3; nz1 = 1; nw1 = 1; nv1 = 4;
ns2 = 2; nz2 = 0; nw2 = 0;
nu  = 2; ny  = 1;

%%  
% channel u parameter dependent

A = rand(ns1,ns1,nv1);
B = rand(ns1,nw1+nu,nv1);
C = rand(nz1+ny,ns1,nv1);
C(nz1+1:end,:,2:end) = repmat(C(nz1+1:end,:,1),1,1,nv1-1);
D = rand(nz1+ny,nw1+nu,nv1);
D(nz1+1:end,:,2:end) = repmat(D(nz1+1:end,:,1),1,1,nv1-1);
D(nz1+1:end,nw1+1:end,:) = 0;
pdG1 = ppss(A,B,C,D,pset.Box([0 2;2 6]));

G2 = rss(ns2,nz2+nu,nw2+ny);

pdGo = lft(pdG1,G2)

Go = ss(pdGo);
Gs = lft(ss(pdG1),G2);

figure
bodemag(Go,'b-',Gs,'r^')
legend('LPV lft','LTi lft')
title('Case u channel parameter dependent')


%% 
% channel y parameter dependent
A = rand(ns1,ns1,nv1);
B = rand(ns1,nw1+nu,nv1);
B(:,nw1+1:end,2:end) = repmat(B(:,nw1+1:end,1),1,1,nv1-1);
C = rand(nz1+ny,ns1,nv1);
D = rand(nz1+ny,nw1+nu,nv1);
D(:,nw1+1:end,2:end) = repmat(D(:,nw1+1:end,1),1,1,nv1-1);
D(nz1+1:end,nw1+1:end,:) = 0;
pdG1 = ppss(A,B,C,D,pset.Box([0 2;2 6]));

G2 = rss(ns2,nz2+nu,nw2+ny);

pdGo = lft(pdG1,G2)

Go = ss(pdGo);
Gs = lft(ss(pdG1),G2);

figure
bodemag(Go,'b-',Gs,'r^')
legend('LPV lft','LTi lft')
title('Case y channel parameter dependent')

%%
% ======================================================================
% pdG1: LPV & pdG2: LPV

ns1 = 3; nz1 = 1; nw1 = 1; nv1 = 4;
ns2 = 2; nz2 = 0; nw2 = 0; nv2 = 4;
nu  = 2; ny  = 1;

A = rand(ns1,ns1,nv1);
B = rand(ns1,nw1+nu,nv1);
B(:,nw1+1:end,2:end) = repmat(B(:,nw1+1:end,1),1,1,nv1-1);
C = rand(nz1+ny,ns1,nv1);
C(nz1+1:end,:,2:end) = repmat(C(nz1+1:end,:,1),1,1,nv1-1);
D = rand(nz1+ny,nw1+nu,nv1);
D(:,nw1+1:end,2:end) = repmat(D(:,nw1+1:end,1),1,1,nv1-1);
D(nz1+1:end,:,2:end) = repmat(D(nz1+1:end,:,1),1,1,nv1-1);
D(nz1+1:end,nw1+1:end,:) = 0;
pdG1 = ppss(A,B,C,D,pset.Box([0 2;2 6]));

Ak = rand(ns2,ns2,nv2);
Bk = rand(ns2,ny,nv2);
Ck = rand(nu,ns2,nv2);
Dk = rand(nu,ny,nv2);
pdG2 = ppss(Ak,Bk,Ck,Dk,pset.Box([0 2;2 6]));

pdGo = lft(pdG1,pdG2)

Go = ss(pdGo);
Gs = lft(ss(pdG1),ss(pdG2));

figure
bodemag(Go,'b-',Gs,'r^')
legend('LPV lft','LTi lft')
title('Both plants LPV')

