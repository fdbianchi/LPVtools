
% Aproximation of LPV models
%
% Academic example: 2nd order system in controllble form

% fbianchi - 2020-12-10


% cleaning
clearvars; 
clc; 
close all


% ========================================================================
% Full LPV model 

% damping coefficient
ep0 =  0.7;
ep1 = -0.6;
ep2 =  0.05;
% natural frequency
wn0 = 1;    
wn1 = 0.2;  
% polynomial coefficients
a1 = 2*[ep0*wn0 ep1*wn0+ep0*wn1 ep2*wn0 ep2*wn1 ep1*wn1];
a2 = [wn0^2 2*wn0*wn1 0 0 wn1^2];

% gains
b0 = 1;
c0 = 1;

% parameter set
n1 = 5; n2 = 4;
range = [ 0 1;
          0 1];
pv = pset.Grid(range,[n1 n2],[]);

% parameter functions
fcnpar =@(p) [p(1);
              p(2);
              p(1).*p(2);
              p(1).^2];

% system matrix   
A(:,:,1) = [0 1;-a2(1) -a1(1)];
for ii = 2:5
    A(:,:,ii) = [0 0;-a2(ii) -a1(ii)];
%     svd(A(:,:,ii))
end
B = [0;1]*b0;
C = [1 0]*c0;
D = zeros(1,1);

% full LPV model
pdG = pgss(A,B,C,D,pv,fcnpar);
pdG.u = 'u';
pdG.y = 'y';

% -----------------------------------------------------------------------
% partial models

% without p1*p2
ind1 = [1 2 3 5];
% fcnpar1 =@(p) E(ind,:)*fcnpar(p);
fcnpar1 = subfunc(fcnpar,ind1(2:end)-1);
pdG1 = pgss(A(:,:,ind1),B,C,D,pv,fcnpar1);
pdG1.u = 'u'; pdG1.y = 'y';

% without p2 & p1*p2
ind2 = [1 3 5];
% fcnpar2 =@(p) E(ind,:)*fcnpar(p);
fcnpar2 = subfunc(fcnpar,ind2(2:end)-1);
pdG2 = pgss(A(:,:,ind2),B,C,D,pv,fcnpar2);
pdG2.u = 'u'; pdG2.y = 'y';

% affine
range3 = [0 1; 0 1; 0 1; 0 1];
pv3 = pset.Box(range3);
pdG3 = pass(A,B,C,D,pv3);
pdG3.u = 'u'; pdG3.y = 'y';

% bodes
figure
p = [0.5 0.75 1;
     0.5 0.75 1];
G  = ss(pdG,p);
G1 = ss(pdG1,p); 
G2 = ss(pdG2,p); 
bodemag(G,G1,G2,{0.1 10})
title('Frequency response models and approximations')
legend('G: Complete model','G1: Model with terms 1, 2, 3 & 5','G2: Model with terms 1, 3 & 5',...
    'Location','Southwest'); 
legend boxoff
% print('-depsc2','-loose','fig_Bode.eps')


% -----------------------------------------------------------------------
% LPV control

s = tf('s');
M = 1/s;%tf(1);
M.u = 'e'; M.y = 'ei';
sb = sumblk('e = r - y');
pdGau = connect(pdG,M,sb,{'r','u'},{'ei','u','ei'});

% weigths
% wc = 10; q = 100; p = wc/q; z = wc*q; k = q;
% W1 = k*(s/z+1)/(s/p+1);
W1 = 3;
wc = 20; q = 10; z = wc/q; p = wc*q; k = 1/q;
W2 = 0.1*k*(s/z+1)/(s/p+1);
% W1 = 5;
% wc = 50; q = 10; z = wc/q; p = wc*q; k = 1/q;
% W2 = 0.5*k*(s/z+1)/(s/p+1);
wout = append(W1,W2,1);

% augmented plant + weigths
pdGaw = wout*pdGau;

% -----------------------------------------------------------------------
% SVDs
A  = pdGaw.A; 
nv = size(A,3);
np = size(pdGaw.parset.points,2);
fprintf('\n')
fprintf('-----------------------------------------------------\n')
fprintf('SVD for each A matrix in A(p)\n')
U = 1:nv;
for ii = U
    if ii == 1
        sg = max(svd(A(:,:,ii)));
    else
        sg = max(svd(A(:,:,ii)));
    end        
    fprintf('svd(A_%1.0f): %6.4f\n',ii,sg)
end

% Sets of terms included in the synthesis
Is = {[1 2 3 4 5];
      [1 2 3 5];
      [1 2 4 5];
      [1];
      };

fprintf('\n-----------------------------------------------------------\n')
fprintf('Testing removing several terms in the matrix A\n')
fprintf('\tI indices of included terms, J indices of removed terms\n')  
fprintf('\tAb: matrix with terms in I, At: matrix with terms in J\n')  
for kk = 1:size(Is,1)
    
    fprintf('\n')
    I = Is{kk}; J = setdiff(U,I);
    auxStr1 = repmat('%1.0f,',1,size(I,2));
    auxStr2 = repmat('%1.0f,',1,size(J,2));
    auxStr3 = ['Test for I={' auxStr1(1:end-1) '}, J={' auxStr2(1:end-1) '}\n'];
    if isempty(J)
        fprintf(auxStr3,I);
    else
        fprintf(auxStr3,I,J);
    end
    sgb = zeros(1,np); sgt = sgb;
    for ii = 1:np
        p  = pv.points(:,ii);
        fpar = pdG.parfcn(p);
        
        Ab = 0;
        for jj = I
            if jj == 1
                Ab = Ab + A(:,:,jj);
            else
                Ab = Ab + A(:,:,jj)*fpar(jj-1);
            end
        end
        At = 0;
        for jj = J
            if jj == 1
                At = At + A(:,:,jj);
            else
                At = At + A(:,:,jj)*fpar(jj-1);
            end
        end
        
        sgb(ii) = max(svd(Ab));     sgt(ii) = max(svd(At));
        
        %     fprintf('svd(Ab): %6.4f - svd(At): %6.4f - p = [%3.2f %3.2f]\n',sgb(ii),st(ii),p)
    end
    % fprintf('Maximum\n')
    fprintf('\t max svd(Ab): %6.4f, max svd(At): %6.4f\n',max(sgb),max(sgt))
    
end
fprintf('-----------------------------------------------------------\n\n')


% design objectives
const(1)   = synConst.Gain(1,1:2);
const(2)   = synConst.Poles('MaxFreq',5000,'MinDecay',0,'MinDamping',0.1);

% ------------------------------------------------------------------------
% Design -> full parameters

[pdKf,constOut] = lpvsyn(pdGaw,3,2,const);
g_f0 = constOut(1).bound;

% full parameters (controller explicit) to check
idx0 = Is{1};
[pdK_th0,constOut] = lpvsyn(pdGaw,3,2,const,[],idx0);
g_th0 = constOut(1).bound;
pdGcl_th0 = lft(pdGaw,pdK_th0);
constOutA = lpvanalysis(pdGcl_th0,const);
g_th0r = constOutA(1).bound;

% ------------------------------------------------------------------------
% Design -> removing: p(1)*p(2)

idx1  = Is{2}; q = 1e-4;

% Procedure 2
[pdK_p21,constOut] = lpvsyn(pdGaw,3,2,const,[],idx1);
g_p21  = constOut(1).bound;
Gaw = ss(pdGaw,pv.points);
K_p21 = ss(pdK_p21,pv.points);
Gcl = lft(Gaw,K_p21);
pdGcl = ppss(Gcl,pv);
constOutA = lpvanalysis(pdGcl,const);
g_p21r = constOutA(1).bound;

% Procedure 1
ctrlfcn.idxfcn  = idx1;
ctrlfcn.XArYbnd = q;
[pdK_p11,constOut] = lpvsyn(pdGaw,3,2,const,[],ctrlfcn);
g_p11  = constOut(1).bound;
Gaw = ss(pdGaw,pv.points);
K_p11 = ss(pdK_p11,pv.points);
Gcl = lft(Gaw,K_p11);
pdGcl = ppss(Gcl,pv);
constOutA = lpvanalysis(pdGcl,const);
g_p11r = constOutA(1).bound;

% ------------------------------------------------------------------------
% Design -> removing: p(2)

idx2 = Is{3};

% Procedure 2
[pdK_p22,constOut] = lpvsyn(pdGaw,3,2,const,[],idx2);
g_p22  = constOut(1).bound;
Gaw = ss(pdGaw,pv.points);
K_p22 = ss(pdK_p22,pv.points);
Gcl = lft(Gaw,K_p22);
pdGcl = ppss(Gcl,pv);
constOutA = lpvanalysis(pdGcl,const);
g_p22r = constOutA(1).bound;

% Procedure 1
ctrlfcn.idxfcn  = idx2;
ctrlfcn.XArYbnd = q;
[pdK_p12,constOut] = lpvsyn(pdGaw,3,2,const,[],ctrlfcn);
g_p12  = constOut(1).bound;
Gaw = ss(pdGaw,pv.points);
K_p12 = ss(pdK_p12,pv.points);
Gcl = lft(Gaw,K_p12);
pdGcl = ppss(Gcl,pv);
constOutA = lpvanalysis(pdGcl,const);
g_p12r = constOutA(1).bound;

% ------------------------------------------------------------------------
% Design -> removing: p(2); p(1)*p(2);

idx3 = Is{4};

% Procedure 2
[pdK_p23,constOut] = lpvsyn(pdGaw,3,2,const,[],idx3);
g_p23  = constOut(1).bound;
Gaw = ss(pdGaw,pv.points);
K_p23 = ss(pdK_p23,pv.points);
Gcl = lft(Gaw,K_p23);
pdGcl = ppss(Gcl,pv);
constOutA = lpvanalysis(pdGcl,const);
g_p23r = constOutA(1).bound;

% Procedure 1
ctrlfcn.idxfcn  = idx3;
ctrlfcn.XArYbnd = q;
[pdK_p13,constOut] = lpvsyn(pdGaw,3,2,const,[],ctrlfcn);
g_p13  = constOut(1).bound;
Gaw = ss(pdGaw,pv.points);
K_p13 = ss(pdK_p13,pv.points);
Gcl = lft(Gaw,K_p13);
pdGcl = ppss(Gcl,pv);
constOutA = lpvanalysis(pdGcl,const);
g_p13r = constOutA(1).bound;


% results
fprintf('\n')
fprintf('-----------------------------------------------\n')
fprintf('Comparison for Hinf synthesis\n')
fprintf('LPV with all terms:\n')
fprintf(' Kf:                g = %6.4f\n',g_f0)
strJs = sprintf('%1.0f, ',idx0);
fprintf('LPV with I = {%s}:\n',strJs(1:end-2));
fprintf(' Kp2f: Procedure 2: g = %6.4f, g_real = %6.4f\n',g_th0,g_th0r)
strJs = sprintf('%1.0f, ',idx1);
fprintf('LPV with I = {%s}:\n',strJs(1:end-2));
fprintf(' Kp11: Procedure 1: g = %6.4f, g_real = %6.4f\n',g_p11,g_p11r)
fprintf(' Kp21: Procedure 2: g = %6.4f, g_real = %6.4f\n',g_p21,g_p21r)
strJs = sprintf('%1.0f, ',idx2);
fprintf('LPV with I = {%s}:\n',strJs(1:end-2));
fprintf(' Kp12: Procedure 1: g = %6.4f, g_real = %6.4f\n',g_p12,g_p12r)
fprintf(' Kp22: Procedure 2: g = %6.4f, g_real = %6.4f\n',g_p22,g_p22r)
strJs = sprintf('%1.0f, ',idx3);
fprintf('LPV with I = {%s}:\n',strJs(1:end-2));
fprintf(' Kp13: Procedure 1: g = %6.4f, g_real = %6.4f\n',g_p13,g_p13r)
fprintf(' Kp23: Procedure 2: g = %6.4f, g_real = %6.4f\n',g_p23,g_p23r)
fprintf('-----------------------------------------------\n')




%% =============================================================
% Step responses

tsim = 5;
t  = linspace(0,tsim,100);

pdGau = connect(pdG,M,sb,{'r','u'},{'y','u','ei'});
points = pgrid(range,[3 3]);
np = size(points,2);
G  = ss(pdGau,points);

% steps responses

% full
Kf = ss(pdKf,points);
Gclf = lft(G,Kf);
yf = zeros(size(t'));
for ii = 1:np
    yf(:,ii) = step(Gclf(1,1,ii),t);
end

% removing: p(1)*p(2)
Kp21 = ss(pdK_p21,points); Gclp21 = lft(G,Kp21);
Kp11 = ss(pdK_p11,points); Gclp11 = lft(G,Kp11);
yp21 = zeros(size(t')); yp11 = zeros(size(t'));
for ii = 1:np
    yp21(:,ii) = step(Gclp21(1,1,ii),t);
    yp11(:,ii) = step(Gclp11(1,1,ii),t);
end

% removing: p(2)
Kp22 = ss(pdK_p22,points); Gclp22 = lft(G,Kp22);
Kp12 = ss(pdK_p12,points); Gclp12 = lft(G,Kp12);
yp22 = zeros(size(t')); yp12 = zeros(size(t'));
for ii = 1:np
    yp22(:,ii) = step(Gclp22(1,1,ii),t);
    yp12(:,ii) = step(Gclp12(1,1,ii),t);
end

% removing: p(2); p(1)*p(2);
Kp23 = ss(pdK_p23,points); Gclp23 = lft(G,Kp23);
Kp13 = ss(pdK_p13,points); Gclp13 = lft(G,Kp13);
yp23 = zeros(size(t')); yp13 = zeros(size(t'));
for ii = 1:np
    yp23(:,ii) = step(Gclp23(1,1,ii),t);
    yp13(:,ii) = step(Gclp13(1,1,ii),t);
end

cl = lines(2);
gray = 0.5*[1 1 1];

options.ax.sepy      = 1.0;
options.ax.offy      = 1.5;
options.ax.indHeight = 6.0;
options.ax.indWidth  = 0.85;
options.ax.FontSize  = 14;
options.ax.Xlim      = [0 tsim];
Ha = subfigures(3,options);

% full
hl = plot(Ha(1), t,yf,...
    'Linewidth',1,...
    'Color',cl(1,:));
ylabel(Ha(1),'Output')
Ha(1).YLim  = [0 1.2];
Ha(1).YTick = 0:0.2:1.2;
legend(hl(1),{'$K_\mathrm{f}$'},'Interpreter','Latex','Location','SouthEast'); 
legend(Ha(1), 'boxoff');
title(Ha(1),'Controller with all terms','FontSize',12,'FontWeight','normal')

% removing: p(1)*p(2)
hl1 = plot(Ha(2), t, yp22,...
    'Linewidth',1,'LineStyle','-',...
    'Color',cl(1,:));
hl2 = plot(Ha(2), t, yp12,...
    'Linewidth',1,'LineStyle','-',...
    'Color',cl(2,:));
ylabel(Ha(2),'Output')
Ha(2).YLim  = [0 1.2];
Ha(2).YTick = 0:0.2:1.2;
legend([hl1(1) hl2(1)],{'$K_{\mathrm{p}2,1}$','$K_{\mathrm{p}1,1}$'},...
    'Interpreter','Latex','Location','SouthEast')
legend(Ha(2), 'boxoff');
auxStr4 = sprintf(repmat('%1.0f,',1,size(Is{2},2)),Is{2});
strTitle = sprintf('Controller with terms {%s}',auxStr4(1:end-1));
title(Ha(2),strTitle,'FontSize',12,'FontWeight','normal')

% removing: p(2); p(1)*p(2);
hl1 = plot(Ha(3), t, yp23,...
    'Linewidth',1,'LineStyle','-',...
    'Color',cl(1,:));
hl2 = plot(Ha(3), t, yp13,...
    'Linewidth',1,'LineStyle','-',...
    'Color',cl(2,:));
ylabel(Ha(3),'Output')
Ha(3).YLim  = [0 1.2];
Ha(3).YTick = 0:0.2:1.2;
legend([hl1(1) hl2(1)],{'$K_{\mathrm{p}1,3}$','$K_{\mathrm{p}2,3}$'},...
    'Interpreter','Latex','Location','SouthEast')
legend(Ha(3), 'boxoff');
auxStr4 = sprintf(repmat('%1.0f,',1,size(Is{4},2)),Is{4});
strTitle = sprintf('Controller with terms {%s}',auxStr4(1:end-1));
title(Ha(3),strTitle,'FontSize',12,'FontWeight','normal')

