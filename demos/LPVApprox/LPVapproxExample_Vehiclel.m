
% Aproximation of LPV models
%
% Example: Torque Vectoring of 4-wheel-drive Formula Student Vehicle
%
% see: F. Bianchi, R. Sánchez-Peña, A method for reducing implementation 
% complexity in linear parameter-varying controllers, Automatica, Volume 146,
% 2022, 110588, https://doi.org/10.1016/j.automatica.2022.110588.
%
% fbianchi - 2021-02-11


% cleaning
clearvars; 
clc; 
close all

% ========================================================================
% Model parameter

% Vehicle Parameters
lf = 0.8;       % [m]
lr = 0.8;       % [m]
l = lr + lf;
m = 300;        % [kg]
Iz = 100;       % [kg*m^2]

% ========================================================================
% Full LPV model 

% parameter grids
nb = 5; nv = 5; 
% nb = 10; nv = 10; % to check in a denser grid
% parameter set         -> p = [beta,v]
betaMn = 0;         betaMx = 0.2;
vMn = 5;            vMx = 28;
range = [betaMn     betaMx;
         vMn        vMx];
pv = pset.Grid(range,[nb nv],[],{'beta','v'});

% parameter functions   -> p = [beta,v]
C_r = 36000;
C_f =@(beta) 14000./(1100*beta.^2 + 1);
% 
fcnpar =@(p) [C_f(p(1));
              C_f(p(1))/p(2);
              1/p(2);
              C_f(p(1))/p(2)^2;
              1/p(2)^2];
% minimum & maximum values of fcnpar
f1 = @(p) C_f(p);               % fplot(f1,[betaMn,betaMx]);
fparMn(1,1) = C_f(betaMx);
fparMx(1,1) = C_f(0);  
f2 = @(p1,p2) C_f(p1)/p2;       % fsurf(f2,range');
fparMn(2,1) = f2(betaMx,vMx);
fparMx(2,1) = f2(0,vMn);
f3 = @(p) 1/p;                  % fplot(f3,[vMn,vMx]);
fparMn(3,1) = f3(vMx);
fparMx(3,1) = f3(vMn);
f4 = @(p1,p2) C_f(p1)/p2^2;     % fsurf(f4,range');
fparMn(4,1) = f4(betaMx,vMx);
fparMx(4,1) = f4(0,vMn);
f5 = @(p) 1/p^2;                % fplot(f5,[vMn,vMx]);
fparMn(5,1) = f5(vMx);
fparMx(5,1) = f5(vMn);
          
% system matrix          
Ap(:,:,1) = [0       -1;          C_r*lr/Iz    0];             % independent term
Ap(:,:,2) = [0        0;         -lf/Iz        0];             % f_1 term
Ap(:,:,3) = [-1/m     0;          0           -lf^2/Iz];       % f_2 term
Ap(:,:,4) = [-C_r/m   0;          0           -C_r*lr^2/Iz];   % f_3 term
Ap(:,:,5) = [0       -lf/m;       0            0];             % f_4 term
Ap(:,:,6) = [0        C_r*lr/m;   0            0];             % f_5 term

Sc = diag(1./[1 1e-4]);
Bp(:,:,1) = [0        0;    0      1/Iz]*Sc;                   % independent term
Bp(:,:,2) = [0        0;    lf/Iz  0]*Sc;                      % f_1 term
Bp(:,:,3) = [1/m      0;    0      0]*Sc;                      % f_2 term
Bp(:,:,4) = [0        0;    0      0]*Sc;                      % f_3 term
Bp(:,:,5) = [0        0;    0      0]*Sc;                      % f_4 term
Bp(:,:,6) = [0        0;    0      0]*Sc;                      % f_5 term
Bp(:,1,:) = [];

Cp = [0 1];
Dp = 0;

% LPV model
pdG = pgss(Ap,Bp,Cp,Dp,pv,fcnpar);
pdG.u = {'u'};    pdG.y = 'y';

% ------------------------------------------------------------------------
% Full Normalized LPV model 

np   = size(pv.points,2);
fpar = zeros(5,np);
for ii = 1:np
    fpar(:,ii) = fcnpar(pv.points(:,ii));
end
fcnparN =@(p) (fcnpar(p) - fparMn)./(fparMx - fparMn);

An = zeros(2,2,6);  An(:,:,1) = Ap(:,:,1);
deltaf = fparMx - fparMn; 
a1 = 1./deltaf;  a2 = -fparMn./deltaf;
for ii = 2:6
    An(:,:,1)  = An(:,:,1) - Ap(:,:,ii)*a2(ii-1)/a1(ii-1);
    An(:,:,ii) = Ap(:,:,ii)/a1(ii-1);
end
Bn = Bp;
Cn = Cp;
Dn = Dp;
pdGn = pgss(An,Bn,Cn,Dn,pv,fcnparN);
pdGn.u = 'u';    pdGn.y = 'y';

% figure; bodemag(pdG,pdGn)

% ------------------------------------------------------------------------
% Augmented plant

% tracking weight
fe = 10;        ke = 2;
We = ke*tf(10*[1/(10*fe) 1],[1/(0.1*fe) 1]);
% control weight
fu = 25;        ku = 0.1;
Wu = tf(ku*[1/(0.1*fu) 1],[1/(10*fu) 1]);
% figure; bodemag(We,Wu)

sb    = sumblk('e = r - y');
% weighted augmented plant
pdGau = connect(pdG,sb,{'r','u'},{'e','u','e'});
pdGaw = balreal(append(We,Wu,1))*pdGau;

% weighted augmented plant (normalized)
pdGau = connect(pdGn,sb,{'r','u'},{'e','u','e'});
pdGawn = balreal(append(We,Wu,1))*pdGau;

% figure; bodemag(pdGaw,pdGawn)


% ========================================================================
% SVDs analysis

fprintf('\n')
fprintf('---------------------------------------------------------------\n')
fprintf('Max SVD for each term in the matrix A of the augmented\n')
fprintf('  (A original plant, A_n normalized plant)\n')
fprintf('\n')
A = pdGaw.A; AN = pdGawn.A;
for ii = 1:nsys(pdGaw)
    sg  = max(svd(A(:,:,ii)));
    sgN = max(svd(AN(:,:,ii)));
    if ii == 1
        deltafi = 1;
    else
        deltafi = deltaf(ii-1);
    end        
    fprintf('svd(A_%1.0f): %8.4f (delta f_i: %8.2f),\tsvd(A_n%1.0f): %8.4f\n',...
        ii,sg,deltafi,ii,sgN)
end
fprintf('\n')

% set of included terms
Is = {[1 2 4 5 6];
      [1 2 3 4 6];
      [1 2 3 4 5];
      [1 2 4 6];
      [1 2 4 5];
      [1 2 3 4];
      [1 2 4];
      };
% set of functions
strFcnpar = {'C_f(p(1))', 'C_f(p(1))/p(2)', '1/p(2)', ...
             'C_f(p(1))/p(2)^2', '1/p(2)^2'};
Us = 1:nsys(pdGaw);

fprintf('\n-----------------------------------------------------------\n')
fprintf('Testing removing several terms in the matrix A\n')
fprintf('\tI indices of included terms, J indices of removed terms\n')  
fprintf('\tAs: matrix with terms in I, Ar: matrix with terms in J\n')  
for kk = 1:size(Is,1)
    
    fprintf('\n')
    I = Is{kk}; J = setdiff(Us,I);
    auxStr1 = repmat('%1.0f,',1,size(I,2));
    auxStr2 = repmat('%1.0f,',1,size(J,2));
    auxStr3 = ['Test for I = {' auxStr1(1:end-1) '}, J = {' auxStr2(1:end-1) '}\n'];
    fprintf(auxStr3,I,J);
    
    auxStr4 = sprintf('%s, ',strFcnpar{I(2:end)-1});
    fprintf('\t fcn = {%s}\n',auxStr4(1:end-2));
    
    sgS  = zeros(1,np); sgR  = sgS;
    sgSn = zeros(1,np); sgRn = sgS;
    for ii = 1:np
        p     = pv.points(:,ii);
        fparN = pdGawn.parfcn(p);
        fpar  = pdGaw.parfcn(p);
        
        As = 0; Asn = 0;
        for jj = I
            if jj == 1
                As  = As + A(:,:,jj);
                Asn = Asn + AN(:,:,jj);
            else
                As  = As + A(:,:,jj)*fpar(jj-1);
                Asn = Asn + AN(:,:,jj)*fparN(jj-1);
            end
        end
        Ar = 0; Arn = 0;
        for jj = J
            if jj == 1
                Ar  = Ar + A(:,:,jj);
                Arn = Arn + AN(:,:,jj);
            else
                Ar  = Ar + A(:,:,jj)*fpar(jj-1);
                Arn = Arn + AN(:,:,jj)*fparN(jj-1);
            end
        end
        sgS(ii)  = max(svd(As));      sgR(ii)  = max(svd(Ar));
        sgSn(ii) = max(svd(Asn));     sgRn(ii) = max(svd(Arn));
        
    end
    fprintf('\t max svd(As): %6.4f, max svd(Ar): %6.4f\n',max(sgS),max(sgR))
    fprintf('\t max svd(As): %6.4f, max svd(Ar): %6.4f (Normalized)\n',max(sgSn),max(sgRn))

end
fprintf('\n')

% ========================================================================
% LPV design

opts = lpvsettings('solver','mosek');

const(1) = synConst.Gain(1,1:2);
const(2) = synConst.Poles('MaxFreq',50000,'MinDamping',0.1);

% -----------------------------------------------------------------------
% Design -> All terms
[pdKf,constOut] = lpvsyn(pdGaw,3,2,const,[],[],opts);
gamf = constOut(1).bound;
pdGcl = lft(pdGaw,pdKf);
constOutA = lpvanalysis(pdGcl,const,[],[],opts);
gamf_r = constOutA(1).bound;

% -----------------------------------------------------------------------
% Design -> with only the terms C_f(p(1)), C_f(p(1))/p(2), 1/p(2)

idx1 = Is{6};

% Procedure 2 (conservative)
[pdK_p21,constOut] = lpvsyn(pdGaw,3,2,const,[],idx1,opts);
g_p21 = constOut(1).bound;
pdGcl = ppss(lft(pdGaw,pdK_p21),pv);
constOutA = lpvanalysis(pdGcl,const,[],[],opts);
g_p21_r = constOutA(1).bound;

% Procedure 1 (optimistic)
ctrlfcn.idxfcn = idx1;
ctrlfcn.XArYbnd = 1e-3;
[pdK_p11,constOut] = lpvsyn(pdGaw,3,2,const,[],ctrlfcn,opts);
g_p11 = constOut(1).bound;
pdGcl = ppss(lft(pdGaw,pdK_p11),pv);
constOutA = lpvanalysis(pdGcl,const,[],[],opts);
g_p11_r = constOutA(1).bound;

% -----------------------------------------------------------------------
% Design -> with only the terms C_f(p(1)), 1/p(2)

idx2 = Is{7};

% Procedure 2 (conservative)
[pdK_p22,constOut] = lpvsyn(pdGaw,3,2,const,[],idx2,opts);
g_p22 = constOut(1).bound;
pdGcl = ppss(lft(pdGaw,pdK_p22),pv);
constOutA = lpvanalysis(pdGcl,const,[],[],opts);
g_p22_r = constOutA(1).bound;

% Procedure 1 (optimistic)
ctrlfcn.idxfcn = idx2;
ctrlfcn.XArYbnd = 1e-3;
const(2) = synConst.Poles('MaxFreq',20000000,'MinDamping',0.1); % to ensure stable CL plant
[pdK_p12,constOut] = lpvsyn(pdGaw,3,2,const,[],ctrlfcn,opts);
g_p12 = constOut(1).bound;
pdGcl = ppss(lft(pdGaw,pdK_p12),pv);
constOutA = lpvanalysis(pdGcl,const,[],[],opts);
g_p12_r = constOutA(1).bound;


% results:
fprintf('\n')
fprintf('-------------------------------------------------\n')
fprintf('Comparison for Hinf synthesis\n\n')
fprintf('LPV with all terms:\n')
fprintf('  Kf:                g = %6.4f, g_real = %6.4f\n',gamf,gamf_r)
strIs = sprintf('%s ',strFcnpar{idx1(2:end)-1});
fprintf('LPV [%s]:\n',strIs(1:end-1));
fprintf('  Kp11: Procedure 1: g = %6.4f, g_real = %6.4f\n',g_p11,g_p11_r)
fprintf('  Kp21: Procedure 2: g = %6.4f, g_real = %6.4f\n',g_p21,g_p21_r)
strIs = sprintf('%s ',strFcnpar{idx2(2:end)-1});
fprintf('LPV [%s]:\n',strIs(1:end-1));
fprintf('  Kp12: Procedure 1: g = %6.4f, g_real = %6.4f\n',g_p12,g_p12_r)
fprintf('  Kp22: Procedure 2: g = %6.4f, g_real = %6.4f\n',g_p22,g_p22_r)
fprintf('-------------------------------------------------\n')


%% ===============================================================
% Step responses

tsim = 1; m = 200;
t = linspace(0,tsim,m);

points = pgrid(pdG.parset.range,[3 3]);
np = size(points,2);
w = logspace(-1,3,150);

pdGau = connect(pdG,sb,{'r','u'},{'e','u','e'});
pdGaw = append(We,Wu,1)*pdGau;
Gaw = ss(pdGaw,points); 
pdGau = connect(pdG,sb,{'r','u'},{'y','u','e'});
G = ss(pdGau,points);

% full
K = ss(pdKf,points);
Gclf = lft(G,K);
% figure, bodemag(lft(Gaw,K),w)
yf = zeros(m,np);
for ii = 1:np
    yf(:,ii) = step(Gclf(1,1,ii),t);
end
% removing: C_f(p(1))/p(2)^2; 1/p(2)^2 
K_p21 = ss(pdK_p21,points); Gcl_p21 = lft(G,K_p21);
K_p11 = ss(pdK_p11,points); Gcl_p11 = lft(G,K_p11);
% figure, bodemag(lft(Gaw,Ka),lft(Gaw,Kb),w)
y_p21 = zeros(m,np); y_p11 = zeros(m,np);
for ii = 1:np
    y_p21(:,ii) = step(Gcl_p21(1,1,ii),t);
    y_p11(:,ii) = step(Gcl_p11(1,1,ii),t);
end
% removing: C_f(p(1))/p(2), C_f(p(1))/p(2)^2; 1/p(2)^2 
K_p22 = ss(pdK_p22,points); Gcl_p22 = lft(G,K_p22);
K_p12 = ss(pdK_p12,points); Gcl_p12 = lft(G,K_p12);
% figure, bodemag(lft(Gaw,Ka),lft(Gaw,Kb),w)
y_p22 = zeros(m,np); y_p12 = zeros(m,np);
for ii = 1:np
    y_p22(:,ii) = step(Gcl_p22(1,1,ii),t);
    y_p12(:,ii) = step(Gcl_p12(1,1,ii),t);
end


% -----------------------------------------------------------------------
%% figures

cl = lines(2);
gray = 0.7*[1 1 1];

options.ax.sepy      = 1.0;
options.ax.offy      = 1.5;
options.ax.indHeight = 4.8;
options.ax.indWidth  = 0.85;
options.ax.FontSize  = 14;
Ha = subfigures(3,options);
set(Ha,'YLim',[0 1.5])

k = 1;
plot(Ha(k),[0 tsim],[1 1],'Linewidth',1,'LineStyle','-','Color',gray);
hl = plot(Ha(k), t,yf,'Linewidth',1,'LineStyle','-','Color',cl(1,:));
ylabel(Ha(k),'$y$','Interpreter','Latex')
legend(hl(1),{'$K_\mathrm{f}$'},'Interpreter','Latex',...
    'Location','NorthEast','Orientation','Horizontal'); 
legend(Ha(k), 'boxoff');
S = stepinfo(yf,t);
rt = [S.RiseTime];      rtmx = max(rt); rtmn = min(rt);
st = [S.SettlingTime];  stmx = max(st); stmn = min(st);
os = [S.Overshoot];     osmx = max(os); osmn = min(os);
es = 1 - yf(end,:);     esmx = max(es); esmn = min(es);
xt = 0.85; yt = 0.15; 
text(Ha(k),xt,yt+3*0.15,sprintf('$%4.2f \\leq t_r \\leq %4.2f$',rtmn,rtmx),...
    'FontName','Times','FontSize',14,'Interpreter','Latex',...
    'HorizontalAlignment','center','Color',cl(1,:))
text(Ha(k),xt,yt+2*0.15,sprintf('$%4.2f \\leq t_s \\leq %4.2f$',stmn,stmx),...
    'FontName','Times','FontSize',14,'Interpreter','Latex',...
    'HorizontalAlignment','center','Color',cl(1,:))
text(Ha(k),xt,yt+1*0.15,sprintf('$%4.2f \\leq M_p \\leq %4.2f$',osmn,osmx),...
    'FontName','Times','FontSize',14,'Interpreter','Latex',...
    'HorizontalAlignment','center','Color',cl(1,:))
text(Ha(k),xt,yt+0*0.15,sprintf('$%4.2f \\leq e_{ss} \\leq %4.2f$',esmn,esmx),...
    'FontName','Times','FontSize',14,'Interpreter','Latex',...
    'HorizontalAlignment','center','Color',cl(1,:))
title(Ha(k),'Controller with all terms','FontSize',12,'FontWeight','normal')

k = 2;
plot(Ha(k),[0 tsim],[1 1],'Linewidth',1,'LineStyle','-','Color',gray);
hl1 = plot(Ha(k), t, y_p21,'Linewidth',1,'LineStyle','-','Color',cl(1,:));
hl2 = plot(Ha(k), t, y_p11,'Linewidth',1,'LineStyle','-','Color',cl(2,:));
ylabel(Ha(k),'$y$','Interpreter','Latex')
legend([hl2(1) hl1(1)],{'$K_{\mathrm{p}1,1}$','$K_{\mathrm{p}2,1}$'},...
    'Interpreter','Latex','Location','NorthEast','Orientation','Horizontal')
legend(Ha(k), 'boxoff');
S = stepinfo(y_p21,t);
rt = [S.RiseTime];      rtmx = max(rt); rtmn = min(rt);
st = [S.SettlingTime];  stmx = max(st); stmn = min(st);
os = [S.Overshoot];     osmx = max(os); osmn = min(os);
es = 1 - y_p21(end,:);    esmx = max(es); esmn = min(es);
xt = 0.85; yt = 0.15; 
text(Ha(k),xt,yt+3*0.15,sprintf('$%4.2f \\leq t_r \\leq %4.2f$',rtmn,rtmx),...
    'FontName','Times','FontSize',14,'Interpreter','Latex',...
    'HorizontalAlignment','center','Color',cl(1,:))
text(Ha(k),xt,yt+2*0.15,sprintf('$%4.2f \\leq t_s \\leq %4.2f$',stmn,stmx),...
    'FontName','Times','FontSize',14,'Interpreter','Latex',...
    'HorizontalAlignment','center','Color',cl(1,:))
text(Ha(k),xt,yt+1*0.15,sprintf('$%4.2f \\leq M_p \\leq %4.2f$',osmn,osmx),...
    'FontName','Times','FontSize',14,'Interpreter','Latex',...
    'HorizontalAlignment','center','Color',cl(1,:))
text(Ha(k),xt,yt+0*0.15,sprintf('$%4.2f \\leq e_{ss} \\leq %4.2f$',esmn,esmx),...
    'FontName','Times','FontSize',14,'Interpreter','Latex',...
    'HorizontalAlignment','center','Color',cl(1,:))
S = stepinfo(y_p11,t);
rt = [S.RiseTime];      rtmx = max(rt); rtmn = min(rt);
st = [S.SettlingTime];  stmx = max(st); stmn = min(st);
os = [S.Overshoot];     osmx = max(os); osmn = min(os);
es = 1 - y_p11(end,:);     esmx = max(es); esmn = min(es);
xt = 0.55;
text(Ha(k),xt,yt+3*0.15,sprintf('$%4.2f \\leq t_r \\leq %4.2f$',rtmn,rtmx),...
    'FontName','Times','FontSize',14,'Interpreter','Latex',...
    'HorizontalAlignment','center','Color',cl(2,:))
text(Ha(k),xt,yt+2*0.15,sprintf('$%4.2f \\leq t_s \\leq %4.2f$',stmn,stmx),...
    'FontName','Times','FontSize',14,'Interpreter','Latex',...
    'HorizontalAlignment','center','Color',cl(2,:))
text(Ha(k),xt,yt+1*0.15,sprintf('$%4.2f \\leq M_p \\leq %4.2f$',osmn,osmx),...
    'FontName','Times','FontSize',14,'Interpreter','Latex',...
    'HorizontalAlignment','center','Color',cl(2,:))
text(Ha(k),xt,yt+0*0.15,sprintf('$%4.2f \\leq e_{ss} \\leq %4.2f$',esmn,esmx),...
    'FontName','Times','FontSize',14,'Interpreter','Latex',...
    'HorizontalAlignment','center','Color',cl(2,:))
auxStr4 = sprintf('%s, ',strFcnpar{Is{6}(2:end)-1});
strTitle = sprintf('Controller without terms {%s}',auxStr4(1:end-2));
title(Ha(k),strTitle,'FontSize',12,'FontWeight','normal')

k = 3;
hr  = plot(Ha(k),[0 tsim],[1 1],'Linewidth',1,'LineStyle','-','Color',gray);
hl1 = plot(Ha(k), t, y_p22,'Linewidth',1,'LineStyle','-','Color',cl(1,:));
hl2 = plot(Ha(k), t, y_p12,'Linewidth',1,'LineStyle','-','Color',cl(2,:));
ylabel(Ha(k),'$y$','Interpreter','Latex')
legend([hl2(1) hl1(1)],{'$K_{\mathrm{p}1,2}$','$K_{\mathrm{p}2,2}$'},...
    'Interpreter','Latex','Location','NorthEast','Orientation','Horizontal')
legend(Ha(k), 'boxoff');
S = stepinfo(y_p22,t);
rt = [S.RiseTime];      rtmx = max(rt); rtmn = min(rt);
st = [S.SettlingTime];  stmx = max(st); stmn = min(st);
os = [S.Overshoot];     osmx = max(os); osmn = min(os);
es = 1 - y_p22(end,:);     esmx = max(es); esmn = min(es);
xt = 0.85; yt = 0.15;
text(Ha(k),xt,yt+3*0.15,sprintf('$%4.2f \\leq t_r \\leq %4.2f$',rtmn,rtmx),...
    'FontName','Times','FontSize',14,'Interpreter','Latex',...
    'HorizontalAlignment','center','Color',cl(1,:))
text(Ha(k),xt,yt+2*0.15,sprintf('$%4.2f \\leq t_s \\leq %4.2f$',stmn,stmx),...
    'FontName','Times','FontSize',14,'Interpreter','Latex',...
    'HorizontalAlignment','center','Color',cl(1,:))
text(Ha(k),xt,yt+1*0.15,sprintf('$%4.2f \\leq M_p \\leq %4.2f$',osmn,osmx),...
    'FontName','Times','FontSize',14,'Interpreter','Latex',...
    'HorizontalAlignment','center','Color',cl(1,:))
text(Ha(k),xt,yt+0*0.15,sprintf('$%4.2f \\leq e_{ss} \\leq %4.2f$',esmn,esmx),...
    'FontName','Times','FontSize',14,'Interpreter','Latex',...
    'HorizontalAlignment','center','Color',cl(1,:))
S = stepinfo(y_p12,t);
rt = [S.RiseTime];      rtmx = max(rt); rtmn = min(rt);
st = [S.SettlingTime];  stmx = max(st); stmn = min(st);
os = [S.Overshoot];     osmx = max(os); osmn = min(os);
es = 1 - y_p12(end,:);    esmx = max(es); esmn = min(es);
xt = 0.55;
text(Ha(k),xt,yt+3*0.15,sprintf('$%4.2f \\leq t_r \\leq %4.2f$',rtmn,rtmx),...
    'FontName','Times','FontSize',14,'Interpreter','Latex',...
    'HorizontalAlignment','center','Color',cl(2,:))
text(Ha(k),xt,yt+2*0.15,sprintf('$%4.2f \\leq t_s \\leq %4.2f$',stmn,stmx),...
    'FontName','Times','FontSize',14,'Interpreter','Latex',...
    'HorizontalAlignment','center','Color',cl(2,:))
text(Ha(k),xt,yt+1*0.15,sprintf('$%4.2f \\leq M_p \\leq %4.2f$',osmn,osmx),...
    'FontName','Times','FontSize',14,'Interpreter','Latex',...
    'HorizontalAlignment','center','Color',cl(2,:))
text(Ha(k),xt,yt+0*0.15,sprintf('$%4.2f \\leq e_{ss} \\leq %4.2f$',esmn,esmx),...
    'FontName','Times','FontSize',14,'Interpreter','Latex',...
    'HorizontalAlignment','center','Color',cl(2,:))
auxStr4 = sprintf('%s, ',strFcnpar{Is{7}(2:end)-1});
strTitle = sprintf('Controller without terms {%s}',auxStr4(1:end-2));
title(Ha(k),strTitle,'FontSize',12,'FontWeight','normal')

