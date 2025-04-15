
% Effects of different location of gamma in the BRL

% fbianchi - 2024-11-20


% cleaning
clearvars
close all
clc

p = 1e-2;
N = 10;     % tests
ns = 4;     
nw = 1;
nz = 1;

for ii = 1:N
    
% randon plant
G = rss(ns,nz,nw);
[A,B,C,D] = ssdata(G);
nG(ii) = norm(G,inf);

% ------------------------------------------------------------------------
% option 1

X = sdpvar(ns);
g = sdpvar(1);

% opt = sdpsettings('solver','sedumi');
opt = sdpsettings('solver','mosek','verbose',0);

brl = blkvar;
brl(1,1) = (X*A) + (X*A)';
brl(1,2) = X*B;
brl(2,2) = -eye(nw);
brl(1,3) = C';
brl(2,3) = D';
brl(3,3) = -g*eye(nz);
brl = sdpvar(brl);

d = optimize(brl <= 0, g + p*norm(X), opt);

g1(ii) = sqrt(value(g));
X1(:,:,ii) = value(X);
nX1(ii) = norm(X1(:,:,ii));


% ------------------------------------------------------------------------
% option 2

X = sdpvar(ns);
g = sdpvar(1);

brl = blkvar;
brl(1,1) = (X*A) + (X*A)';
brl(1,2) = X*B;
brl(2,2) = -g*eye(nw);
brl(1,3) = C';
brl(2,3) = D';
brl(3,3) = -g*eye(nz);
brl = sdpvar(brl);

d = optimize(brl <= 0, g + p*norm(X), opt);

g2(ii) = value(g);
X2(:,:,ii) = value(X);
nX2(ii) = norm(X2(:,:,ii));

% ------------------------------------------------------------------------
% option 3

X = sdpvar(ns);
g = sdpvar(1);

brl = blkvar;
brl(1,1) = (X*A) + (X*A)';
brl(1,2) = X*B;
brl(2,2) = -g*eye(nw);
brl(1,3) = C';
brl(2,3) = D';
brl(3,3) = -eye(nz);
brl = sdpvar(brl);

d = optimize(brl <= 0, g + p*norm(X), opt);

g3(ii) = sqrt(value(g));
X3(:,:,ii) = value(X);
nX3(ii) = norm(X3(:,:,ii));

% ------------------------------------------------------------------------
% checking results

fprintf('\n')
fprintf('-----------------------------------------------------\n')
fprintf(' ||G||_inf = %6.3f\n',nG(ii))
fprintf('-----------------------------------------------------\n')
fprintf('Comparison of gamma locations in BRL\n')
fprintf(' (3,3)         : gamma = %6.3f, ||X|| = %6.3f\n',g1(ii),nX1(ii))
fprintf(' (2,2) & (3,3) : gamma = %6.3f, ||X|| = %6.3f\n',g2(ii),nX2(ii))
fprintf(' (2,2)         : gamma = %6.3f, ||X|| = %6.3f\n',g3(ii),nX3(ii))
fprintf('-----------------------------------------------------\n')
fprintf('\n')

end


figure

c = 1:N;

subplot(2,1,1)
plot(c,[g1;g2;g3;nG],'.-')
ylabel('gamma')
legend('(3,3)','(2,2) & (3,3)','(2,2)','||G||')

subplot(2,1,2)
ylabel('||X||')
plot(c,[nX1;nX2;nX3],'.-')

