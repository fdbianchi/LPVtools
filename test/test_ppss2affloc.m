
% test ppss2affloc

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
Points = pgrid([2, 4; 0, 20], N);
rate = [-1.5 1.5;-100 100];
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
% checking each simplex

for idx = 1:size(Pset.simplices,1)
 
    simplex = Pset.simplices(idx, :);
    nv = length(simplex);

    pdGloc = ppss2affloc(pdG, idx);

    % a point in the simplex
    alpha = rand(nv,1);
    alpha = alpha/sum(alpha);
    ptest = Pset.points(:,simplex)*alpha;
    
    p0 = Pset.points(:,simplex)*ones(3,1)/3;
    ptestloc = ptest - p0;

    % original model
    G1 = ss(pdG, ptest);

    % new mocel
    G2 = ss(pdGloc, ptestloc);

    figure
    bode(G1, G2);
    legend('Original','Local')
    title(sprintf('Comparison for simplex %d',idx))
    
end
