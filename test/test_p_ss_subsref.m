
% -------------------------------------------------------------------------
% Test for subrefs
%
% -------------------------------------------------------------------------

% fbianchi - 2024-12-13

% cleaning
clc
clear all
close all

fprintf('\n------------------------------------------------------------\n')
fprintf('pass subsref\n')


% -----------------------------------------------------------
% LPV outputs

G = rss(5,3,4);
G.y = 'y';
G.u = 'u';

G.InputGroup.u1 = [1 2];
G.InputGroup.u2 = [2 3 4];
G.InputGroup.u3 = [3 4];

G.OutputGroup.y1 = [1 2];
G.OutputGroup.y2 = [2 3];

pdG = ppss(G);

% -----------------------------------------------------------
% case 1:
fprintf('\ncase 1: \n')

idx1 = 1:2;
idx2 = 2:4;
G1 = G(idx1,idx2);
pdG1 = pdG(idx1,idx2);

figure
bodemag(pdG1,G1)

% checking input group
sg1 = G1.InputGroup;
sg2 = pdG1.InputGroup;
if isequal(sg1,sg2)
    fprintf('InputGroup -> ok\n')
else
    fprintf('InputGroup -> ERROR\n')
end

% checking output group
sg1 = G1.OutputGroup;
sg2 = pdG1.OutputGroup;
if isequal(sg1,sg2)
    fprintf('OutputGroup -> ok\n')
else
    fprintf('OutputGroup -> ERROR\n')
end

% checking output group
sg1 = G1.y;
sg2 = pdG1.y;
if isequal(sg1,sg2)
    fprintf('Input -> ok\n')
else
    fprintf('Input -> ERROR\n')
end

% checking output group
sg1 = G1.u;
sg2 = pdG1.u;
if isequal(sg1,sg2)
    fprintf('Output -> ok\n')
else
    fprintf('Output -> ERROR\n')
end

% -----------------------------------------------------------
% case 1: names
fprintf('\ncase 1 (names): \n')

idx1 = {'y(1)','y(2)'};
idx2 = {'u(2)','u(3)','u(4)'};
G1 = G(idx1,idx2);
pdG1 = pdG(idx1,idx2);

figure
bodemag(pdG1,G1)

% checking input group
sg1 = G1.InputGroup;
sg2 = pdG1.InputGroup;
if isequal(sg1,sg2)
    fprintf('InputGroup -> ok\n')
else
    fprintf('InputGroup -> ERROR\n')
end

% checking output group
sg1 = G1.OutputGroup;
sg2 = pdG1.OutputGroup;
if isequal(sg1,sg2)
    fprintf('OutputGroup -> ok\n')
else
    fprintf('OutputGroup -> ERROR\n')
end

% checking output group
sg1 = G1.y;
sg2 = pdG1.y;
if isequal(sg1,sg2)
    fprintf('Input -> ok\n')
else
    fprintf('Input -> ERROR\n')
end

% checking output group
sg1 = G1.u;
sg2 = pdG1.u;
if isequal(sg1,sg2)
    fprintf('Output -> ok\n')
else
    fprintf('Output -> ERROR\n')
end

% -----------------------------------------------------------
% case 2:
fprintf('\ncase 2: \n')

idx1 = 3;
idx2 = 4;
G2 = G(idx1,idx2);
pdG2 = pdG(idx1,idx2);

figure
bodemag(pdG2,G2)

% checking input group
sg1 = G2.InputGroup;
sg2 = pdG2.InputGroup;
if isequal(sg1,sg2)
    fprintf('InputGroup -> ok\n')
else
    fprintf('InputGroup -> ERROR\n')
end

% checking output group
sg1 = G2.OutputGroup;
sg2 = pdG2.OutputGroup;
if isequal(sg1,sg2)
    fprintf('OutputGroup -> ok\n')
else
    fprintf('OutputGroup -> ERROR\n')
end

% checking output
sg1 = G2.y;
sg2 = pdG2.y;
if isequal(sg1,sg2)
    fprintf('Output -> ok\n')
else
    fprintf('Output -> ERROR\n')
end

% checking input
sg1 = G2.u;
sg2 = pdG2.u;
if isequal(sg1,sg2)
    fprintf('Input -> ok\n')
else
    fprintf('Input -> ERROR\n')
end


% -----------------------------------------------------------
% case 3:
fprintf('\ncase 3: \n')

idx2 = 4;
G2 = G(:,idx2);
pdG2 = pdG(1:3,idx2);

figure
bodemag(pdG2,G2)


% -----------------------------------------------------------
% case 4:
fprintf('\ncase 4: using inputgruops\n')

G2 = G('y1','u2');
pdG2 = pdG('y1','u2');

figure
bodemag(pdG2,G2)
