
% -------------------------------------------------------------------------
% How to use and Examples for PSET.GRAL
%
% -------------------------------------------------------------------------

% fbianchi - 2021-03-26

% *****************
clc; 
clear all; 
close all
% *****************

rng(23656);

%% ------------------------------------------------------------------------
% 1d example

fprintf('\n')
fprintf('---------------------------------------------------------------\n')
fprintf('Checking PSET.GRAL with 1 parameter\n')

points = rand(1,5);

% only points
setGral1a = pset.Gral(points);
disp('Set defition with points only')
disp(setGral1a)
% points + rate
setGral1b = pset.Gral(points,[-1 1]);
disp('Set defition with points and rates')
disp(setGral1b)
% point + rate + names
setGral1c = pset.Gral(points,[-1 1],{'r'});
disp('Set defition with points, rates and names')
disp(setGral1c)

% check if the set is empty
bool = isempty(setGral1a);
if ~bool
    disp('ISEMPTY working correctly')
else
    error('error in ISEMPTY')
end

% get set size (np = # parameters, nv = # of points)
[np,nv] = size(setGral1a);
if (np == 1) && (nv == 5)
    disp('SIZE working correctly')
else
    error('error in SIZE')
end

% check if a value is in the set
p1 = 0.5;    % in the set
[bool1,msg1] = checkval(setGral1a,p1);
p2 = -1;    % out the set
[bool2,msg2] = checkval(setGral1a,p2);
if (bool1 == 1) && (bool2 == 0)
    disp('CHECKVAL working correctly')
else
    error('error in CHECKVAL')
end

% plot
figure
plot(setGral1a)

% testing cvxdec
p = 0.5;
[alpha,idx] = cvxdec(setGral1a,p);
ps = setGral1a.points(:,idx);
pc = ps*alpha;
plot(setGral1a.points(1,idx([1 2])),[0 0],'b-s')
title('Example of 1D set')
plot(p(1,:),0,'rx')
plot(pc(1,:),0,'ms')
if all(abs(p - pc) <= eps)
    disp('CVXDEC working correctly')
else
    error('error in CVXDEC')
end

%% ------------------------------------------------------------------------
% 2d example

fprintf('\n')
fprintf('---------------------------------------------------------------\n')
fprintf('Checking PSET.GRAL with 2 parameters\n')

points = rand(2,8);

% only points
setGral2a = pset.Gral(points);
disp('Set defition with points only')
disp(setGral2a)
% points + rate
setGral2b = pset.Gral(points,[-1 1;-1 1]);
disp('Set defition with points and rates')
disp(setGral2b)
% point + rate + names
setGral2c = pset.Gral(points,[-1 1],{'r','e'});
disp('Set defition with points, rates and names')
disp(setGral2c)

% check if the set is empty
bool = isempty(setGral2a);
if ~bool
    disp('ISEMPTY working correctly')
else
    error('error in ISEMPTY')
end

% get set size (np = # parameters, nv = # of points)
[np,nv] = size(setGral2a);
if (np == 2) && (nv == 8)
    disp('SIZE working correctly')
else
    error('error in SIZE')
end

% check if a value is in the set
p1 = [0.5;0.7];    % in the set
[bool1,msg1] = checkval(setGral2a,p1);
p2 = [0;-1];    % out the set
[bool2,msg2] = checkval(setGral2a,p2);
if (bool1 == 1) && (bool2 == 0)
    disp('CHECKVAL working correctly')
else
    error('error in CHECKVAL')
end

% plot
figure
plot(setGral2a)
title('Example of 2D set')

% testing cvxdec
p = [0.5;0.4];
[alpha,idx] = cvxdec(setGral2a,p);
ps = setGral2a.points(:,idx);
pc = ps*alpha;
plot(setGral2a.points(1,idx([1 2 3 1])),setGral2a.points(2,idx([1 2 3 1])),'b-s')
plot(p(1,:),p(2,:),'rx')
plot(pc(1,:),pc(2,:),'ms')
if all(abs(p - pc) <= eps)
    disp('CVXDEC working correctly')
else
    error('error in CVXDEC')
end


%% ------------------------------------------------------------------------
% 3d example

fprintf('\n')
fprintf('---------------------------------------------------------------\n')
fprintf('Checking PSET.GRAL with 3 parameters\n')

points = rand(3,10);

setGral3 = pset.Gral(points);
disp('Set defition with points only')
disp(setGral3)

% check if the set is empty
bool = isempty(setGral3);
if ~bool
    disp('ISEMPTY working correctly')
else
    error('error in ISEMPTY')
end

% get set size (np = # parameters, nv = # of points)
[np,nv] = size(setGral3);
if (np == 3) && (nv == 10)
    disp('SIZE working correctly')
else
    error('error in SIZE')
end

% check if a value is in the set
p1 = [0.5;0.4;0.6];    % in the set
[bool1,msg1] = checkval(setGral3,p1);
p2 = [0;-1;0];    % out the set
[bool2,msg2] = checkval(setGral3,p2);
if (bool1 == 1) && (bool2 == 0)
    disp('CHECKVAL working correctly')
else
    error('error in CHECKVAL')
end

% testing cvxdec
p = [0.5;0.4;0.6];
[alpha,idx] = cvxdec(setGral3,p);
ps = setGral3.points(:,idx);
pc = ps*alpha;
if all(abs(p - pc) <= eps)
    disp('CVXDEC working correctly')
else
    error('error in CVXDEC')
end

% plot
figure
plot(setGral3)
plot3(ps(1,:),ps(2,:),ps(3,:),'bs')
plot3(p(1,:),p(2,:),p(3,:),'rx')
plot3(pc(1,:),pc(2,:),pc(3,:),'ms')
title('Example of 3D set')


%% ------------------------------------------------------------------------
% conversion from lmitool object

fprintf('\n')
fprintf('---------------------------------------------------------------\n')
fprintf('Conversion from lmitool object\n')

% from box lmitool obj to pset
setLMI1 = pvec('pol',rand(2,4));
setGral5 = pset.Gral('pvec',setLMI1);
disp('from PVEC to PSET.BOX')
disp(setGral5)

setLMI2 = pvec(setGral5);
disp('from PSET.BOX to PVEC')
pvinfo(setLMI2)

% ------------------------------------------------------------------------
% subset

fprintf('\n')
fprintf('---------------------------------------------------------------\n')
fprintf('Subsets\n')

disp('Original set')
P = [1 1 4 3;2 2 3 4;2 3 7 8];
setGral9 = pset.Gral(P);
disp(setGral9)

disp('Subset with parameters 1 and 3')
setGral9a = setGral9([1 2]);
disp(setGral9a)

figure
plot(setGral9)
hold on
xa = setGral9a.points(1,:);
ya = setGral9a.points(2,:);
patch(xa,ya,P(3,1)*[1 1 1],0.5*[1 1 1],'EdgeColor',0.5*[1 1 1],'FaceAlpha',0.3)


%% ------------------------------------------------------------------------
% extreme cases

fprintf('\n')
fprintf('---------------------------------------------------------------\n')
fprintf('Extreme cases\n')

fprintf('\nEmpty sets\n')
setGral4 = pset.Gral();
bool = isempty(setGral4);
if (bool == 1)
    disp('ISEMPTY working correctly')
else
    error('error in ISEMPTY')
end
[np,nv] = size(setGral4);
if (np == 0) && (nv == 0)
    disp('SIZE working correctly')
else
    error('error in SIZE')
end
[bool1,msg1] = checkval(setGral4,[]);
[bool2,msg2] = checkval(setGral4,1);
if (bool1 == 1) && (bool2 < 0)
    disp('CHECKVAL working correctly')
else
    error('error in CHECKVAL')
end
try 
    [alpha,idx] = cvxdec(setGral4,1);
catch ME
    disp(ME.message)
    disp('CVXDEC working correctly')
end

fprintf('\nSet 1D with only one point\n')
setGral5 = pset.Gral(5);
[bool1,msg1] = checkval(setGral5,5);
[bool2,msg2] = checkval(setGral5,1.5);
if (bool1 == 1) && (bool2 == 0)
    disp('CHECKVAL working correctly')
else
    error('error in CHECKVAL')
end
p1 = 5;
[alpha,idx] = cvxdec(setGral5,p1);
p1r = setGral5.points(:,idx)*alpha;
if all(abs(p1 - p1r) <= 10*eps)
    disp('CVXDEC working correctly')
else
    error('error in CVXDEC')
end

fprintf('\nSet 3D with only one point\n')
setGral6 = pset.Gral([1;3;5]);
p1 = [1;3;5];
[bool1,msg1] = checkval(setGral6,p1);
p2 = [1;2;5];
[bool2,msg2] = checkval(setGral6,p2);
if (bool1 == 1) && (bool2 == 0)
    disp('CHECKVAL working correctly')
else
    error('error in CHECKVAL')
end
[alpha,idx] = cvxdec(setGral6,p1);
p1r = setGral6.points(:,idx)*alpha;
if all(abs(p1 - p1r) <= 10*eps)
    disp('CVXDEC working correctly')
else
    error('error in CVXDEC')
end


fprintf('\nChecking plots\n')

setGral8 = pset.Gral(1);
figure; plot(setGral8); title('Set 1D with 1 vertex')
setGral8 = pset.Gral([1 1]);
figure; plot(setGral8); title('Set 1D with coincident vertices')

setGral8 = pset.Gral([2;1]);
figure; plot(setGral8); title('Set 2D with 1 vertex')
setGral8 = pset.Gral([1 2;2 3]);
figure; plot(setGral8); title('Set 2D with 2 vertices')
setGral8 = pset.Gral([1 2;2 2]);
figure; plot(setGral8); title('Set 2D with 2 vertices, with one dimension not changing')

setGral8 = pset.Gral([1;2;5]);
figure; plot(setGral8); title('Set 3D with 1 vertex')
setGral8 = pset.Gral([0 1;2 2;3 5]);
figure; plot(setGral8); title('Set 3D with 2 vertices')
setGral8 = pset.Gral([0 1 5;2 2 7;5 6 9]);
figure; plot(setGral8); title('Set 3D with 3 vertices')
setGral8 = pset.Gral([0 1 5;2 2 2;5 6 9]);
figure; plot(setGral8); title('Set 3D with 3 vertices, with one dimension not changing')
