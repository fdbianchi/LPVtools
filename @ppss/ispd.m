function [dgn,msg] = ispd(pdG)

% ISPD(pdG) checks if pdG is parameter dependent
%   
% Use:
%   dgn = ISPD(pdG,idxo,idxi)
%
% Inputs:
%   - pdG: pass object
%   - idxo, idxi: IO map
%
% Output:
%   - dgn(1) = true: the matrix A is parameter dependent
%   - dgn(2) = true: Input channel is parameter dependent
%   - dgn(3) = true: Output channel is parameter dependent
%
% See also ppss

% fbianchi - 2021-03-30


tol = 1e-6; 
dgn = false(1,3);
if nargin < 2
    [no,ni] = iosize(pdG);
    idxo = 1:no;
    idxi = 1:ni;
end

%
for ii = 2:nsys(pdG)
    % i-th system
    if (norm(pdG.A(:,:,ii) - pdG.A(:,:,1),1) > tol)
        dgn(1) = true;
    end
    if (norm(pdG.B(:,idxi,ii) - pdG.B(:,idxi,1),1) > tol || ...
        norm(pdG.D(:,idxi,ii) - pdG.D(:,idxi,1),1) > tol) 
        dgn(2) = true;
    end
    if (norm(pdG.C(idxo,:,ii) - pdG.C(idxo,:,1),1) > tol || ...
        norm(pdG.D(idxo,:,ii) - pdG.D(idxo,:,1),1) > tol)
        dgn(3) = true;
    end
end

if all(dgn(2:3) == [true false])
    msg = 'Input channel is parameter dependent';
elseif all(dgn(2:3) == [false true])
    msg = 'Output channel is parameter dependent';
elseif all(dgn(2:3) == [true true])
    msg = 'Input and output channels are parameter dependent';
else
    msg = 'Input and output channels are parameter INdependent';
end