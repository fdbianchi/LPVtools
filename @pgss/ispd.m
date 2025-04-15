function [dgn,msg] = ispd(pdG)

% ISPD(pdG) checks if pdG is parameter dependent
%   
% Use:
%   dgn = ISPD(pdG)
%   dgn = ISPD(pdG,idxo,idxi)
%
% Inputs:
%   - pdG: pgss object
%   - idxo, idxi: IO map
%
% Output:
%   - dgn(1) = true: the matrix A is parameter dependent
%   - dgn(2) = true: Input channel is parameter dependent
%   - dgn(3) = true: Output channel is parameter dependent
%
% See also pgss

% fbianchi - 2021-03-31


% tol = 1e-6; 
if nargin < 2
    [no,ni] = iosize(pdG);
    idxo = 1:no;
    idxi = 1:ni;
end
%
dgn = false(1,3);
for ii = 2:nsys(pdG)
    % i-th system
    if (norm(pdG.A(:,:,ii),1) > 0)
        dgn(1) = true;
    end
    if (norm(pdG.B(:,idxi,ii),1) > 0 || norm(pdG.D(:,idxi,ii),1) > 0) 
        dgn(2) = true;
    end
    if (norm(pdG.C(idxo,:,ii),1) > 0 || norm(pdG.D(idxo,:,ii),1) > 0)
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