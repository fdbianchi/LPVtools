function pdGint = connect(varargin)

% CONNECT connects the block diagram elements based on i/o names
%
% Use:
%   sysc = CONNECT(sys1,...,sysN,inputs,outputs)
%
% where ONE system can be a pgss object.
%
% See also pgss, connect

% fbianchi - 2023-07-07

% ------------------------------------------------------------------------
% check if there is more than 1 lpv model to interconnect
types = cellfun(@class, varargin, 'uni', false);
typePass = ismember(types,'pgss');  
if (sum(typePass) > 1)
    error('PGSS:CONNECT','Only 1 pgss object can be interconnected');
end

% pass plant
pdG = varargin{typePass};
[ny,nu] = iosize(pdG);
u = pdG.u';
y = pdG.y';

% auxiliary model
I = ss(eye(nu)); 
I.u = u;        
I.y = sprintf('in_%6.0f',rand(1)*1e6);
O = ss(eye(ny)); 
O.u = sprintf('out_%6.0f',rand(1)*1e6); 
O.y = y;
pdG.u = I.y;
pdG.y = O.u;

% LTI part
inputs  = [O.u' varargin{end-1}];
outputs = [I.y' varargin{end}];
sys = varargin(1:end-2);

% removing the LPV model from the list of models
sys(typePass) = [];

% adding the new auxiliary models
sys = [sys, {I, O}];

% interconnecting all model but the LPV one but leaving the first inputs
% and output to do the LFT connection
M = connect(sys{:},inputs,outputs);

% now build the new interconnected model as the LFT between pdG and M
[no,ni] = iosize(M);
nw = ni - ny;
nz = no - nu;
[E,F1,F2,G1,G2,H11,H12,H21,H22] = parsysdata(M,[nz nw]);
if any(H11)
    error('PGSS:CONNECT','H11 must be 0')
end

A = pdG.A;
B = pdG.B;
C = pdG.C;
D = pdG.D;
n = nsys(pdG);

% independent term
ii = 1;
At(:,:,ii) = [A(:,:,ii)       B(:,:,ii)*G1;
              F1*C(:,:,ii)    E + F1*D(:,:,ii)*G1];
Bt(:,:,ii) = [B(:,:,ii)*H12;
              F1*D(:,:,ii)*H12 + F2];
Ct(:,:,ii) = [H21*C(:,:,ii), G2 + H21*D(:,:,ii)*G1];
Dt(:,:,ii) = H22 + H21*D(:,:,ii)*H12;
% remaining terms
for ii = 2:n
    At(:,:,ii) = [A(:,:,ii)       B(:,:,ii)*G1;
                  F1*C(:,:,ii)    F1*D(:,:,ii)*G1];
    Bt(:,:,ii) = [B(:,:,ii)*H12;
                  F1*D(:,:,ii)*H12];
    Ct(:,:,ii) = [H21*C(:,:,ii), H21*D(:,:,ii)*G1];
    Dt(:,:,ii) = H21*D(:,:,ii)*H12;
end

% resulting model
warning('off','Control:ltiobject:RepeatedChannelNames')

pdGint = pgss(At, Bt, Ct, Dt, pdG.parset, pdG.parfcn, ...
            'InputName', M.u(ny+1:end)', ...
            'OutputName', M.y(nu+1:end)');

warning('on','Control:ltiobject:RepeatedChannelNames')



    
