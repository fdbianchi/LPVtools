function Go = uminus(Gi)

% UMINUS unary minus for LPV models.
%
%   MM = UMINUS(M) is invoked by MM = -M.
%
% See also pass

% fbianchi - 2020-02-20

Go = pass(Gi.A, Gi.B, -Gi.C, -Gi.D, Gi.parset,...
        'StateName',Gi.StateName,...
        'InputName',Gi.InputName,...
        'OutputName',Gi.OutputName,...
        'InputGroup',Gi.InputGroup,...
        'OutputGroup',Gi.OutputGroup);

