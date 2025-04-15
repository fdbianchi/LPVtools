
**LPVtool** is a small toolbox aimed to design gain-scheduled Linear Parameter Varying (LPV) controllers. It provides tools for modelling, setting specifications, synthesizing controllers and evaluating results. 

# Getting started

**LPVtool** implements several synthesis procedures for linear parameter varying (LPV) systems. The toolbox provides a set of objects for modelling and analyzing LPV systems and the function !lpvsyn! for designing gain-scheduled LPV controllers. The basic use is shown in the following code:

```matlab
% definition of the parameter set
pv = pset.Box([0 1;0 2]);
% definition of the LPV model (affine in this case)
A(:,:,1) = [-2 1; -1 -2];
A(:,:,2) = [-1.3 0.1; -0.1 -1.3];
A(:,:,3) = [-0.5 0; 0 -0.5];
B        = [1;1];
C(:,:,1) = [10 10];
C(:,:,2) = [2 2];
C(:,:,3) = [5 5];
D = zeros(1,1);
pdG = pass(A,B,C,D,pv);
pdG.u = 'u';    pdG.y = 'y';
% definition of the augmented plant
sb    = sumblk('e = r - y');
pdGau = connect(pdG,sb,{'r','u'},{'e','u','e'});
% weights
We = tf([0.2 20],[1 1]);
Wu = tf([0.04 0.1],[0.004 1]);
Wout = append(W1,W2,1);
% augmented plant + weights
Gaw = wout*Gau;
% constraints
const(1) = synConst.Gain(1,1:2);
const(2) = synConst.Poles('MaxFreq',1000);
% synthesis
pdK = lpvsyn(pdGaw,3,2,const1);
```

In the folder **demos**, it can be found several examples of the use of the toolbox in more complex applications. A documentation report can be found in the folder **doc**.
