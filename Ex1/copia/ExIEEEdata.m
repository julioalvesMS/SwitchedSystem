%System matrices
F{1} = [0 1 0; 0 0 1; -1 -1 -1];
F{2} = [0 1 0 ;0 0 1; 0 -1 -1];
g{1} = [1 0 0]';
g{2} = [0 1 0]';
%System dimensions
n=3;
N=2;
%Lambda for finding xe
lambda = [0.6;0.4];

T=1;% [s]

%Pi matrices
p = 0;
q = 0;
MpiT1 = [p  1-q;
       1-p  q];

MpiC1 = [p  1-q;
       1-p  q];
%Betas for Theorem 1
betas = [1 0.4720];
betasOF = [1.3316 1.0184];
p = 0;
q = 0.3140;
MpiOF =  [p  1-q;
       1-p  q];

%[0 0.5450 1 1.5250 0.0331 1.0000e-04]


C{1} = [1 0 0; 0 0 1];%eye(n);
C{2} = [1 0 0; 0 0 1];%eye(n);