
%System matrices

        F{1} = [-3 -6 3; 2 2 -3;  1.6  0 -2];
        F{2} = [1 3 3 ;-0.2 -3 -3; 0 0 -2];
g{1} = [.5 0 0]';
g{2} = [0 0 .5]';
%System dimensions
n=3;
N=2;
%Lambda for finding xe
lambda = [0.56;0.44];

T=1;% [s]

%Pi matrices
p = 0;
q =  0
MpiT1 = [p  1-q;
       1-p  q];

MpiC1 = [p  1-q;
       1-p  q];
%Betas for Theorem 1
betas = [1 .1684];
betasOF = [1 1];
p = 0;
q = 0;
MpiOF =  [p  1-q;
       1-p  q];

%[0 0.5450 1 1.5250 0.0331 1.0000e-04]


C{1} = [1 0 0; 0 0 1];%eye(n);
C{2} = [1 0 0; 0 0 1];%eye(n);
