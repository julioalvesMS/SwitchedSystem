
%System matrices

        F{1} = [ -247.3498         0
         0   -4.5914];
        F{2} = [-247.3498 -504.7956
  444.4444   -4.5914];
g{1} = 1.0e+04*[3.3317
         0];
g{2} = [0 0]';
%System dimensions
n = 2;
N = 2;
%Lambda for finding xe
lambda = [0.6102    0.3898]';

T = 2.5e-5;% [s]

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
