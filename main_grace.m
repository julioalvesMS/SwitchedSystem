%% Initial Setup
clear; clc; close all;
global A b P C xe Q x0 opt

addpath(genpath('functions'))
addpath(genpath('models'))

image_folder = 'images';

sim_theorem_1 = 'theorem_1.slx';


[~,~]=mkdir(image_folder);
image_folder = strcat(image_folder, '\');

%% Desired DC-DC converter to use
opt = 2
% Options can be found in the system directory:
%   buck
%   boost
%   buck_boost

circuit = boost;

%% System specifications

R  = 2; % [Ohm]
Ro = 50; % [Ohm]
Co = 470e-6; % [F]
L  = 500e-6; % [H]

Vs = 100; % [V]

x0 = [0; 0];


A{1} = [
    -R/L  0
    0     -1/(Ro*Co)
    ];


A{2} = [
    -R/L  -1/L
    1/Co -1/(Ro*Co)
    ];

b{1} = [1/L; 0]*Vs;
b{2} = b{1};

C{1} = [0 1/sqrt(Ro)];
C{2} = C{1};

D{1} = 0;
D{2} = D{1};

Q{1} = [
    0   0
    0   1/Ro
    ];
Q{2} = Q{1};

%% Common calculations

lamb1 = 0:0.1:1;
vet = [];
for i=1:length(lamb1)
    
    Al = lamb1(i)*A{1} + (1-lamb1(i))*A{2};
    bl = lamb1(i)*b{1} + (1-lamb1(i))*b{2};
    
    if (max(real(eig(Al))) <0)
        
        
        xe = -inv(Al)*bl;
        
        vet = [vet;lamb1(i), xe(2)];
    end
end

plot(vet(:,1), vet(:,2));

lambda = [0.2164; 1-0.2164];

        Alamb = lambda(1)*A{1}+(1-lambda(1))*A{2};
        blamb = lambda(1)*b{1}+(1-lambda(1))*b{2};
        
xe = -Alamb\blamb;

switch(opt)
    case 1
        Qlamb = lambda(1)*Q{1}+(1-lambda(1))*Q{2};
        P = teorema1(Alamb, Qlamb);
    case 2
        P = teorema2(A, Q);
end
sim('simu_theo1_converter.slx',[0,0.06])

figure
plot(x(:,1), x(:,2:end));

figure
plot(x(:,2), x(:,3));

figure
plot(sig(:,1), sig(:,2));

figure
plot(custo(:,1), custo(:,2));


function xdata = switching_rule(u)
    global A b P Q C xe opt
    
    x = u;
    N = size(A,2);
    qsi = x-xe;
    
    switch(opt)
        case 1
            for i=1:N
                v(i) = qsi'*(Q{i}*qsi + 2*P*(A{i}*x + b{i}));
            end
        case 2
            
            for i=1:N
                v(i) = qsi'*P*(A{i}*xe + b{i});
            end
    end
    [~,idx] = min(v);
    xdot = A{idx}*x+b{idx};
    ze = C{idx}*qsi;

    xdata = [idx;ze'*ze;xdot];
end

function Po = teorema1(Alamb, Qlamb)
global x0 xe
    nx = size(Alamb,1);

    % Descreve a LMI a ser projetada
    setlmis([])

    % declaracao de variveis
    P = lmivar(1,[nx 1]);
    
    ct = newlmi;
    lmiterm([ct,1,1,P],Alamb',1,'s')
    lmiterm([ct,1,1,0],Qlamb)

    
    ct = newlmi;
    lmiterm([-ct,1,1,P],1,1)
    
    lmisys = getlmis;


    % Declaracao funcao objetivo
    options = [1e-4,2000,1e9,200,0];
    %===============

    np = decnbr(lmisys);

    c = zeros(np,1);
    for i=1:np
        Pi = defcx(lmisys,i,P);
        c(i) = (x0-xe)'*Pi*(x0-xe);
    end

    [copt,xopt] = mincx(lmisys,c,options);


    if (isempty(copt))
        Po = [];
        return
    end

    Po = dec2mat(lmisys,xopt,P);
end

function Po = teorema2(A, Q)
global x0 xe
    N = size(A, 2);
    nx = size(A{1},1);

    % Descreve a LMI a ser projetada
    setlmis([])

    % declaracao de variveis
    P = lmivar(1,[nx 1]);
    
    for i=1:N
        ct = newlmi;
        lmiterm([ct,1,1,P],A{i}',1,'s')
        lmiterm([ct,1,1,0],Q{i})
        
    end
    
    ct = newlmi;
    lmiterm([-ct,1,1,P],1,1)
    
    lmisys = getlmis;


    % Declaracao funcao objetivo
    options = [1e-4,2000,1e9,200,0];
    %===============

    np = decnbr(lmisys);

    c = zeros(np,1);
    for i=1:np
        Pi = defcx(lmisys,i,P);
        c(i) = (x0-xe)'*Pi*(x0-xe);
    end

    [copt,xopt] = mincx(lmisys,c,options);


    if (isempty(copt))
        Po = [];
        return
    end

    Po = dec2mat(lmisys,xopt,P);
end