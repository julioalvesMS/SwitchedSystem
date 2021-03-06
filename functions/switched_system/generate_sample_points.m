function [test_lambdas,equilibrium] = generate_sample_points(sys)
%GENERATE_SAMPLE_POINTS Summary of this function goes here
%   Detailed explanation goes here

    % Sample points
    test_lambdas = generate_lambda_2d(0.0001);
    
    Nl = size(test_lambdas, 1);
    
    if sys.N == 3
        lambdas = [zeros(Nl,1) test_lambdas];
    else
        lambdas = test_lambdas;
    end

    % Calculate the equilibrium point for each lambda sample
    for i=Nl:-1:1
        [Al, Bl, ~] = calc_sys_lambda(sys, lambdas(i,:));
        if sys.discrete
            xe = -(Al-eye(2))\Bl*sys.U;
        else
            xe = -Al\Bl*sys.U;
        end
        equilibrium(i,:) = xe;
    end
end

