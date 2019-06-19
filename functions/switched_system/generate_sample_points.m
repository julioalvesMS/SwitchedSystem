function [test_lambdas,equilibrium] = generate_sample_points(sys)
%GENERATE_SAMPLE_POINTS Summary of this function goes here
%   Detailed explanation goes here

    % Sample points
    test_lambdas = generate_lambda_2d(0.001);

    % Calculate the equilibrium point for each lambda sample
    for i=size(test_lambdas, 1):-1:1
        [Al, Bl, ~] = calc_sys_lambda(sys, test_lambdas(i,:));
        xe = -Al\Bl*sys.U;

        equilibrium(i,:) = xe;
    end
end

