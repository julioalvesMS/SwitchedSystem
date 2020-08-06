function lambdas = generate_lambda_voltage(sys, voltages)
%GENERATE_LAMBDA_VOLTAGE Summary of this function goes here
%   Detailed explanation goes here

    [sample_lambdas, equilibrium] = generate_sample_points(sys);

    % To be able to interpolate the data, we will use the only the
    % first lambda and the output voltage
    lamb1 = sample_lambdas(:,2);
    Vlamb = equilibrium(:,2);

    ns = length(voltages);
    lambdas = zeros(ns, sys.N);
    for i=1:ns
        % Only use the lambdas before the inflection point from the
        % voltage
        [~, inflection] = max(Vlamb);

        % interpolate the data to discover the needed lambda to reach
        % the desired voltage
        lamb = interp1(Vlamb(1:inflection), lamb1(1:inflection), voltages(i));
        
        if sys.N == 3
            lambdas(i,:) = [0 lamb 1-lamb];
        else
            lambdas(i,:) = [1-lamb lamb];
        end
    end

end

