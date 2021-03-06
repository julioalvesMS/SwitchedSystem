function lambda = generate_lambda_2d(step, last)
%GENERATE_LAMBDA Generates generic a vector for a 2 dimensional system
    if nargin == 0
        step = 0.05;
    end
    
    if nargin < 2
        last = 1;
    end
    
    sequence = 0:step:last;
    
    lambda = zeros(length(sequence) ,2);
    lambda(:,1) = sequence;
    lambda(:,2) = 1-sequence;
end
