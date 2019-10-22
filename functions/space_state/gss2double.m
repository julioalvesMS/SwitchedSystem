function [A, B, C, D, Q, L] = gss2double(sys)
%GSS2DOUBLE General Space State system
    A = catcell(sys.A);
    B = catcell(sys.B);
    C = catcell(sys.C);
    D = catcell(sys.D);
    Q = catcell(sys.Q);
    
    if isfield(sys, 'L')
        L = catcell(sys.L);
    end
end

function matrix = catcell(cell_array)
    for i=length(cell_array):-1:1
        matrix(:,:,i) = cell_array{i};
    end
end
