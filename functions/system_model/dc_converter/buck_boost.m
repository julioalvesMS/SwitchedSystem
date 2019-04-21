classdef buck_boost
    
    properties (Constant = true) 
        name = 'Buck-Boost';
        
        simulink = 'ideal_buck_boost.slx'
    end
    
    methods(Static)
        function sys = get_sys(R, Ro, Co, L)
        %SYS_BUCK_BOOST Space State from Buck-Boost DC-DC converter
            
            A{1} = [
                -R/L  0
                0     -1/(Ro*Co)
            ];


            A{2} = [
                -R/L  -1/L
                 1/Co -1/(Ro*Co)
            ];

            B{1} = [1/L; 0];
            B{2} = [0; 0];

            C{1} = [0 1/sqrt(Ro)];
            C{2} = C{1};

            D{1} = 0;
            D{2} = D{1};

            Q{1} = [
                0   0
                0   1/Ro
            ];
            Q{2} = Q{1};

            sys = gss(A, B, C, D, Q);
        end
    end
end

