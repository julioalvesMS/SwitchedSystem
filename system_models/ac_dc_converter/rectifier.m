classdef rectifier
    
    properties (Constant = true) 
        name = 'Buck';
        class_name = 'Buck';
        
        simulink = 'sim_rectifier.slx'
        
        single_voltage = 400;
        limit_cycle_voltage = 30;
        limit_cycle_gamma = [0 3];
        limit_cycle_kappa = 10;
        
        pwm_pid_kp = 0.00676;
        pwm_pid_ki = 3.13;
        
        pwm_pid_vc_vp = 0.316;
        pwm_pid_vc_vi = 3.23;
        pwm_pid_vc_cp = 0.0203;
        pwm_pid_vc_ci = 4.77;
    end
    
    properties
        R;
        L;
        Co;
    end
    
    methods
        function converter = rectifier(R, L, Co)
            converter.R = R;
            converter.L = L;
            converter.Co = Co;
        end
        
        function sys = get_sys(self)
            % Subsistemas 1-6: Quando corrente na terceira fase é maior ou
            % igual a zero
            % Subsistemas 7-12: Quando a corrente na terceira fase é menor
            % que zero
            
            % ================= Matrizes A =================
            A{1} = [
                0               0               -1/self.Co
                -1/(3*self.L)   -self.R/self.L  0
                2/(3*self.L)    0               -self.R/self.L
            ];
            A{2} = [
                0               1/self.Co       1/self.Co
                -1/(3*self.L)   -self.R/self.L  0
                -1/(3*self.L)   0               -self.R/self.L
            ];
            A{3} = [
                0               -1/self.Co      0
                2/(3*self.L)    -self.R/self.L  0
                -1/(3*self.L)   0               -self.R/self.L
            ];
            A{4} = [
                0               1/self.Co       1/self.Co
                -1/(3*self.L)   -self.R/self.L  0
                -1/(3*self.L)   0               -self.R/self.L
            ];
            A{5} = [
                0               -1/self.Co      0
                2/(3*self.L)    -self.R/self.L  0
                -1/(3*self.L)   0               -self.R/self.L
            ];
            A{6} = [
                0               -0              -1/self.Co
                -1/(3*self.L)   -self.R/self.L  0
                2/(3*self.L)    0               -self.R/self.L
            ];
        
            A{7} = [
                0               -2/self.Co      -1/self.Co
                -2/(3*self.L)   -self.R/self.L  0
                1/(3*self.L)    0               -self.R/self.L
            ];
            A{8} = [
                0               1/self.Co       -1/self.Co
                -2/(3*self.L)   -self.R/self.L  0
                1/(3*self.L)    0               -self.R/self.L
            ];
            A{9} = [
                0               -1/self.Co      -2/self.Co
                1/(3*self.L)    -self.R/self.L  0
                -2/(3*self.L)   0               -self.R/self.L
            ];
            A{10} = [
                0               -1/self.Co      1/self.Co
                1/(3*self.L)    -self.R/self.L  0
                -2/(3*self.L)   0               -self.R/self.L
            ];
            A{11} = [
                0               -1/self.Co      -2/self.Co
                1/(3*self.L)    -self.R/self.L  0
                1/(3*self.L)    0               -self.R/self.L
            ];
            A{12} = [
                0               -2/self.Co      -1/self.Co
                1/(3*self.L)    -self.R/self.L  0
                1/(3*self.L)    0               -self.R/self.L
            ];

            
            % ================= Matrizes B =================
            B{1} = [
                0               0               0
                2/(3*self.L)    -1/(3*self.L)   -1/(3*self.L)
                -1/(3*self.L)    2/(3*self.L)   -1/(3*self.L)
            ];
            B{2} = [
                0               0               0
                2/(3*self.L)    -1/(3*self.L)   -1/(3*self.L)
                -1/(3*self.L)   2/(3*self.L)    -1/(3*self.L)
            ];
            B{3} = [
                0               0               0
                2/(3*self.L)    -1/(3*self.L)   -1/(3*self.L)
                -1/(3*self.L)   2/(3*self.L)    -1/(3*self.L)
            ];
            B{4} = [
                0               0               0
                2/(3*self.L)    -1/(3*self.L)   -1/(3*self.L)
                -1/(3*self.L)   2/(3*self.L)    -1/(3*self.L)
            ];
            B{5} = [
                0               0               0
                2/(3*self.L)    -1/(3*self.L)   -1/(3*self.L)
                -1/(3*self.L)   2/(3*self.L)    -1/(3*self.L)
            ];
            B{6} = [
                0               0               0
                2/(3*self.L)    -1/(3*self.L)   -1/(3*self.L)
                -1/(3*self.L)   2/(3*self.L)    -1/(3*self.L)
            ];
            B{7} = B{1};
            B{8} = B{2};
            B{9} = B{3};
            B{10} = B{4};
            B{11} = B{5};
            B{12} = B{6};
            

            C{1} = [1   0   0];
            C{2} = C{1};
            C{3} = C{1};
            C{4} = C{1};
            C{5} = C{1};
            C{6} = C{1};
            C{7} = C{1};
            C{8} = C{1};
            C{9} = C{1};
            C{10} = C{1};
            C{11} = C{1};
            C{12} = C{1};

            D{1} = [0   0   0];
            D{2} = D{1};
            D{3} = D{1};
            D{4} = D{1};
            D{5} = D{1};
            D{6} = D{1};
            D{7} = D{1};
            D{8} = D{1};
            D{9} = D{1};
            D{10} = D{1};
            D{11} = D{1};
            D{12} = D{1};

            Q{1} = [
                1e0     0       0
                0       0       0
                0       0       0
            ];
            Q{2} = Q{1};
            Q{3} = Q{1};
            Q{4} = Q{1};
            Q{5} = Q{1};
            Q{6} = Q{1};
            Q{7} = Q{1};
            Q{8} = Q{1};
            Q{9} = Q{1};
            Q{10} = Q{1};
            Q{11} = Q{1};
            Q{12} = Q{1};

            E{1} = eye(3);
            E{2} = E{1};
            E{3} = E{1};
            E{4} = E{1};
            E{5} = E{1};
            E{6} = E{1};

            G{1} = 0;
            G{2} = G{1};
            G{3} = G{1};
            G{4} = G{1};
            G{5} = G{1};
            G{6} = G{1};

            H{1} = 0;
            H{2} = H{1};
            H{3} = H{1};
            H{4} = H{1};
            H{5} = H{1};
            H{6} = H{1};

            sys = gss(A, B, C, D, Q, E, G, H);
        end
        
        function [lower, upper] = get_pwm_control_limits(self)
            
            lower = 0;
            upper = 1;
        end
        
        function [lower, upper] = get_reference_ve_limits(self, Vin)
            upper = Vin;
            lower = 0;
        end
        
        function vi = get_reference_initial(self, Vs)
            vi = 0;
        end
    end
end