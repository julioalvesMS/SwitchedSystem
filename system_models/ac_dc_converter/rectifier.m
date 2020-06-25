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
            % Subsistemas 1-2: Quando corrente na terceira fase é maior ou
            % igual a zero
            % Subsistemas 3-4: Quando a corrente na terceira fase é menor
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
            B{3} = B{1};
            B{4} = B{2};
            

            C{1} = [1   0   0];
            C{2} = C{1};
            C{3} = C{1};
            C{4} = C{1};

            D{1} = [0   0   0];
            D{2} = D{1};
            D{3} = D{1};
            D{4} = D{1};

            Q{1} = [
                1e0     0       0
                0       0       0
                0       0       0
            ];
            Q{2} = Q{1};
            Q{3} = Q{1};
            Q{4} = Q{1};

            E{1} = eye(3);
            E{2} = E{1};
            E{3} = E{1};
            E{4} = E{1};

            G{1} = 0;
            G{2} = G{1};
            G{3} = G{1};
            G{4} = G{1};

            H{1} = 0;
            H{2} = H{1};
            H{3} = H{1};
            H{4} = H{1};

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