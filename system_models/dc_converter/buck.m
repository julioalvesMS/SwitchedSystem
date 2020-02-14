classdef buck
    
    properties (Constant = true) 
        name = 'Buck';
        class_name = 'Buck';
        
        simulink = 'ideal_buck.slx'
        discrete_simulink = 'discrete_buck.slx'
        limit_cycle_simulink = 'limit_cycle_buck.slx'
        
        test_voltages = 5:5:50;
        
        single_voltage = 30;
        limit_cycle_voltage = 30;
        limit_cycle_gamma = [0 3];
        limit_cycle_kappa = 10;
        
%         pwm_pid_kp = 0.00919;
%         pwm_pid_ki = 3.28;
%         pwm_pid_kp = 0.00675;
%         pwm_pid_ki = 3.85;
        pwm_pid_kp = 0.00676;
        pwm_pid_ki = 3.13;
        
        pwm_pid_vc_vp = 0.316;
        pwm_pid_vc_vi = 3.23;
        pwm_pid_vc_cp = 0.0203;
        pwm_pid_vc_ci = 4.77;

%         pwm_pid_kp = 5e-2;
%         pwm_pid_ki = 100e-2;
        
        reference_pid_kp = 1;
        reference_pid_ki = 100;
        
        current_correction_gain = 2;
    end
    
    properties
        R;
        L;
        Ro;
        Co;
    end
    
    methods
        function converter = buck(R, Ro, Co, L)
            converter.R = R;
            converter.L = L;
            converter.Ro = Ro;
            converter.Co = Co;
        end
        
        function sys = get_sys(self)
        %SYS_BUCK Space State from Buck DC-DC converter
            
            A{1} = [
                -self.R/self.L  -1/self.L
                 1/self.Co -1/(self.Ro*self.Co)
            ];

            A{2} = A{1};

            B{1} = [1/self.L; 0];
            B{2} = [0; 0];

            C{1} = [0 1/sqrt(self.Ro)];
            C{2} = C{1};

            D{1} = 0;
            D{2} = D{1};

            Q{1} = [
                1e-1   0
                0   1/self.Ro
            ];
            Q{2} = Q{1};

            E{1} = [
                0   0
                0   0.1/self.Ro
            ];
            E{2} = E{1};

            G{1} = [0;0];
            G{2} = G{1};

            H{1} = [0;0];
            H{2} = H{1};

            sys = gss(A, B, C, D, Q, E, G, H);
        end
        
        function [lower, upper] = get_pwm_control_limits(self)
            
            lower = 0;
            upper = 1;
        end
        
        function fnc = get_converter_Ie_fnc(self)
            fnc = @(Ve, Vin) (Ve/self.Ro);
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