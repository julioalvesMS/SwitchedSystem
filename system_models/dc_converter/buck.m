classdef buck
    
    properties (Constant = true) 
        name = 'Buck';
        
        simulink = 'ideal_buck.slx'
        discrete_simulink = 'discrete_buck.slx'
        
        test_voltages = 5:5:50;
        
        single_voltage = 60;
        
        pwm_pid_kp = 0.15;
        pwm_pid_ki = 15;
        pwm_pid_kd = 0;

%         pwm_pid_kp = 5e-2;
%         pwm_pid_ki = 100e-2;
        
        reference_pid_kp = 5e-1;
        reference_pid_ki = 100e0;
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
                1e-2   0
                0   1/self.Ro
            ];
            Q{2} = Q{1};

            sys = gss(A, B, C, D, Q);
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
    end
end