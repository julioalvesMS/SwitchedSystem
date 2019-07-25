classdef buck
    
    properties (Constant = true) 
        name = 'Buck';
        
        simulink = 'ideal_buck.slx'
        
        test_voltages = 10:10:90;
        
        single_voltage = 60;
        
        pwm_pid_kp = 0.5786;
        pwm_pid_ki = 14.24;
        pwm_pid_kd = 0.0119;
        
        reference_pid_kp = 5;
        reference_pid_ki = 140;
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
                0   0
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
            
            error = 1;
            
            upper = Vin;
            lower = -Vin;
        end
    end
end