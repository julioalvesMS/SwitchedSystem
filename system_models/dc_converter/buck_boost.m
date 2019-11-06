classdef buck_boost
    
    properties (Constant = true) 
        name = 'Buck-Boost';
        class_name = 'BuckBoost';
        
        simulink = 'ideal_buck_boost.slx'
        discrete_simulink = 'discrete_buck_boost.slx'
        
        test_voltages = 10:20:170;
        
        single_voltage = 80;
        
        pwm_pid_kp = 0.123;
        pwm_pid_ki = 28.9;
        pwm_pid_kd = 0.0;
        
        reference_pid_kp = 3;
        reference_pid_ki = 15.3;
    end
    
    properties
        R;
        L;
        Ro;
        Co;
    end
    
    methods
        function converter = buck_boost(R, Ro, Co, L)
            converter.R = R;
            converter.L = L;
            converter.Ro = Ro;
            converter.Co = Co;
        end
        
        function sys = get_sys(self)
        %SYS_BUCK_BOOST Space State from Buck-Boost DC-DC converter
            
            A{1} = [
                -self.R/self.L  0
                0     -1/(self.Ro*self.Co)
            ];


            A{2} = [
                -self.R/self.L  -1/self.L
                 1/self.Co -1/(self.Ro*self.Co)
            ];

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
            
            Rratio = self.R/self.Ro;
            
            lower = 0;
            upper = Rratio + 1 - sqrt(Rratio*(1 + Rratio));
        end
        
        function fnc = get_converter_Ie_fnc(self)
            fnc = @(Ve, Vin) (Vin/(2*self.R) - sqrt( Vin^2/(4*self.R^2) - Ve*(Ve+Vin)/(self.R*self.Ro) ) );
        end
        
        function [lower, upper] = get_reference_ve_limits(self, Vin)
            
            error = 1;
            
            upper = -Vin/2 * ( 1 - sqrt(1 + self.Ro/self.R)) - error;
            lower = -Vin/2 * ( 1 + sqrt(1 + self.Ro/self.R)) + error;
        end
    end
end

