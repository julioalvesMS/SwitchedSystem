classdef buck_boost_non_inverting
    
    properties (Constant = true) 
        name = 'Buck-Boost Bidirecional';
        
        simulink = 'ideal_buck_boost_non_inverter.slx'
        
        test_voltages = 10:20:170;
        
        single_voltage = 120;
        
        pwm_pid_kp = 0.5;
        pwm_pid_ki = 0.1;
        pwm_pid_kd = 0;

%         pwm_pid_kp = 5e-2;
%         pwm_pid_ki = 100e-2;
        
        reference_pid_kp = 5e-1;
        reference_pid_ki = 100e-1;
    end
    
    properties
        R;
        L;
        Ro;
        Co;
    end
    
    methods
        function converter = buck_boost_non_inverting(R, Ro, Co, L)
            converter.R = R;
            converter.L = L;
            converter.Ro = Ro;
            converter.Co = Co;
        end
        
        function sys = get_sys(self)
        %SYS_BUCK_BOOST_NON_INVERTING Space State from Buck-Boost DC-DC converter
            
            A{1} = [
                -self.R/self.L  -1/self.L
                 1/self.Co -1/(self.Ro*self.Co)
            ];

            A{2} = A{1};
            
            A{3} = [
                -self.R/self.L  0
                0     -1/(self.Ro*self.Co)
            ];

            B{1} = [1/self.L; 0];
            B{2} = [0; 0];
            B{3} = B{1};

            C{1} = [0 1/sqrt(self.Ro)];
            C{2} = C{1};
            C{3} = C{1};

            D{1} = 0;
            D{2} = D{1};
            D{3} = D{1};

            Q{1} = [
                0   0
                0   1/self.Ro
            ];
            Q{2} = Q{1};
            Q{3} = Q{1};

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
            
            upper = Vin/2 * sqrt(self.Ro/self.R) - error;
            lower = 0;
        end
    end
end

