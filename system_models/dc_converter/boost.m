classdef boost
    
    properties (Constant = true) 
        name = 'Boost';
        class_name = 'Boost';
        
        simulink = 'ideal_boost.slx'
        discrete_simulink = 'discrete_boost.slx'
        limit_cycle_simulink = 'limit_cycle_boost.slx'
        
        test_voltages = 80:20:200;
        
        single_voltage = 110;
        limit_cycle_voltage = 100;
        limit_cycle_gamma = [0 0.25];
        
        pwm_pid_kp = 0.0020;
        pwm_pid_ki = 0.2;
        
        pwm_pid_vc_vp = 0.316;
        pwm_pid_vc_vi = 3.23;
        pwm_pid_vc_cp = 0.0203;
        pwm_pid_vc_ci = 4.77;
        
        reference_pid_kp = 1;
        reference_pid_ki = 15.3;
        
        current_correction_gain = 2;
        
        operation_range_voltage_min = 70;
        operation_range_voltage_max = 120;
        
        % Gi = ((Ve/L)*s + (1-d)*Ie/(L*Co))/(s^2 - R*s/L + (1-d)^2/L*Co);
    end
    
    properties
        R;
        L;
        Ro;
        Co;
    end
    
    methods
        function converter = boost(R, Ro, Co, L)
            converter.R = R;
            converter.L = L;
            converter.Ro = Ro;
            converter.Co = Co;
        end
        
        function sys = get_sys(self)
        %SYS_BOOST Space State from Boost DC-DC converter
            
            A{1} = [
                -self.R/self.L  0
                0     -1/(self.Ro*self.Co)
            ];


            A{2} = [
                -self.R/self.L  -1/self.L
                 1/self.Co -1/(self.Ro*self.Co)
            ];

            B{1} = [1/self.L; 0];
            B{2} = B{1};

            C{1} = [0 1/sqrt(self.Ro)];
            C{2} = C{1};

            D{1} = 0;
            D{2} = D{1};

            Q{1} = [
                self.R   0
                0   3e1/self.Ro
            ];
            Q{2} = Q{1};

            E{1} = [
                1   0
                0   1/self.Ro
            ];
            E{2} = E{1};

            G{1} = [0;0];
            G{2} = G{1};

            H{1} = [0;0];
            H{2} = H{1};

            sys = gss(A, B, C, D, Q, E, G, H);
        end
        
        function [lower, upper] = get_pwm_control_limits(self)
            
            Rratio = self.R/self.Ro;
            
            lower = 0;
            upper = 1 - sqrt(Rratio);
        end
        
        function fnc = get_converter_Ie_fnc(self)
            fnc = @(Ve, Vin) (Vin/(2*self.R) - sqrt( Vin^2/(4*self.R^2) - Ve^2/(self.R*self.Ro) ) );
        end
    
        function [lower, upper] = get_reference_ve_limits(self, Vin)
            
            error = 1;
            
            upper = Vin/2 * sqrt(self.Ro/self.R) - error;
            lower = -Vin/2 * sqrt(self.Ro/self.R) + error;
        end
        
        function vi = get_reference_initial(self, Vs)
            vi = Vs*self.Ro/(self.R+self.Ro);
        end
    end
end
