classdef buck_boost
    
    properties (Constant = true) 
        name = 'Buck-Boost';
        class_name = 'BuckBoost';
        
        simulink = 'ideal_buck_boost.slx'
        discrete_simulink = 'discrete_buck_boost.slx'
        limit_cycle_simulink = 'limit_cycle_buck_boost.slx'
        
        test_voltages = 10:20:170;
        
        single_voltage = 100;
        limit_cycle_voltage = 80;
        limit_cycle_gamma = [0 0.3];
        
%         pwm_pid_kp = 0.00483;
%         pwm_pid_ki = 0.147;
%         pwm_pid_kp = 0.00256;
%         pwm_pid_ki = 0.388;
        pwm_pid_kp = 0.00283;
        pwm_pid_ki = 0.312;
        
        pwm_pid_vc_vp = 0.316;
        pwm_pid_vc_vi = 3.23;
        pwm_pid_vc_cp = 0.0203;
        pwm_pid_vc_ci = 4.77;
        
        reference_pid_kp = 1;
        reference_pid_ki = 15.3;
        
        current_correction_gain = 0.8;
        
        operation_range_voltage_min = 0;
        operation_range_voltage_max = 120;
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
                self.R   0
                0       3e1/self.Ro
            ];
            Q{2} = Q{1};

            E{1} = [
                self.R   0
                0       1e3/self.Ro
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
        
        function vi = get_reference_initial(self, Vs)
            vi = 0;
        end
    end
end

