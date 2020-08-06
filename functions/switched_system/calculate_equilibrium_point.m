function [xe] = calculate_equilibrium_point(circuit, Vin, Vref)
%CALCULATE_EQUILIBRIUM_POINT Summary of this function goes here
%   Detailed explanation goes here
    fnc = circuit.get_converter_Ie_fnc();
    
    xe(2,1) = Vref;
    xe(1,1) = fnc(Vref, Vin);
end

