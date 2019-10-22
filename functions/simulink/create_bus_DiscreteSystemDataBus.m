function SystemDataBus = create_bus_DiscreteSystemDataBus(A, B, P, Q, N, h, L)
%BUS_ADAPT Summary of this function goes here
%   Detailed explanation goes here

    elems(1) = create_element_general(A, 'A');
    elems(2) = create_element_general(B, 'B');
    elems(3) = create_element_general(P, 'P');
    elems(4) = create_element_general(Q, 'Q');
    elems(5) = create_element_general(N, 'N');
    elems(6) = create_element_general(h, 'h');
    elems(7) = create_element_general(L, 'L');

    SystemDataBus = Simulink.Bus;
    SystemDataBus.Elements = elems;
end

