function SystemDataBus = create_bus_LimitCycleDataBus(A, B, P, Q, N, xe, ell, rho, E, G, H)
%BUS_ADAPT Summary of this function goes here
%   Detailed explanation goes here

    elems(1)  = create_element_general(A, 'A');
    elems(2)  = create_element_general(B, 'B');
    elems(3)  = create_element_general(P, 'P');
    elems(4)  = create_element_general(Q, 'Q');
    elems(5)  = create_element_general(N, 'N');
    elems(6)  = create_element_general(xe, 'xe');
    elems(7)  = create_element_general(ell, 'ell');
    elems(8)  = create_element_general(rho, 'rho');
    elems(9)  = create_element_general(E, 'E');
    elems(10) = create_element_general(G, 'G');
    elems(11) = create_element_general(H, 'H');

    SystemDataBus = Simulink.Bus;
    SystemDataBus.Elements = elems;
end

