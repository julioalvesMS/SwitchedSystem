function [K,M] = project_state_feedback_h2(sys, T)
    
    z = tf('z', T);
    dsys = c2d(sys, T, 'tustin');
    [~,K,~] = calc_P_norm_h2(dsys);

    Ad = dsys.A;
    Bd = dsys.B;
    Cd = sys.C;
    Dd = sys.D;
    I = eye(size(Ad));
    S = (Cd-Dd*K)/(z*I-(Ad-Bd*K))*Bd+Dd;

    M = 1/(evalfr(S, 1));
end

