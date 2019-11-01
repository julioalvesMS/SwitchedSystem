    void Boost::GetP(double P[SYSTEM_ORDER][SYSTEM_ORDER])
    {
        switch(controlStrategy)
        {
        case CS_CONTINUOUS_THEOREM_1:
            //
            // Boost Converter - Rule 1
            //
            P[0][0] = 4.48515e-05;
            P[0][1] = 2.22111e-05;
            P[1][0] = 2.22111e-05;
            P[1][1] = 7.63392e-05;
            break;

        case CS_CONTINUOUS_THEOREM_2:
            //
            // Boost Converter - Rule 2
            //
            P[0][0] = 0.000971581;
            P[0][1] = 4.96527e-07;
            P[1][0] = 4.96527e-07;
            P[1][1] = 0.001125;
            break;

        case CS_DISCRETE_THEOREM_1:
            //
            // Boost Converter - Discrete Rule 1
            //
            P[0][0] = 2.91295e-06;
            P[0][1] = 1.6538e-06;
            P[1][0] = 1.6538e-06;
            P[1][1] = 8.80047e-06;
            break;

        default:
            break;
        }
    }


    void Boost::GetH(double h[SYSTEM_ORDER])
    {
        switch(controlStrategy)
        {
        case CS_DISCRETE_THEOREM_1:
            //
            // Boost Converter - Discrete Rule 1
            //
            h[0] = 6.01166e-06;
            h[1] = 9.7411e-06;
            break;

        default:
            break;
        }
    }


    double Buck::GetD(double P[SYSTEM_ORDER][SYSTEM_ORDER], double h[SYSTEM_ORDER])
    {
        double d = 0;

        switch(controlStrategy)
        {
        case CS_DISCRETE_THEOREM_1:
            //
            // Boost Converter - Discrete Rule 1
            //
            d = 1.75003e-05;
            break;

        default:
            break;
        }

        return d;
    }


    void DefineDiscreteSystem()
    {
        SubSystem* subSys;

        //
        // =============== Subsystem 1 ===============
        //
        subSys = &(discreteSystem.subSystems[0]);

        //
        // Subsystem 1 -- Matrix A
        //
        subSys->A[0][0] = 0.989817;
        subSys->A[0][1] = 0;
        subSys->A[1][0] = 0;
        subSys->A[1][1] = 0.99977;
        //
        // Subsystem 1 -- Matrix B
        //
        subSys->L[0] = 0.0247173;
        subSys->L[1] = -0.000526048;
        //
        // Subsystem 1 -- Matrix Q
        //
        subSys->Q[0][0] = 0.01;
        subSys->Q[0][1] = 0;
        subSys->Q[1][0] = 0;
        subSys->Q[1][1] = 0.0103306;


        //
        // =============== Subsystem 2 ===============
        //
        subSys = &(discreteSystem.subSystems[0]);

        //
        // Subsystem 2 -- Matrix A
        //
        subSys->A[0][0] = 0.989534;
        subSys->A[0][1] = -0.0254527;
        subSys->A[1][0] = 0.0221042;
        subSys->A[1][1] = 0.999487;
        //
        // Subsystem 2 -- Matrix B
        //
        subSys->L[0] = -0.0336363;
        subSys->L[1] = 0.000715866;
        //
        // Subsystem 2 -- Matrix Q
        //
        subSys->Q[0][0] = 0.01;
        subSys->Q[0][1] = 0;
        subSys->Q[1][0] = 0;
        subSys->Q[1][1] = 0.0103306;

    }


