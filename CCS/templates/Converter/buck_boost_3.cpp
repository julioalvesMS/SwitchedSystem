#include <src/Converter/buck_boost_3.h>

using namespace SwitchedSystem;
using namespace BaseConverter;
using namespace Controller;

extern ControlStrategy controlStrategy;

namespace ConverterBuckBoost3
{
    static System system;
    static System discreteSystem;


    System* BuckBoost3::GetSys()
    {
        DefineSystem();

        system.N = 3;

        return &(system);
    }


    System* BuckBoost3::GetDiscreteSys()
    {
        DefineDiscreteSystem();

        discreteSystem.N = 3;

        return &(discreteSystem);
    }


    void BuckBoost3::GetP(double P[SYSTEM_ORDER][SYSTEM_ORDER])
    {
        switch(controlStrategy)
        {
        case CS_CONTINUOUS_THEOREM_1:
            //
            // Buck-Boost Converter - Rule 1
            //
            P[0][0] = 4.4977e-05;
            P[0][1] = 1.9091e-05;
            P[1][0] = 1.9091e-05;
            P[1][1] = 6.9922e-05;
            break;

        case CS_CONTINUOUS_THEOREM_2:
            //
            // Buck-Boost Converter - Rule 2
            //
            P[0][0] = 9.7491e-04;
            P[0][1] = 7.8654e-07;
            P[1][0] = 7.8654e-07;
            P[1][1] = 0.0011;
            break;

        case CS_DISCRETE_THEOREM_1:
            //
            // Buck-Boost Converter - Discrete Rule 1
            //
            P[0][0] = 3.7319e-06;
            P[0][1] = 2.4064e-06;
            P[1][0] = 2.4064e-06;
            P[1][1] = 1.0160e-05;
            break;

        default:
            break;
        }
    }

    void BuckBoost3::GetH(double h[SYSTEM_ORDER])
    {
        switch(controlStrategy)
        {
        case CS_DISCRETE_THEOREM_1:
            //
            // Buck-Boost Converter - Discrete Rule 1
            //
            h[0] = 5.9279e-06;
            h[1] = 9.3472e-06;
            break;

        default:
            break;
        }
    }


    double BuckBoost3::GetD(double P[SYSTEM_ORDER][SYSTEM_ORDER], double h[SYSTEM_ORDER])
    {
        double d = 0;

        switch(controlStrategy)
        {
        case CS_DISCRETE_THEOREM_1:
            //
            // Buck-Boost Converter - Discrete Rule 1
            //
            d = 1.2962e-05;
            break;

        default:
            break;
        }

        return d;
    }


    int BuckBoost3::SubSystem2SwitchState(int SubSystem)
    {
        return SubSystem;
    }


    void DefineSystem()
    {
        SubSystem* subSys;

        //
        // =============== Subsystem 1 ===============
        //
        subSys = &(system.subSystems[0]);

        //
        // Subsystem 1 -- Matrix A
        //
        subSys->A[0][0] = -R/L;
        subSys->A[0][1] = -1/L;
        subSys->A[1][0] =  1/Co;
        subSys->A[1][1] = -1/(Ro*Co);
        //
        // Subsystem 1 -- Matrix B
        //
        subSys->B[0] = 1/L;
        subSys->B[1] = 0;
        //
        // Subsystem 1 -- Matrix Q
        //
        subSys->Q[0][0] = 0;
        subSys->Q[0][1] = 0;
        subSys->Q[1][0] = 0;
        subSys->Q[1][1] = 1/Ro;


        //
        // =============== Subsystem 2 ===============
        //
        subSys = &(system.subSystems[1]);

        //
        // Subsystem 2 -- Matrix A
        //
        subSys->A[0][0] = -R/L;
        subSys->A[0][1] = -1/L;
        subSys->A[1][0] =  1/Co;
        subSys->A[1][1] = -1/(Ro*Co);
        //
        // Subsystem 2 -- Matrix B
        //
        subSys->B[0] = 0;
        subSys->B[1] = 0;
        //
        // Subsystem 2 -- Matrix Q
        //
        subSys->Q[0][0] = 0;
        subSys->Q[0][1] = 0;
        subSys->Q[1][0] = 0;
        subSys->Q[1][1] = 1/Ro;


        //
        // =============== Subsystem 3 ===============
        //
        subSys = &(system.subSystems[2]);

        //
        // Subsystem 1 -- Matrix A
        //
        subSys->A[0][0] = -R/L;
        subSys->A[0][1] =  0;
        subSys->A[1][0] =  0;
        subSys->A[1][1] = -1/(Ro*Co);
        //
        // Subsystem 1 -- Matrix B
        //
        subSys->B[0] = 1/L;
        subSys->B[1] = 0;
        //
        // Subsystem 1 -- Matrix Q
        //
        subSys->Q[0][0] = 0;
        subSys->Q[0][1] = 0;
        subSys->Q[1][0] = 0;
        subSys->Q[1][1] = 1/Ro;


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
        subSys->A[0][0] =  0.9786;
        subSys->A[0][1] = -0.0506;
        subSys->A[1][0] =  0.0440;
        subSys->A[1][1] =  0.9984;
        //
        // Subsystem 1 -- Matrix B
        //
        subSys->L[0] = 0.0193;
        subSys->L[1] = 4.3142e-04;
        //
        // Subsystem 1 -- Matrix Q
        //
        subSys->Q[0][0] = 0;
        subSys->Q[0][1] = 0;
        subSys->Q[1][0] = 0;
        subSys->Q[1][1] = 1/Ro;


        //
        // =============== Subsystem 2 ===============
        //
        subSys = &(discreteSystem.subSystems[1]);

        //
        // Subsystem 2 -- Matrix A
        //
        subSys->A[0][0] =  0.9786;
        subSys->A[0][1] = -0.0506;
        subSys->A[1][0] =  0.0440;
        subSys->A[1][1] =  0.9984;
        //
        // Subsystem 2 -- Matrix B
        //
        subSys->L[0] = -0.0313;
        subSys->L[1] = -6.9774e-04;
        //
        // Subsystem 2 -- Matrix Q
        //
        subSys->Q[0][0] = 0;
        subSys->Q[0][1] = 0;
        subSys->Q[1][0] = 0;
        subSys->Q[1][1] = 1/Ro;


        //
        // =============== Subsystem 3 ===============
        //
        subSys = &(discreteSystem.subSystems[2]);

        //
        // Subsystem 1 -- Matrix A
        //
        subSys->A[0][0] =  0.9786;
        subSys->A[0][1] = -0.0506;
        subSys->A[1][0] =  0.0440;
        subSys->A[1][1] =  0.9984;
        //
        // Subsystem 1 -- Matrix B
        //
        subSys->L[0] = 0.0193;
        subSys->L[1] = 4.3142e-04;
        //
        // Subsystem 1 -- Matrix Q
        //
        subSys->Q[0][0] = 0;
        subSys->Q[0][1] = 0;
        subSys->Q[1][0] = 0;
        subSys->Q[1][1] = 1/Ro;
    }
}
