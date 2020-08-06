#include <src/Converter/buck_boost_3.h>

using namespace SwitchedSystem;
using namespace BaseConverter;
using namespace Controller;

extern ControlStrategy controlStrategy;

namespace ConverterBuckBoost3
{
    static System system;
    static System discreteSystem;
    static Cycle limitCycle;


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


    Cycle* BuckBoost3::GetLimitCycle()
    {
        switch(controlStrategy)
        {
        case CS_LIMIT_CYCLE_COST:
            DefineLimitCycleCost();
            break;
        case CS_LIMIT_CYCLE_H2:
            DefineLimitCycleH2();
            break;
        case CS_LIMIT_CYCLE_Hinf:
            DefineLimitCycleHinf();
        default:
            break;
        }

        return &(limitCycle);
    }


{{getP}}


{{getReferenceController}}


{{getCurrentCorrectionController}}


    int BuckBoost3::SubSystem2SwitchState(int SubSystem)
    {
        return SubSystem;
    }


{{DefineSystem}}


{{DefineDiscreteSystem}}


{{DefineLimitCycleCost}}


{{DefineLimitCycleH2}}


{{DefineLimitCycleHinf}}
}
