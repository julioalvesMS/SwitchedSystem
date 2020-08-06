#include <src/Converter/buck.h>

using namespace SwitchedSystem;
using namespace BaseConverter;
using namespace Controller;

extern ControlStrategy controlStrategy;

namespace ConverterBuck
{
    static System system;
    static System discreteSystem;
    static Cycle limitCycle;


    System* Buck::GetSys()
    {
        DefineSystem();

        system.N = 2;

        return &(system);
    }


    System* Buck::GetDiscreteSys()
    {
        DefineDiscreteSystem();

        discreteSystem.N = 2;

        return &(discreteSystem);
    }


    Cycle* Buck::GetLimitCycle()
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


{{getClassicVoltageController}}


{{getClassicVoltageCurrentController}}


{{getStateFeedbackH2Controller}}


{{getReferenceController}}


{{getCurrentCorrectionController}}


    int Buck::SubSystem2SwitchState(int SubSystem)
    {
        return SubSystem;
    }


{{DefineSystem}}


{{DefineDiscreteSystem}}


{{DefineLimitCycleCost}}


{{DefineLimitCycleH2}}


{{DefineLimitCycleHinf}}
}
