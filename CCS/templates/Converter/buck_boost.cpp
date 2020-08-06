#include <src/Converter/buck_boost.h>

using namespace SwitchedSystem;
using namespace BaseConverter;
using namespace Controller;

extern ControlStrategy controlStrategy;

namespace ConverterBuckBoost
{
    static System system;
    static System discreteSystem;
    static Cycle limitCycle;


    System* BuckBoost::GetSys()
    {
        DefineSystem();

        system.N = 2;

        return &(system);
    }


    System* BuckBoost::GetDiscreteSys()
    {
        DefineDiscreteSystem();

        discreteSystem.N = 2;

        return &(discreteSystem);
    }


    Cycle* BuckBoost::GetLimitCycle()
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


    int BuckBoost::SubSystem2SwitchState(int SubSystem)
    {
        int switchState;

        switch(SubSystem)
        {
        case 0:
            switchState = 2;
            break;
        case 1:
            switchState = 1;
            break;
        default:
            switchState = -1;
            break;
        }

        return switchState;
    }


{{DefineSystem}}


{{DefineDiscreteSystem}}


{{DefineLimitCycleCost}}


{{DefineLimitCycleH2}}


{{DefineLimitCycleHinf}}
}
