#include "fn_OPT_create_problem.hpp"
#include "cl_OPT_Problem_User_Defined.hpp"

namespace moris
{
    namespace opt
    {
        std::shared_ptr<Problem> create_problem(ParameterList aProblemParameterList, std::shared_ptr<Criteria_Interface> aInterface)
        {
            std::string tProblemType = aProblemParameterList.get<std::string>("problem");

            sint tReinitializeInterfaceIter = aProblemParameterList.get<sint>("reinitialize_interface_iter");
            sint tFirstReinitializeInterfaceIter = aProblemParameterList.get<sint>("first_reinitialize_interface_iter");

            aInterface->set_reinitialize_iter( tReinitializeInterfaceIter );

            aInterface->set_first_reinitialize_iter( tFirstReinitializeInterfaceIter );

            if (!tProblemType.compare("user_defined"))
            {
                return std::make_shared<Problem_User_Defined>(aProblemParameterList, aInterface);
            }
            else
            {
                MORIS_ERROR(false, tProblemType.append(" is not recognized as a valid Problem in fn_OPT_create_problem.").c_str());
                return nullptr;
            }
        }
    }
}
