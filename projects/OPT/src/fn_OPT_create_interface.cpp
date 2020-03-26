//
// Created by christopherson on 3/4/20.
//

#include "fn_OPT_create_interface.hpp"
#include "cl_OPT_Interface_User_Defined.hpp"

namespace moris
{
    namespace opt
    {
        std::shared_ptr<Interface> create_interface(ParameterList aParameterList)
        {
            std::string tInterfaceType = aParameterList.get<std::string>("interface");
            if (!tInterfaceType.compare("user_defined"))
            {
                return std::make_shared<Interface_User_Defined>(aParameterList);
            }
            else
            {
                MORIS_ERROR(false, tInterfaceType.append(" is not recognized as a valid Problem in fn_OPT_create_problem.").c_str());
                return nullptr;
            }
        }
    }
}
