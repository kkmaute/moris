//
// Created by christopherson on 3/4/20.
//

#include "fn_OPT_create_interface.hpp"
#include "cl_OPT_Interface_User_Defined.hpp"
#include "cl_OPT_Interface_Manager.hpp"

namespace moris
{
    namespace opt
    {
        //--------------------------------------------------------------------------------------------------------------

        std::shared_ptr<Interface> create_interface(Cell<ParameterList> aParameterLists)
        {
            // Get number of interfaces
            uint tNumInterfaces = aParameterLists.size() - 1;

            // Single interface without manager
            if (tNumInterfaces == 0)
            {
                return create_interface(aParameterLists(0));
            }

            // Multiple interfaces, create interface manager
            else
            {
                Cell<std::shared_ptr<Interface>> tInterfaces(tNumInterfaces);
                for (uint tInterfaceIndex = 0; tInterfaceIndex < tNumInterfaces; tInterfaceIndex++)
                {
                    tInterfaces(tInterfaceIndex) = create_interface(aParameterLists(tInterfaceIndex + 1));
                }
                return std::make_shared<Interface_Manager>(aParameterLists(0), tInterfaces);
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        std::shared_ptr<Interface> create_interface(ParameterList aParameterList)
        {
            std::string tInterfaceType = aParameterList.get<std::string>("type");
            if (!tInterfaceType.compare("user_defined"))
            {
                return std::make_shared<Interface_User_Defined>(aParameterList);
            }
            else
            {
                MORIS_ERROR(false, tInterfaceType.append(" is not recognized as a valid Interface type in fn_OPT_create_interface.").c_str());
                return nullptr;
            }
        }

        //--------------------------------------------------------------------------------------------------------------
    }
}
