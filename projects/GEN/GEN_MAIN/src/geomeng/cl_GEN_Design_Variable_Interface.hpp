/*
 * cl_GEN_Design_Variable_Interface.hpp
 *
 *  Created on: Jan 15, 2020
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_SRC_GEOMENG_CL_GEN_DESIGN_VARIABLE_INTERFACE_HPP_
#define PROJECTS_GEN_SRC_GEOMENG_CL_GEN_DESIGN_VARIABLE_INTERFACE_HPP_

#include "cl_GEN_Pdv_Host_Manager.hpp"

#include "cl_MSI_Design_Variable_Interface.hpp"

namespace moris
{
namespace ge
{
    class GEN_Design_Variable_Interface : MSI::Design_Variable_Interface
    {

    private:

        // pdv host manager pointer
        Pdv_Host_Manager* mPdvHostManager;

    public:

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
        /**
         * get requested dv types for sensitivity analysis
         * @param[ in ] aDvTypes list of dv types to fill
         */
        void get_ip_requested_dv_types(Cell<GEN_DV> & aDvTypes)
        {
            // total number of unique DV types
            uint tNumDvTypes = mPdvHostManager->get_ip_pdv_type_list().size();

            // list of all DV types which are changing
            aDvTypes.reserve(tNumDvTypes);

            // total number of unchanging DV types
            uint tNumUnchangingDvTypes = mPdvHostManager->get_ip_unchanging_type_list().size();

            for (uint i=0; i<tNumDvTypes; i++)
            {
                bool tIsChanging;
                for (uint j=0; j<tNumUnchangingDvTypes; j++)
                {
                    if (mPdvHostManager->get_ip_pdv_type_list()(i) == mPdvHostManager->get_ip_unchanging_type_list()(j))
                    {
                        // this type is not changing and therefore is not added to the list
                        tIsChanging = false;
                        break;
                    }
                    else
                    {
                        // this type is changing and so we add it to the list
                        tIsChanging = true;
                    }
                }
                // use flag to make assignment
                if (tIsChanging)
                {
                    aDvTypes.push_back(mPdvHostManager->get_ip_pdv_type_list()(i));
                }
            }

            // shrink list to appropriate size
            aDvTypes.shrink_to_fit();

            // get rid of redundant entries
            unique(aDvTypes);
        }
//------------------------------------------------------------------------------
        /**
         * Returns a cell of all GEN_DV types which are "changing" on the IG mesh
         */
        void get_ig_requested_dv_types(Cell<GEN_DV> & aDvTypes)
        {
            // total number of unique DV types
            uint tNumDvTypes = mPdvHostManager->get_ig_pdv_type_list().size();

            // list of all DV types which are changing
            aDvTypes.reserve(tNumDvTypes);

            // total number of unchanging DV types
            uint tNumUnchangingDvTypes = mPdvHostManager->get_ig_unchanging_type_list().size();

            for (uint i=0; i<tNumDvTypes; i++)
            {
                bool tIsChanging;
                for (uint j=0; j<tNumUnchangingDvTypes; j++)
                {
                    if (mPdvHostManager->get_ig_pdv_type_list()(i) == mPdvHostManager->get_ig_unchanging_type_list()(j))
                    {
                        // this type is not changing and therefore is not added to the list
                        tIsChanging = false;
                        break;
                    }
                    else
                    {
                        // this type is changing and so we add it to the list
                        tIsChanging = true;
                    }
                }
                // use flag to make assignment
                if (tIsChanging)
                {
                    aDvTypes.push_back(mPdvHostManager->get_ig_pdv_type_list()(i));
                }
            }

            // shrink list to appropriate size
            aDvTypes.shrink_to_fit();

            // get rid of redundant entries
            unique(aDvTypes);
        }
    };

}   // end ge namespace
}   // end moris namespace

#endif /* PROJECTS_GEN_SRC_GEOMENG_CL_GEN_DESIGN_VARIABLE_INTERFACE_HPP_ */
