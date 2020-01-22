/*
 * cl_MSI_Design_Variable_Interface.hpp
 *
 *  Created on: Jan 10, 2020
 *      Author: schmidt
 */
#ifndef SRC_FEM_CL_MSI_DESIGN_VARIABLE_INTERFACE_HPP_
#define SRC_FEM_CL_MSI_DESIGN_VARIABLE_INTERFACE_HPP_

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Map.hpp"

#include "cl_GEN_Dv_Enums.hpp"

#include "cl_MSI_Dof_Type_Enums.hpp"

namespace moris
{
class Dist_Vector;

namespace mdl
{
    class Model;
}

    namespace MSI
    {
        class Design_Variable_Interface
        {
        private:
//                moris::MSI::Model_Solver_Interface * mMSI = nullptr;
//                moris::MSI::Dof_Manager            * mDofMgn = nullptr;

                Matrix< DDRMat>  mTime;

                mdl::Model * mModel = nullptr;

        public:
            Design_Variable_Interface( )
            {

            };

//------------------------------------------------------------------------------

            ~Design_Variable_Interface() {};

//------------------------------------------------------------------------------
            /**
             * @brief Set model pointeer
             *
             * @param[in] aModel  Model pointer
             *
             */
            void set_model( mdl::Model * aModel )
            {
                mModel = aModel;
            }

//------------------------------------------------------------------------------
            /**
             * @brief Set time
             *
             * @param[in] aTime  Time
             *
             */
            void set_time( const Matrix< DDRMat> & aTime )
            {
                mTime = aTime;
            };

//------------------------------------------------------------------------------
            /**
             * @brief Function providing pdv values for requested vertex indices and Dv types
             *
             * @param[in] aIntegrationMeshSetIndex  Integration Mesh index
             * @param[in] aDvTypes                  List of Dv types
             *
             */
            virtual void get_dv_types_for_set( const moris::moris_index    aIntegrationMeshSetIndex,
                                                     Cell< enum GEN_DV > & aDvTypes ) = 0;

//------------------------------------------------------------------------------
            /**
             * @brief Function providing pdv values for requested vertex indices and Dv types
             *
             * @param[in] aVertexIndex   List of vertex indices
             * @param[in] aDvType        List of Dv types
             * @param[in] aDvValues      List of Dv values
             *
             */
            virtual void get_pdv_value( const moris::Cell< moris::moris_index > & aNodeIndices,
                                        const Cell< enum GEN_DV >               & aDvTypes,
                                              Cell< moris::Matrix< DDRMat > >   & aDvValues,
                                              Cell< moris::Matrix< DDSMat > >   & aIsActiveDv ) = 0;

//------------------------------------------------------------------------------
            /**
             * @brief Retunr local to global dv type map
             *
             */
            virtual moris::Matrix< DDSMat > get_my_local_global_map() = 0;

//------------------------------------------------------------------------------
            /**
             * @brief Retunr local to global dv type map
             *
             * @param[in] aVertexIndex   List of vertex indices
             * @param[in] aDvType        List of Dv types
             * @param[in] aDvIds         List of Dv Ids
             *
             */
            virtual void get_dv_ids_for_type_and_ind( const moris::Cell< moris::moris_index > & aNodeIndices,
                                                      const Cell< enum GEN_DV >               & aDvTypes,
                                                            Cell< moris::Matrix< IdMat > >    & aDvIds ) = 0;

//------------------------------------------------------------------------------
        };
    }
}

#endif /* SRC_FEM_CL_MSI_SOLVER_INTERFACE_HPP_ */
