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

#include "cl_Matrix_Vector_Factory.hpp"
#include "cl_Map_Class.hpp"

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

        protected:
                Map_Class * mVectorMap = nullptr;

                //! Full Vector
                Dist_Vector * mVector = nullptr;

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
             * @brief retunr dRdp pointer
             *
             */
            Dist_Vector * get_dRdp()
            {
                return mVector;
            };

            virtual void get_unique_dv_types_for_set( const moris::moris_index    aIntegrationMeshSetIndex,
                                                            Cell< enum GEN_DV > & aDvTypes )
            {
                //FIXME: need to change everywhere so that we have everything using GEN_DV list
                MORIS_ASSERT(false, "Design_Variable_Interface::get_unique_dv_types_for_set() - not implemented on GE side");
            }

            virtual void get_unique_dv_types_for_set( const moris::moris_index          aIntegrationMeshSetIndex,
                                                            Cell< enum MSI::Dv_Type > & aDvTypes )
            {
                //FIXME: need to change everywhere so that we have everything using GEN_DV list
                MORIS_ASSERT(false, "Design_Variable_Interface::get_unique_dv_types_for_set() - not implemented on GE side");
            }

//------------------------------------------------------------------------------
            /**
             * @brief Function providing pdv values for requested vertex indices and Dv types
             *
             * @param[in] aIntegrationMeshSetIndex  Integration Mesh index
             * @param[in] aDvTypes                  List of Dv types
             *
             */
            virtual void get_dv_types_for_set( const moris::moris_index          aIntegrationMeshSetIndex,
                                                     Cell< Cell< enum GEN_DV > > & aDvTypes )
            {
                //FIXME: need to change everywhere so that we have everything using GEN_DV list
                MORIS_ASSERT(false, "Design_Variable_Interface::get_dv_types_for_set() - not implemented on GE side");
            }

            virtual void get_dv_types_for_set( const moris::moris_index          aIntegrationMeshSetIndex,
                                                     Cell< Cell< enum MSI::Dv_Type > > & aDvTypes )
            {
                //FIXME: need to change everywhere so that we have everything using GEN_DV list
                MORIS_ASSERT(false, "Design_Variable_Interface::get_dv_types_for_set() - not implemented on GE side");
            }

//------------------------------------------------------------------------------
            /**
             * @brief Function providing pdv values for requested vertex indices and Dv types
             *
             * @param[in] aNodeIndices    - List of vertex indices
             * @param[in] aDvTypes        - List of Dv types
             * @param[in/out] aDvValues   - List of Dv values
             * @param[in/out] aIsActiveDv - List of active whether or not DV is active
             *
             */
            virtual void get_pdv_value( const Matrix< IndexMat >                & aNodeIndices,
                                        const moris::Cell< enum GEN_DV >        & aDvTypes,
                                        moris::Cell< moris::Matrix< DDRMat > >  & aDvValues,
                                        moris::Cell< moris::Matrix< DDSMat > >  & aIsActiveDv )
            {
                //FIXME: these functions [get_pdv_values()] need to be pure virtual
            }

            virtual void get_pdv_value( const Matrix< IndexMat >                & aNodeIndices,
                                        const moris::Cell< enum MSI::Dv_Type >  & aDvTypes,
                                        moris::Cell< moris::Matrix< DDRMat > >  & aDvValues,
                                        moris::Cell< moris::Matrix< DDSMat > >  & aIsActiveDv )
            {
                //FIXME: need to change everywhere so that we have everything using GEN_DV list
                MORIS_ASSERT(false, "Design_Variable_Interface::get_pdv_value() - not implementd on GE side");
            }

            virtual void get_pdv_value( const Matrix< IndexMat >                & aNodeIndices,     // temporary, once the interface is finalized, need to make FEM play nice with it
                                        const moris::Cell< enum MSI::Dv_Type >  & aDvTypes,
                                        moris::Cell< moris::Matrix< DDRMat > >  & aDvValues )
            {
                MORIS_ASSERT(false, "Design_Variable_Interface::get_pdv_value() - not implementd on GE side");
            }
//------------------------------------------------------------------------------
            /**
             * @brief Retunr local to global dv type map
             *
             */
            virtual moris::Matrix< DDSMat > get_my_local_global_map() = 0;

//------------------------------------------------------------------------------
            /**
             * @brief Return local to global dv type map
             *
             * @param[in] aVertexIndex   List of vertex indices
             * @param[in] aDvType        List of Dv types
             * @param[in] aDvIds         List of Dv Ids
             *
             */
            virtual void get_dv_ids_for_type_and_ind( const moris::Cell< moris::moris_index > & aNodeIndices,
                                                      const Cell< enum GEN_DV >               & aDvTypes,
                                                            Cell< moris::Matrix< IdMat > >    & aDvIds )
            {
                //FIXME: need to change everywhere so that we have everything using GEN_DV list
                MORIS_ASSERT(false, "Design_Variable_Interface::get_dv_ids_for_type_and_ind() - not implemented on GE side");
            }

            virtual void get_dv_ids_for_type_and_ind( const moris::Cell< moris::moris_index > & aNodeIndices,
                                                      const Cell< enum MSI::Dv_Type >         & aDvTypes,
                                                            Cell< moris::Matrix< IdMat > >    & aDvIds )
            {
                //FIXME: need to change everywhere so that we have everything using GEN_DV list
                MORIS_ASSERT(false, "Design_Variable_Interface::get_dv_ids_for_type_and_ind() - not implemented on GE side");
            }

//------------------------------------------------------------------------------
        };
    }
}

#endif /* SRC_FEM_CL_MSI_SOLVER_INTERFACE_HPP_ */
