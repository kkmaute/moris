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

//#include "cl_Matrix_Vector_Factory.hpp"
//#include "cl_SOL_Dist_Map.hpp"

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

//                Dist_Map * mVectorMap = nullptr;
//
//                //! Full Vector
//                Dist_Vector * mVector = nullptr;

        public:

//------------------------------------------------------------------------------
            /**
             * trivial constructor
             */
            Design_Variable_Interface( ){};

//------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~Design_Variable_Interface(){};

//------------------------------------------------------------------------------
            /**
             * set model pointer
             * @param[ in ] aModel Model pointer
             */
            void set_model( mdl::Model * aModel )
            {
                mModel = aModel;
            }

//------------------------------------------------------------------------------
            /**
             * set time
             * @param[ in ] aTime Time
             */
            void set_time( const Matrix< DDRMat> & aTime )
            {
                mTime = aTime;
            };

//------------------------------------------------------------------------------
            /**
             * get unique dv types for set
             * @param[ in ] aIntegrationMeshSetIndex
             * @param[ in ] aDvTypes
             */
            virtual void get_unique_dv_types_for_set( const moris::moris_index    aIntegrationMeshSetIndex,
                                                            Cell< enum GEN_DV > & aDvTypes ) = 0;

//------------------------------------------------------------------------------
            /**
             * get requested dv types for sensitivity analysis
             * @param[ in ] aDvTypes list of dv types to fill
             */
            virtual void get_requested_dv_types( Cell< enum GEN_DV > & aDvTypes )
            {
                MORIS_ERROR( false, "get_requested_dv_types - not implemented in msi base class.");
            }

//------------------------------------------------------------------------------
            /**
             * get pdv values for requested vertex indices and dv types
             * @param[ in ] aIntegrationMeshSetIndex  integration Mesh index
             * @param[ in ] aDvTypes                  list of group of dv types
             */
            virtual void get_dv_types_for_set( const moris::moris_index            aIntegrationMeshSetIndex,
                                                     Cell< Cell< enum GEN_DV > > & aDvTypes ) = 0;

//------------------------------------------------------------------------------
            /**
             * get pdv values for requested vertex indices and dv types
             * @param[ in ]     aNodeIndices list of vertex indices
             * @param[ in ]     aDvTypes     list of dv types
             * @param[ in/out ] aDvValues    list of dv values
             * @param[ in/out ] aIsActiveDv  list of active whether or not dv is active
             */
            virtual void get_pdv_value( const Matrix< IndexMat >                & aNodeIndices,
                                        const moris::Cell< enum GEN_DV >        & aDvTypes,
                                        moris::Cell< moris::Matrix< DDRMat > >  & aDvValues,
                                        moris::Cell< moris::Matrix< DDSMat > >  & aIsActiveDv ) = 0;
//------------------------------------------------------------------------------
            /**
             * get pdv values for requested vertex indices and dv types
             * @param[ in ]     aNodeIndices list of vertex indices
             * @param[ in ]     aDvTypes     list of dv types
             * @param[ in/out ] aDvValues    list of dv values
             */
            virtual void get_pdv_value( const Matrix< IndexMat >                & aNodeIndices,     // TODO: temporary, once the interface is finalized, need FEM to play nice with this
                                        const moris::Cell< enum GEN_DV >        & aDvTypes,
                                        moris::Cell< moris::Matrix< DDRMat > >  & aDvValues ) = 0;

//------------------------------------------------------------------------------
            /**
             * reshape pdv values
             * i.e. reshape a cell of matrix to a matrix
             * @param[ in ] aPdvValues         a cell of matrices with pdv values
             * @param[ in ] aReshapedPdvValues a matrix of pdv values
             */
            virtual void reshape_pdv_values( const moris::Cell< moris::Matrix< DDRMat > > & aPdvValues,
                                                   moris::Matrix< DDRMat >                & aReshapedPdvValues )
            {
                MORIS_ERROR( false, "Design_Variable_Interface::reshape_pdv_values - not implemented for base class." );
            }

//------------------------------------------------------------------------------
            /**
             * return local to global dv type map
             */
            virtual moris::Matrix< DDSMat > get_my_local_global_map() = 0;

//------------------------------------------------------------------------------
            /**
             * return local to global dv type map
             * @param[ in ] aVertexIndex   List of vertex indices
             * @param[ in ] aDvType        List of Dv types
             * @param[ in ] aDvIds         List of Dv Ids
             */
            virtual void get_dv_ids_for_type_and_ind( const moris::Cell< moris::moris_index > & aNodeIndices,
                                                      const Cell< enum GEN_DV >               & aDvTypes,
                                                            Cell< moris::Matrix< IdMat > >    & aDvIds ) = 0;
//------------------------------------------------------------------------------
            /**
             * get QI assembly map
             */
            virtual moris::Cell< moris::Cell< moris_index > > & get_QI_assembly_map() = 0;


        };
    }
}

#endif /* SRC_FEM_CL_MSI_SOLVER_INTERFACE_HPP_ */
