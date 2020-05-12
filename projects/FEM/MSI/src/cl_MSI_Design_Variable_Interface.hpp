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

#include "cl_GEN_Pdv_Enums.hpp"
#include "cl_FEM_Enums.hpp"

#include "cl_MSI_Dof_Type_Enums.hpp"

//#include "cl_Matrix_Vector_Factory.hpp"
//#include "cl_SOL_Dist_Map.hpp"

namespace moris
{
    namespace sol
    {
        class Dist_Vector;
    }

namespace mdl
{
    class Model;
}

    namespace MSI
    {
        class Design_Variable_Interface
        {
        private:

                Matrix< DDRMat>  mTime;

        protected:
                mdl::Model * mModel = nullptr;

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
            virtual void get_ip_unique_dv_types_for_set( const moris::moris_index    aIntegrationMeshSetIndex,
                                                               Cell< enum PDV_Type > & aDvTypes ) = 0;

            virtual void get_ig_unique_dv_types_for_set( const moris::moris_index    aIntegrationMeshSetIndex,
                                                               Cell< enum PDV_Type > & aDvTypes ) = 0;

//------------------------------------------------------------------------------
            /**
             * get pdv values for requested vertex indices and dv types
             * @param[ in ] aIntegrationMeshSetIndex  integration Mesh index
             * @param[ in ] aDvTypes                  list of group of dv types
             */
            virtual void get_ip_dv_types_for_set( const moris::moris_index            aIntegrationMeshSetIndex, //**********
                                                        Cell< Cell< enum PDV_Type > > & aDvTypes ) = 0;

            virtual void get_ig_dv_types_for_set( const moris::moris_index            aIntegrationMeshSetIndex,
                                                        Cell< Cell< enum PDV_Type > > & aDvTypes ) = 0;

//------------------------------------------------------------------------------
            /**
             * get pdv values for requested vertex indices and dv types
             * @param[ in ]     aNodeIndices list of vertex indices
             * @param[ in ]     aDvTypes     list of dv types
             * @param[ in/out ] aDvValues    list of dv values
             * @param[ in/out ] aIsActiveDv  list of active whether or not dv is active
             */
            virtual void get_ip_pdv_value( const Matrix< IndexMat >                & aNodeIndices,
                                           const moris::Cell< enum PDV_Type >        & aDvTypes,
                                           moris::Cell< moris::Matrix< DDRMat > >  & aDvValues,
                                           moris::Cell< moris::Matrix< DDSMat > >  & aIsActiveDv ) = 0;

            virtual void get_ig_pdv_value( const Matrix< IndexMat >                & aNodeIndices,
                                           const moris::Cell< enum PDV_Type >        & aDvTypes,
                                           moris::Cell< moris::Matrix< DDRMat > >  & aDvValues,
                                           moris::Cell< moris::Matrix< DDSMat > >  & aIsActiveDv ) = 0;
//------------------------------------------------------------------------------
            /**
             * get pdv values for requested vertex indices and dv types
             * @param[ in ]     aNodeIndices list of vertex indices
             * @param[ in ]     aDvTypes     list of dv types
             * @param[ in/out ] aDvValues    list of dv values
             */
            virtual void get_ip_pdv_value( const Matrix< IndexMat >                & aNodeIndices,     // TODO: does this need to be overloaded?
                                           const moris::Cell< enum PDV_Type >        & aDvTypes,
                                           moris::Cell< moris::Matrix< DDRMat > >  & aDvValues ) = 0;

            virtual void get_ig_pdv_value( const Matrix< IndexMat >                & aNodeIndices,     // TODO: does this need to be overloaded?
                                           const moris::Cell< enum PDV_Type >        & aDvTypes,
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
             * return local to global dv map
             */
            virtual moris::Matrix< DDSMat > get_my_local_global_map() = 0;

//------------------------------------------------------------------------------
            /**
             * return local to global dv type map
             * @param[ in ] aVertexIndex   List of vertex indices
             * @param[ in ] aDvType        List of Dv types
             * @param[ in ] aDvIds         List of Dv Ids
             */
            virtual void get_ip_dv_ids_for_type_and_ind( const Matrix<IndexMat> & aNodeIndices,
                                                         const Cell< enum PDV_Type >               & aDvTypes,
                                                               Cell< moris::Matrix< IdMat > >    & aDvIds ) = 0;

            virtual void get_ig_dv_ids_for_type_and_ind( const Matrix<IndexMat> & aNodeIndices,
                                                         const Cell< enum PDV_Type >               & aDvTypes,
                                                               Cell< moris::Matrix< IdMat > >    & aDvIds ) = 0;
//------------------------------------------------------------------------------
            /**
             * get requested dv types for sensitivity analysis
             * @param[ in ] aDvTypes list of dv types to fill
             */
            virtual void get_ip_requested_dv_types( Cell< enum PDV_Type > & aDvTypes ) = 0;

//------------------------------------------------------------------------------
            /**
             * get requested dv types for sensitivity analysis
             * @param[ in ] aDvTypes list of dv types to fill
             */
            virtual void get_ig_requested_dv_types( Cell< enum PDV_Type > & aDvTypes ) = 0;

//------------------------------------------------------------------------------
            /**
             * set requested IQI type for sensitivity analysis
             * @param[ in ] aRequestedIQIType
             */
            virtual void set_requested_IQI_type( const moris::Cell< moris::Cell< enum fem::IQI_Type > > & aRequestedIQIType )
            {
                MORIS_ERROR( false, "Design_Variable_Interface::set_requested_IQI_type - not implemented for base class." );
            };

//------------------------------------------------------------------------------
            /**
             * set requested IQI type for sensitivity analysis
             * @param[ in ] aRequestedIQIType
             */
            virtual void set_requested_IQIs( const moris::Cell< std::string > & aRequestedIQIs )
            {
                MORIS_ERROR( false, "Design_Variable_Interface::set_requested_IQIs - not implemented for base class." );
            };

        };
    }
}

#endif /* SRC_FEM_CL_MSI_SOLVER_INTERFACE_HPP_ */
