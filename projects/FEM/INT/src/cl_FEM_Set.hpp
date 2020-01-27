/*
 * cl_FEM_Set.hpp
 *
 *  Created on: Mar 10, 2019
 *      Author: schmidt/noel
 */

#ifndef SRC_FEM_CL_FEM_SET_HPP_
#define SRC_FEM_CL_FEM_SET_HPP_

#include "assert.h"
#include "cl_Communication_Tools.hpp"

#include "cl_MTK_Enums.hpp"            //MTK
#include "cl_MSI_Equation_Set.hpp"     //FEM/MSI/src
#include "cl_FEM_Enums.hpp"            //FEM/INT/src
#include "cl_FEM_Node_Base.hpp"        //FEM/INT/src
#include "cl_FEM_Property.hpp"                       //FEM/INT/src
#include "cl_FEM_Constitutive_Model.hpp"             //FEM/INT/src
#include "cl_FEM_Stabilization_Parameter.hpp"                     //FEM/INT/src
#include "cl_FEM_Set_User_Info.hpp" //FEM/INT/src
#include "cl_FEM_IQI.hpp"

#include "cl_MTK_Cell_Cluster.hpp"
#include "cl_MTK_Side_Cluster.hpp"
#include "cl_MTK_Double_Side_Cluster.hpp"

#include "cl_VIS_Output_Enums.hpp"

namespace moris
{
namespace mtk
{
   class Cell;
   class Set;
}
namespace MSI
{
   class Model_Solver_Interface;
}
    namespace fem
    {
    class IWG;
    class IQI;
    class Field_Interpolator;
    class Geometry_Interpolator;
    class Field_Interpolator_Manager;

//------------------------------------------------------------------------------
    /**
     * @brief element block class that communicates with the mesh interface
     */
    class Set : public MSI::Equation_Set
    {
    private:

        // pointer to the corresponding mesh set
        moris::mtk::Set * mMeshSet = nullptr;

        // list of mesh cluster pointers
        moris::Cell< mtk::Cluster const* > mMeshClusterList;

        // space dimension
        uint mSpaceDim;

        // interpolation mesh geometry type
        mtk::Geometry_Type mIPGeometryType = mtk::Geometry_Type::UNDEFINED;

        // integration mesh geometry type
        mtk::Geometry_Type mIGGeometryType = mtk::Geometry_Type::UNDEFINED;

        // space interpolation order for IP cells
        mtk::Interpolation_Order mIPSpaceInterpolationOrder = mtk::Interpolation_Order::UNDEFINED;

        // space interpolation order for IG cells
        mtk::Interpolation_Order mIGSpaceInterpolationOrder = mtk::Interpolation_Order::UNDEFINED;

        // time interpolation order for IP cells
        mtk::Interpolation_Order mIPTimeInterpolationOrder = mtk::Interpolation_Order::UNDEFINED;

        // space interpolation order for IG cells
        mtk::Interpolation_Order mIGTimeInterpolationOrder = mtk::Interpolation_Order::UNDEFINED;

        // list of fem node pointers for IP nodes
        moris::Cell< Node_Base* > mNodes;

        // field interpolator manager pointers
        Field_Interpolator_Manager * mMasterFIManager = nullptr;
        Field_Interpolator_Manager * mSlaveFIManager  = nullptr;

        // cell of pointers to IWG objects
        moris::Cell< std::shared_ptr< IWG > > mIWGs;
        moris::Cell< std::shared_ptr< IWG > > mRequestedIWGs;

        // cell of pointer to IQI objects
        moris::Cell< std::shared_ptr< IQI > > mIQIs;

        // enum for element type
        enum fem::Element_Type mElementType;

        // integration points
        Matrix< DDRMat > mIntegPoints;

        // integration weights
        Matrix< DDRMat > mIntegWeights;

        // map for the dof type
        Matrix< DDSMat > mDofTypeMap;
        Matrix< DDSMat > mDvTypeMap;

        // map visualization cell id to position in vector
        moris::Cell< moris::Matrix< DDSMat > > mCellAssemblyMap;
        moris::Cell< uint >                    mMtkIgCellOnSet;

        bool mIsTrivialMaster = false;
        bool mIsTrivialSlave  = false;
        bool mIsResidual = false;

        friend class MSI::Equation_Object;
        friend class Cluster;
        friend class Element_Bulk;
        friend class Element_Sideset;
        friend class Element_Double_Sideset;
        friend class Element_Time_Sideset;
        friend class Element;
        friend class Field_Interpolator_Manager;
        friend class Interpolation_Element;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------
        /**
         * constructor
         * @param[ in ] aMeshSet a set from the mesh
         * @param[ in ] aSetInfo user defined info for set
         * @param[ in ] aIPNodes cell of node pointers
         */
        Set( moris::mtk::Set           * aMeshSet,
             fem::Set_User_Info        & aSetInfo,
             moris::Cell< Node_Base* > & aIPNodes );

        /**
         * trivial constructor
         */
        Set()
        {
            mIsEmptySet = true;    //FIXME this flag is a hack. find better solution
        };

//------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~Set();

//------------------------------------------------------------------------------
        /**
         * delete the pointers built on the set
         */
        void delete_pointers();

//------------------------------------------------------------------------------
        /**
         * create and set the field interpolators
         * param[ in ] aModelSolverInterface model solver interface pointer
         */
        void finalize( MSI::Model_Solver_Interface * aModelSolverInterface );

//------------------------------------------------------------------------------
        /**
         * initialize the set
         * @param[ in ] aIsResidual bool true if ???
         */
        void initialize_set( const bool aIsResidual );

//------------------------------------------------------------------------------
        /**
         * free the memory on the set
         */
        void free_memory();

//------------------------------------------------------------------------------
        /**
         * set visualization mesh set
         * @param[ in ] aVisMeshSet a mesh set pointer for visualization
         */
        void set_visualization_set( const uint              aMeshIndex,
                                          moris::mtk::Set * aVisMeshSet,
                                    const bool              aOnlyPrimayCells);

//------------------------------------------------------------------------------
        /**
         * get the element on the set
         * param[ out ] aElementType element type for the set
         */
        enum fem::Element_Type get_element_type() const
        {
            return mElementType;
        }

//------------------------------------------------------------------------------
        /**
         * get the clusters on the set
         * param[ out ] aClusters cell of mesh cluster pointers
         */
        moris::Cell< mtk::Cluster const* > const & get_clusters_on_set() const
        {
            return mMeshClusterList;
        }

//------------------------------------------------------------------------------
        /**
         * create a unique dof type list for the solver
         * Cell< MSI::Dof_Type >, no group of dof type
         * one for both master and slave
         */
        void create_unique_dof_type_list();

//------------------------------------------------------------------------------
        /**
         * create a unique group of dof type list for the set
         * Cell< Cell< MSI::Dof_Type > >, list of groups of dof type
         * one for the master, one for the slave
         */
        void create_dof_type_list();

//------------------------------------------------------------------------------
        /**
         * get dof type list
         * @param[ in ] aIsMaster enum for master or slave
         */
        moris::Cell< moris::Cell< MSI::Dof_Type > > & get_dof_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
        {
            switch ( aIsMaster )
            {
                case ( mtk::Master_Slave::MASTER ):
                {
                    return mMasterDofTypes;
                }
                case( mtk::Master_Slave::SLAVE ):
                {
                    return mSlaveDofTypes;
                }
                default:
                {
                    MORIS_ERROR( false, "Set::get_dof_type_list - can only be MASTER or SLAVE");
                    return mMasterDofTypes;
                }
            }
        }

//------------------------------------------------------------------------------
        /**
         * get dv type list
         * @param[ in ] aIsMaster enum for master or slave
         */
        moris::Cell< moris::Cell< MSI::Dv_Type > > & get_dv_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
        {
            switch ( aIsMaster )
            {
                case ( mtk::Master_Slave::MASTER ):
                {
                    return mMasterDvTypes;
                }
                case( mtk::Master_Slave::SLAVE ):
                {
                    return mSlaveDvTypes;
                }
                default:
                {
                    MORIS_ERROR( false, "Set::get_dv_type_list - can only be MASTER or SLAVE.");
                    return mMasterDvTypes;
                }
            }
        }

//------------------------------------------------------------------------------
        /**
         * create a map of the dof type for the set
         * one for the master, one for the slave
         */
        void create_dof_type_map();

//------------------------------------------------------------------------------
        /**
         * get dof type map
         * @param[ in ] aIsMaster enum for master or slave
         */
        Matrix< DDSMat > & get_dof_type_map( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
        {
            switch ( aIsMaster )
            {
                case ( mtk::Master_Slave::MASTER ):
                {
                    return mMasterDofTypeMap;
                }
                case( mtk::Master_Slave::SLAVE ):
                {
                    return mSlaveDofTypeMap;
                }
                default:
                {
                    MORIS_ERROR(false, "Set::get_dof_type_map - can only be MASTER or SLAVE");
                    return mMasterDofTypeMap;
                }
            }
        }

//------------------------------------------------------------------------------
        /**
         * get residual dof assembly map
         */
        moris::Cell< moris::Matrix< DDSMat > > & get_res_dof_assembly_map()
        {
            return mResDofAssemblyMap;
        }

//------------------------------------------------------------------------------
        /**
         * get jacobian dof assembly map
         */
        moris::Cell< moris::Matrix< DDSMat > > & get_jac_dof_assembly_map()
        {
            return mJacDofAssemblyMap;
        }

//------------------------------------------------------------------------------
        /**
         * get dv assembly map
         */
        moris::Cell< moris::Matrix< DDSMat > > & get_dv_assembly_map()
        {
            return mDvAssemblyMap;
        }

//------------------------------------------------------------------------------
         /**
          * create field interpolator managers for the set
          * @param[ in ] aModelSolverInterface model solver interface
          * ( only used to set the time levels )
          */
         void create_field_interpolator_managers( MSI::Model_Solver_Interface * aModelSolverInterface );

//------------------------------------------------------------------------------
        /**
         * get number of field interpolators
         * @param[ in ] aIsMaster enum for master or slave
         */
        uint get_number_of_field_interpolators( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
        {
            switch ( aIsMaster )
            {
                case ( mtk::Master_Slave::MASTER ):
                {
                    return mMasterDofTypes.size();
                }
                case( mtk::Master_Slave::SLAVE ):
                {
                    return mSlaveDofTypes.size();
                }
                default:
                {
                    MORIS_ERROR(false, "Set::get_number_of_field_interpolators - can only be MASTER or SLAVE.");
                    return mMasterDofTypes.size();
                }
            }
        }

//------------------------------------------------------------------------------
        /**
         * get number of IWGs
         */
        uint get_number_of_IWGs()
        {
            return mIWGs.size();
        }

//------------------------------------------------------------------------------
        /**
         * get number of requested IWGs
         */
        uint get_number_of_requested_IWGs()
        {
            return mRequestedIWGs.size();
        }

//------------------------------------------------------------------------------
        /**
         * get IWGs
         * param[ out ] aIWGs cell of IWG pointers
         */
        moris::Cell< std::shared_ptr< IWG > > & get_IWGs()
        {
            return mIWGs;
        }

//------------------------------------------------------------------------------
        /**
         * get IWGs
         * param[ out ] aIWGs cell of IWG pointers
         */
        moris::Cell< std::shared_ptr< IWG > > & get_requested_IWGs()
        {
            return mRequestedIWGs;
        }

//------------------------------------------------------------------------------
        /**
         * get IQI from type
         * param[ out ] aIWGs cell of IWG pointers
         */
        std::shared_ptr< IQI > get_requested_IQI( enum vis::Output_Type aOutputType )
        {
            // FIXME use a map

            // init the IQI pointer for return
            std::shared_ptr< IQI > tIQI = nullptr;

            // select the IQI based on type
            for ( uint iIQI = 0; iIQI < mIQIs.size(); iIQI++ )
            {
                // if the output type is the same as the IQI
                if ( aOutputType == mIQIs( iIQI )->get_IQI_type() )
                {
                    // fill the return pointer to the IQI
                    tIQI = mIQIs( iIQI );
                }
            }

            // return the selected IQI
            return tIQI;
        }

//------------------------------------------------------------------------------
        /**
         * set the field interpolator managers for the IWGs
         */
        void set_IWG_field_interpolator_managers();

//------------------------------------------------------------------------------
        /**
         * set the field interpolator managers for the IQIs
         */
        void set_IQI_field_interpolator_managers();

//------------------------------------------------------------------------------
        /*!
         * Sets the cluster in the stabilization parameter associated with
         * this set
         */
        void
        set_cluster_in_stabilization_params(fem::Cluster * aCluster);
//------------------------------------------------------------------------------
        /**
         * create the dof assembly map for the residual/rows
         */
        void create_residual_dof_assembly_map();

//------------------------------------------------------------------------------
        /**
         * create the dof assembly map for the jacobian/cols
         */
        void create_dof_assembly_map( const bool aIsResidual );

//------------------------------------------------------------------------------
        /**
         * create the dof assembly map for the jacobian/cols
         */
        void create_jacobian_dof_assembly_map();

//------------------------------------------------------------------------------
        /**
         * create the dof assembly map
         * for the off-diagonal requested jacobian/cols
         * for R = R_0 - A_{01} x_{1}
         */
        void create_staggered_jacobian_dof_assembly_map();

//------------------------------------------------------------------------------
        /**
         * create the dv assembly map
         */
        void create_dv_assembly_map();

//------------------------------------------------------------------------------
        /**
         * create a list of IWGs requested by the solver
         */
        void create_requested_IWG_list();

//------------------------------------------------------------------------------
        /**
         * create a dof type list for the list of IWGs requested by the solver
         */
        void build_requested_IWG_dof_type_list( const bool aItResidual );

//------------------------------------------------------------------------------
        /**
         * get the IP geometry type
         * @param[ out ] aGeometryType
         */
        mtk::Geometry_Type get_IP_geometry_type()
        {
            return mIPGeometryType;
        }

//------------------------------------------------------------------------------
        /**
         * get the IG geometry type
         * @param[ out ] aGeometryType
         */
        mtk::Geometry_Type get_IG_geometry_type()
        {
            return mIGGeometryType;
        }

//------------------------------------------------------------------------------
        /**
         * FIXME we should not be able to ask this question
         * as the interpolator order is not unique on the set
         * get the IG space interpolation order
         */
        mtk::Interpolation_Order get_IG_space_interpolation_order()
        {
            return mIGSpaceInterpolationOrder;
        }

//------------------------------------------------------------------------------
        /**
         * get the number of integration points
         * @param[ out ] aNumIntegPoints number of integration points
         */
        uint get_number_of_integration_points()
        {
            return mIntegWeights.numel();
        }

//------------------------------------------------------------------------------
        /**
         * get the integration points
         * @param[ out ] aIntegPoints integration points
         */
        const Matrix< DDRMat > & get_integration_points()
        {
            return mIntegPoints;
        }

//------------------------------------------------------------------------------
        /**
         * get the integration weights
         * @param[ out ] aIntegWeights integration weights
         */
        const Matrix< DDRMat > & get_integration_weights()
        {
            return mIntegWeights;
        }

//------------------------------------------------------------------------------
        /**
         * get the field interpolator manager
         * @param[ in ]  aIsMaster an enum for master or slave
         * @param[ out ] mFIManger a field interpolator manager pointer
         */
        Field_Interpolator_Manager * get_field_interpolator_manager( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
        {
            switch ( aIsMaster )
            {
                case ( mtk::Master_Slave::MASTER ) :
                    return mMasterFIManager;

                case ( mtk::Master_Slave::SLAVE ) :
                    return mSlaveFIManager;

                default :
                {
                    MORIS_ERROR( false, "Set::get_field_interpolator_manager - can only be master or slave.");
                    return mMasterFIManager;
                }
            }
        };

//------------------------------------------------------------------------------
        /**
         * auto detect full integration scheme
         * @param[ in ] aGeometryType a geometry type
         */
        fem::Integration_Order get_auto_integration_order( const mtk::Geometry_Type aGeometryType );

//------------------------------------------------------------------------------
        /**
         * auto detect space interpolation scheme
         * @param[ in ] aNumVertices  number of vertices for the considered geometry type
         * @param[ in ] aGeometryType geometry type
         * FIXME: works for Lagrange only
         */
        mtk::Interpolation_Order get_auto_interpolation_order( const moris::uint        aNumVertices,
                                                               const mtk::Geometry_Type aGeometryType );

//------------------------------------------------------------------------------
        /**
          * auto detect time interpolation scheme
          * @param[ in ] aNumVertices number of vertices for a line
          */
        fem::Interpolation_Type get_auto_time_interpolation_type( const moris::uint aNumVertices );

//------------------------------------------------------------------------------
        /**
         * set size and reset values for Jacobian
         */
        void initialize_mJacobian();

//------------------------------------------------------------------------------
        /**
         * set size and reset values for residual
         */
        void initialize_mResidual();

//------------------------------------------------------------------------------
        /*
         * get the set index for the specified dof type
         *@param[ in ] aDofType a dof type enum
         */
        moris::sint get_dof_index_for_type_1( enum MSI::Dof_Type aDofType );

//------------------------------------------------------------------------------
        /**
         * get the set index for the specified dv type
         *@param[ in ] aDvType a dv type enum
         * FIXME should this info be stored in the set?
         */
        moris::sint get_dv_index_for_type_1( enum MSI::Dv_Type aDvType );

//------------------------------------------------------------------------------
        /**
         * get number of dof types on the set
         */
        moris::uint get_num_dof_types();

//------------------------------------------------------------------------------
        /**
         * get number of dv types on the set
         * FIXME where should this be stored???
         */
        moris::uint get_num_dv_types()
        {
            return 1;
        }

//------------------------------------------------------------------------------

        moris::Cell < enum MSI::Dof_Type > get_requested_dof_types();

//------------------------------------------------------------------------------

        moris::Cell< moris::Cell< enum MSI::Dof_Type > > get_secundary_dof_types();

//------------------------------------------------------------------------------

        void create_dof_type_map_unique()
        {
            // Create temporary dof type list
            moris::Cell< enum MSI::Dof_Type > tDofType = get_unique_dof_type_list();

            //Get number of unique adofs of this equation object
            moris::uint tNumUniqueDofTypes = tDofType.size();

            // Get maximal dof type enum number
            moris::sint tMaxDofTypeEnumNumber = 0;

            // Loop over all pdof types to get the highest enum index
            for ( moris::uint Ii = 0; Ii < tNumUniqueDofTypes; Ii++ )
            {
                tMaxDofTypeEnumNumber = std::max( tMaxDofTypeEnumNumber, static_cast< int >( tDofType( Ii ) ) );
            }

            // +1 because c++ is 0 based
            tMaxDofTypeEnumNumber++;

            // Set size of mapping matrix
            mDofTypeMap       .set_size( tMaxDofTypeEnumNumber, 1, -1 );

            // Loop over all pdof types to create the mapping matrix
            for ( moris::uint Ii = 0; Ii < tNumUniqueDofTypes; Ii++ )
            {
                mDofTypeMap( static_cast< int >( tDofType( Ii ) ), 0 ) = Ii;
            }
        }

//------------------------------------------------------------------------------
        /**
         * compute a quantity of interest
         * @param[ in ] aElementFieldValues
         * @param[ in ] aNodelFieldValues
         * @param[ in ] aGlobalScalar
         * @param[ in ] aOutputType
         * @param[ in ] aFieldType
         */
        void compute_quantity_of_interest( const uint            aMeshIndex,
                                           Matrix< DDRMat >      * aElementFieldValues,
                                           Matrix< DDRMat >      * aNodalFieldValues,
                                           moris::real           * aGlobalScalar,
                                           enum vis::Output_Type   aOutputType,
                                           enum vis::Field_Type    aFieldType );

//------------------------------------------------------------------------------
        /**
         * determine set type from mtk set type
         */
        void determine_set_type();

    };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_SET_HPP_ */
