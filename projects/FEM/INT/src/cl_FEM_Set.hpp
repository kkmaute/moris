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
#include "cl_FEM_CM_Factory.hpp"                     //FEM/INT/src
#include "cl_FEM_Set_User_Info.hpp" //FEM/INT/src

#include "cl_MTK_Cell_Cluster.hpp"
#include "cl_MTK_Side_Cluster.hpp"
#include "cl_MTK_Double_Side_Cluster.hpp"

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
    class Field_Interpolator;
    class Geometry_Interpolator;

//------------------------------------------------------------------------------
    /**
     * \brief element block class that communicates with the mesh interface
     */
    class Set : public MSI::Equation_Set
    {
    private:
        moris::mtk::Set * mMeshSet = nullptr;

        //! Mesh cluster
        moris::Cell< mtk::Cluster const* > mMeshClusterList;

        // space dimension
        uint mSpaceDim;

        // interpolation mesh geometry type
        mtk::Geometry_Type mIPGeometryType;

        // integration mesh geometry type
        mtk::Geometry_Type mIGGeometryType;

        // space interpolation order for IP cells
        mtk::Interpolation_Order mIPSpaceInterpolationOrder;

        // space interpolation order for IG cells
        mtk::Interpolation_Order mIGSpaceInterpolationOrder;

        // time interpolation order for IP cells
        mtk::Interpolation_Order mIPTimeInterpolationOrder;

        // space interpolation order for IG cells
        mtk::Interpolation_Order mIGTimeInterpolationOrder;

        // List of fem node pointers
        moris::Cell< Node_Base* >           mNodes;

        // geometry interpolator pointer for the interpolation cells
        Geometry_Interpolator             * mMasterIPGeometryInterpolator = nullptr;
        Geometry_Interpolator             * mSlaveIPGeometryInterpolator  = nullptr;

        // geometry interpolator pointer for the integration cells
        Geometry_Interpolator             * mMasterIGGeometryInterpolator = nullptr;
        Geometry_Interpolator             * mSlaveIGGeometryInterpolator  = nullptr;

        // list of field interpolator pointers
        moris::Cell< Field_Interpolator* >  mMasterFI;
        moris::Cell< Field_Interpolator* >  mSlaveFI;

        // cell of pointers to IWG objects
        moris::Cell< std::shared_ptr< IWG > > mIWGs;
        moris::Cell< moris::Matrix < DDUMat > > mIWGJacDofAssemblyMap;
        moris::Cell< moris::Matrix < DDUMat > > mIWGResDofAssemblyMap;

        enum fem::Element_Type mElementType;

        // map of master and slave dof types for assembly
        moris::Matrix< DDSMat > mDofAssemblyMap;
        uint                    mTotalDof;

        // lists of master and slave groups of dof types
        moris::Cell< moris::Cell< enum MSI::Dof_Type > > mMasterDofTypes;
        moris::Cell< moris::Cell< enum MSI::Dof_Type > > mSlaveDofTypes;

        // maps for the master and slave dof type
        moris::Matrix< DDSMat > mMasterDofTypeMap;
        moris::Matrix< DDSMat > mSlaveDofTypeMap;

        // integration points
        Matrix< DDRMat > mIntegPoints;

        // integration weights
        Matrix< DDRMat > mIntegWeights;

        bool mIsTrivialMaster = false;
        bool mIsTrivialSlave  = false;

        friend class MSI::Equation_Object;
        friend class Cluster;
        friend class Element_Bulk;
        friend class Element_Sideset;
        friend class Element_Double_Sideset;
        friend class Element_Time_Sideset;
        friend class Element;

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
        Set(){};

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
         * create a unique property type list for the set
         * one for the master, one for the slave
         */
        void create_property_type_list();

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
                    MORIS_ERROR(false, "Set::get_dof_type_list - can only be MASTER or SLAVE");
                    return mMasterDofTypes;
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
         * create dof assembly map
         */
        void create_dof_assembly_map();

//------------------------------------------------------------------------------
        /**
         * get dof assembly map
         */
        moris::Matrix< DDSMat > & get_dof_assembly_map()
        {
            return mDofAssemblyMap;
        }

//------------------------------------------------------------------------------
        /**
         * get the total number of dof
         * param [ out ] aTotalNumDof total number of dof
         */
        uint get_total_number_of_dofs()
        {
            return mTotalDof;
        }

//------------------------------------------------------------------------------
         /**
          * create field interpolators for the set
          * @param[ in ] aModelSolverInterface model solver interface
          * ( only used to set the time levels )
          */
         void create_field_interpolators( MSI::Model_Solver_Interface * aModelSolverInterface );

//------------------------------------------------------------------------------
        /**
         * get field interpolators for the set
         * @param[ in ] aIsMaster enum for master or slave
         */
        moris::Cell< Field_Interpolator* > & get_field_interpolators( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
        {
            switch ( aIsMaster )
            {
                case ( mtk::Master_Slave::MASTER ):
                {
                    return mMasterFI;
                }
                case( mtk::Master_Slave::SLAVE ):
                {
                    return mSlaveFI;
                }
                default:
                {
                    MORIS_ERROR(false, "Set::get_field_interpolators - can only be MASTER or SLAVE.");
                    return mMasterFI;
                }
            }
        }

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
         * get IWGs
         * param[ out ] aIWGs cell of IWG pointers
         */
        moris::Cell< std::shared_ptr< IWG > > & get_IWGs()
        {
            return mIWGs;
        }

//------------------------------------------------------------------------------
        /**
         * set the field interpolators for the IWGs
         */
        void set_IWG_field_interpolators();

//------------------------------------------------------------------------------
        /**
         * set the field interpolators for the IWGs
         */
        void set_IWG_geometry_interpolators();

//------------------------------------------------------------------------------
        /**
         * create dof assembly maps for each IWG
         */
        void create_IWG_dof_assembly_map();

//------------------------------------------------------------------------------
        /**
         * get residual dof assembly maps for all IWG
         * param[ out ] aIWGResDofAssemblyMap a map for the residual dof on the IWG
         */
        moris::Cell< Matrix< DDUMat > > & get_IWG_res_dof_assembly_map()
        {
            return mIWGResDofAssemblyMap;
        }

//------------------------------------------------------------------------------
        /**
         * get jacobian dof assembly maps for all IWG
         * param[ out ] aIWGJacDofAssemblyMap a map for the jacobian dof on the IWG
         */
        moris::Cell< Matrix< DDUMat > > & get_IWG_jac_dof_assembly_map()
        {
            return mIWGJacDofAssemblyMap;
        }

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
         * get the IP geometry interpolator
         * @param[ in ]  aIsMaster             enum for master or slave
         * @param[ out ] aGeometryInterpolator geometry interpolator pointer for IP cells
         */
        Geometry_Interpolator * get_IP_geometry_interpolator( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
        {
            switch ( aIsMaster )
            {
                case ( mtk::Master_Slave::MASTER ):
                {
                    return mMasterIPGeometryInterpolator;
                }
                case( mtk::Master_Slave::SLAVE ):
                {
                    return mSlaveIPGeometryInterpolator;
                }
                default:
                {
                    MORIS_ERROR(false, "Set::get_IP_geometry_interpolator - can only be MASTER or SLAVE");
                    return nullptr;
                }
            }
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
         * get the IG geometry interpolator
         * @param[ in ] aIsMaster              enum for master or slave
         * @param[ out ] aGeometryInterpolator geometry interpolator pointer for IG cells
         */
        Geometry_Interpolator * get_IG_geometry_interpolator( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
        {
            switch ( aIsMaster )
            {
                case ( mtk::Master_Slave::MASTER ):
                {
                    return mMasterIGGeometryInterpolator;
                }
                case( mtk::Master_Slave::SLAVE ):
                {
                    return mSlaveIGGeometryInterpolator;
                }
                default:
                {
                    MORIS_ERROR(false, "Set::get_IG_geometry_interpolator - can only be MASTER or SLAVE");
                    return nullptr;
                }
            }
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
         * get the field interpolators for a dof type
         * param[ in ] aDofType  dof type of the field interpolator to grab
         * param[ in ] aIsMaster enum for master or slave
         */
         Field_Interpolator* get_dof_type_field_interpolators( enum MSI::Dof_Type aDofType,
                                                               mtk::Master_Slave  aIsMaster = mtk::Master_Slave::MASTER );

//------------------------------------------------------------------------------
        /**
         * auto detect full integration scheme
         * @param[ in ] aGeometryType geometry type
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

    };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_SET_HPP_ */
