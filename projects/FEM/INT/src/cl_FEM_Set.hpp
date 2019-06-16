/*
 * cl_FEM_Set.hpp
 *
 *  Created on: Mar 10, 2019
 *      Author: schmidt
 */

#ifndef SRC_FEM_CL_FEM_SET_HPP_
#define SRC_FEM_CL_FEM_SET_HPP_

#include "assert.h"
#include "cl_MTK_Enums.hpp"            //MTK
#include "cl_MSI_Equation_Set.hpp"     //FEM/MSI/src
#include "cl_FEM_Enums.hpp"            //FEM/INT/src
#include "cl_FEM_Node_Base.hpp"        //FEM/INT/src
#include "cl_Communication_Tools.hpp"

#include "cl_MTK_Cell_Cluster.hpp"
#include "cl_MTK_Side_Cluster.hpp"
#include "cl_MTK_Double_Side_Cluster.hpp"

namespace moris
{
namespace mtk
{
   class Cell;
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
        // if block-set
        moris::Cell< mtk::Cluster const* > mMeshClusterList;
//        // if block-set
//        moris::Cell< mtk::Cell_Cluster const* > mMeshCellClusterList;
//        // if side-set
//        moris::Cell< mtk::Side_Cluster const* > mMeshSideClusterList;
//        // if double side-set
//        moris::Cell< mtk::Double_Side_Cluster > const mMeshDoubleSideClusterList;

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
        Geometry_Interpolator             * mIPGeometryInterpolator = nullptr;

        // geometry interpolator pointer for the integration cells
        Geometry_Interpolator             * mIGGeometryInterpolator = nullptr;

        // list of field interpolator pointers
        moris::Cell< Field_Interpolator* >  mFieldInterpolators;

        // list of field interpolator pointers for double side-set
        moris::Cell< Field_Interpolator* >  mLeftFieldInterpolators;
        moris::Cell< Field_Interpolator* >  mRightFieldInterpolators;

        // geometry interpolator pointer for the interpolation cells for double side-set
        Geometry_Interpolator             * mLeftIPGeometryInterpolator = nullptr;
        Geometry_Interpolator             * mRightIPGeometryInterpolator = nullptr;

        // geometry interpolator pointer for the integration cells
        Geometry_Interpolator             * mLeftIGGeometryInterpolator = nullptr;
        Geometry_Interpolator             * mRightIGGeometryInterpolator = nullptr;

        // cell of pointers to IWG objects
        moris::Cell< IWG* > mIWGs;

        enum fem::Element_Type mElementType;

        // map of the element active dof types
        moris::Matrix< DDSMat >                   mInterpDofTypeMap;
        moris::Cell< moris::Cell< enum MSI::Dof_Type > > mInterpDofTypeList;
        moris::Matrix< DDSMat >                   mInterpDofAssemblyMap;
        uint                                      mTotalDof;
        uint                                      mNumOfInterp;

        // list of model parameter type for the set
        moris::Cell< fem::Property_Type > mPropertyTypeList;

        // number of integration points
        uint mNumOfIntegPoints;

        // integration points
        Matrix< DDRMat > mIntegPoints;

        // integration weights
        Matrix< DDRMat > mIntegWeights;

        friend class MSI::Equation_Object;
        friend class Element_Bulk;
        friend class Element_Sideset;
        friend class Element_Double_Sideset;
        friend class Element_Time_Sideset;
        friend class Element;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------
                /**
                 * constructor for block-set
                 * @param[ in ]     List of mtk::Cell_Cluster
                 */
          Set( moris::Cell< mtk::Cluster const * > & aMeshClusterList,
               Element_Type                               aElementType,
               moris::Cell< IWG* >                      & aIWGs,
               moris::Cell< Node_Base* >                & aIPNodes);
//        /**
//         * constructor for block-set
//         * @param[ in ]     List of mtk::Cell_Cluster
//         */
//        Set( moris::Cell< mtk::Cell_Cluster const * > & aMeshClusterList,
//             Element_Type                               aElementType,
//             moris::Cell< IWG* >                      & aIWGs,
//             moris::Cell< Node_Base* >                & aIPNodes);
//
//        /**
//         * constructor for side-set
//         * @param[ in ]     List of mtk::Side_Cluster
//         */
//        Set( moris::Cell< mtk::Side_Cluster const * > & aMeshClusterList,
//             Element_Type                               aElementType,
//             moris::Cell< IWG* >                      & aIWGs,
//             moris::Cell< Node_Base* >                & aIPNodes);
//
//        /**
//         * constructor for double side-set
//         * @param[ in ]     List of mtk::Double_Side_Cluster
//         */
//        Set( moris::Cell< mtk::Double_Side_Cluster > const & aMeshClusterList,
//             Element_Type                                    aElementType,
//             moris::Cell< IWG* >                           & aIWGs,
//             moris::Cell< Node_Base* >                     & aIPNodes);
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

        void delete_pointers();

//------------------------------------------------------------------------------

        void finalize( MSI::Model_Solver_Interface * aModelSolverInterface );

//------------------------------------------------------------------------------

        void create_dof_type_lists();

//------------------------------------------------------------------------------

        void create_unique_dof_type_lists();

//------------------------------------------------------------------------------

        void create_unique_property_type_list();

//------------------------------------------------------------------------------

        void create_dof_assembly_map();

        void create_dof_assembly_map_double();

//------------------------------------------------------------------------------

        uint get_num_equation_objects()
        {
            return mEquationObjList.size();
        };

//------------------------------------------------------------------------------

        moris::Cell< MSI::Equation_Object * > & get_equation_object_list()
        {
            return mEquationObjList;
        };

//------------------------------------------------------------------------------

        moris::Cell< Field_Interpolator* > & get_field_interpolator()
        {
            return mFieldInterpolators;
        }

        moris::Cell< Field_Interpolator* > & get_left_field_interpolator()
        {
            return mLeftFieldInterpolators;
        }
        moris::Cell< Field_Interpolator* > & get_right_field_interpolator()
        {
            return mRightFieldInterpolators;
        }
//------------------------------------------------------------------------------

        Geometry_Interpolator * get_IP_geometry_interpolator()
        {
            return mIPGeometryInterpolator;
        }

        Geometry_Interpolator * get_IG_geometry_interpolator()
        {
            return mIGGeometryInterpolator;
        }

        Geometry_Interpolator * get_left_IP_geometry_interpolator()
        {
            return mLeftIPGeometryInterpolator;
        }
        Geometry_Interpolator * get_right_IP_geometry_interpolator()
        {
            return mRightIPGeometryInterpolator;
        }
        Geometry_Interpolator * get_left_IG_geometry_interpolator()
        {
            return mLeftIGGeometryInterpolator;
        }
        Geometry_Interpolator * get_right_IG_geometry_interpolator()
        {
            return mRightIGGeometryInterpolator;
        }

//------------------------------------------------------------------------------

        mtk::Geometry_Type get_IG_geometry_type()
        {
            return mIGGeometryType;
        }

//------------------------------------------------------------------------------
        mtk::Interpolation_Order get_IG_space_interpolation_order()
        {
            return mIGSpaceInterpolationOrder;
        }
//------------------------------------------------------------------------------

        uint get_num_IWG()
        {
            return mIWGs.size();
        }

//------------------------------------------------------------------------------

        moris::Cell< IWG* > & get_IWGs()
        {
            return mIWGs;
        }

//------------------------------------------------------------------------------

        moris::Cell< enum MSI::Dof_Type > & get_unique_dof_type_list()
        {
            return mEqnObjDofTypeList;
        }

//------------------------------------------------------------------------------

        moris::Matrix< DDSMat > & get_interpolator_dof_type_map()
        {
            return mInterpDofTypeMap;
        }

//------------------------------------------------------------------------------

        moris::Cell< moris::Cell< enum MSI::Dof_Type > > & get_interpolator_dof_type_list()
        {
            return mInterpDofTypeList;
        }

//------------------------------------------------------------------------------

        moris::Matrix< DDSMat > & get_interpolator_dof_assembly_map()
        {
            return mInterpDofAssemblyMap;
        }
//------------------------------------------------------------------------------

        uint get_total_number_of_dofs()
        {
            return mTotalDof;
        }

//------------------------------------------------------------------------------

        uint & get_num_interpolators()
        {
            return mNumOfInterp;
        }

//------------------------------------------------------------------------------

        uint get_num_integration_points()
        {
            return mNumOfIntegPoints;
        }

//------------------------------------------------------------------------------

        const Matrix< DDRMat > & get_integration_points()
        {
            return mIntegPoints;
        }

//------------------------------------------------------------------------------

        const Matrix< DDRMat > & get_integration_weights()
        {
            return mIntegWeights;
        }

//------------------------------------------------------------------------------

        /**
         * get the field interpolators for an IWG
         */
        moris::Cell< Field_Interpolator* > get_IWG_field_interpolators ( IWG*                               & aIWG,
                                                                         moris::Cell< Field_Interpolator* > & aFieldInterpolators );

//------------------------------------------------------------------------------

        /**
         * auto detect full integration scheme
         */
        //FIXME: works for Lagrange only
        mtk::Interpolation_Order get_auto_interpolation_order( const moris::uint aNumVertices,
                                                               const mtk::Geometry_Type aGeometryType );

//------------------------------------------------------------------------------

        fem::Interpolation_Type get_auto_time_interpolation_type( const moris::uint aNumVertices );

//------------------------------------------------------------------------------

        mtk::Geometry_Type get_auto_side_geometry_type( const mtk::Geometry_Type aGeometryType );

//------------------------------------------------------------------------------

        void create_field_interpolators( MSI::Model_Solver_Interface * aModelSolverInterface );

        void create_field_interpolators_double( MSI::Model_Solver_Interface * aModelSolverInterface );

//------------------------------------------------------------------------------

        /**
          * auto detect interpolation scheme
          */
        fem::Integration_Order get_auto_integration_order( const mtk::Geometry_Type aGeometryType );

//------------------------------------------------------------------------------

        void initialize_mJacobian();

//------------------------------------------------------------------------------

        void initialize_mResidual();

//------------------------------------------------------------------------------

    };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_SET_HPP_ */
