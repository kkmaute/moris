/*
 * cl_FEM_Set.hpp
 *
 *  Created on: Mar 10, 2019
 *      Author: schmidt
 */

#ifndef SRC_FEM_CL_FEM_ELEMENT_BLOCK_HPP_
#define SRC_FEM_CL_FEM_ELEMENT_BLOCK_HPP_

#include "assert.h"
#include "cl_MSI_Equation_Set.hpp"               //FEM/INT/src
#include "cl_FEM_Enums.hpp"               //FEM/INT/src
#include "cl_MTK_Enums.hpp"               //FEM/INT/src
#include "cl_FEM_Node_Base.hpp"               //FEM/INT/src

#include "cl_Communication_Tools.hpp"               //FEM/INT/src


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
        moris::Cell< mtk::Cell const* >     mMeshElementPointer;

        moris::Cell< Node_Base* >           mNodes;

        Geometry_Interpolator             * mGeometryInterpolator = nullptr;

        moris::Cell< Field_Interpolator* >  mFieldInterpolators;

        // cell of pointers to IWG objects
        moris::Cell< IWG* > mIWGs;

        enum fem::Element_Type mElementType;

        // map of the element active dof types
        moris::Matrix< DDSMat >                   mInterpDofTypeMap;
        moris::Cell< Cell< enum MSI::Dof_Type > > mInterpDofTypeList;
        moris::Matrix< DDSMat >                   mInterpDofAssemblyMap;
        uint                                      mTotalDof;
        uint                                      mNumOfInterp;

        // number of integration points
        uint mNumOfIntegPoints;

        // integration points
        Matrix< DDRMat > mSurfRefIntegPoints;

        // integration weights
        Matrix< DDRMat > mIntegWeights;

        friend class MSI::Equation_Object;
        friend class Element_Bulk;
        friend class Element_Sideset;
        friend class Element_Time_Sideset;
        friend class Element;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------
        /**
         * constructor
         *
         * @param[ in ]     List of mtk::Cell pointer
         */
        Set( moris::Cell< mtk::Cell const * > & aCell,
                       Element_Type                aElementType,
                       Cell< IWG* >              & aIWGs,
                       Cell< Node_Base* >        & aNodes);

        /**
         * trivial constructor
         */
        Set()
        {};

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

        void create_dof_assembly_map();

//------------------------------------------------------------------------------

        uint get_num_equation_objects()
        {
            return mEquationObjList.size();
        };

//------------------------------------------------------------------------------

        Cell< MSI::Equation_Object * > & get_equation_object_list()
        {
            return mEquationObjList;
        };

//------------------------------------------------------------------------------

        moris::Cell< Field_Interpolator* > & get_block_field_interpolator()
        {
            return mFieldInterpolators;
        }

//------------------------------------------------------------------------------

        Geometry_Interpolator * get_block_geometry_interpolator()
        {
            return mGeometryInterpolator;
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

        moris::Cell< Cell< enum MSI::Dof_Type > > & get_interpolator_dof_type_list()
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
            return mSurfRefIntegPoints;
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

        void create_field_interpolators( MSI::Model_Solver_Interface * aModelSolverInterface );

//------------------------------------------------------------------------------

        /**
          * auto detect interpolation scheme
          */
        fem::Integration_Order get_auto_integration_order( const mtk::Geometry_Type aGeometryType );

//------------------------------------------------------------------------------

        void initialize_mJacobianElement();

//------------------------------------------------------------------------------

        void initialize_mResidualElement();

//------------------------------------------------------------------------------

    };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_ELEMENT_BLOCK_HPP_ */
