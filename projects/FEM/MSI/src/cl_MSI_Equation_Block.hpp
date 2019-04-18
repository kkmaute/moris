/*
 * cl_MSI_Equation_Block_Block.hpp
 *
 *  Created on: Apr 10, 2019
 *      Author: schmidt
 */

#ifndef SRC_FEM_CL_MSI_EQUATION_BLOCK_HPP_
#define SRC_FEM_CL_MSI_EQUATION_BLOCK_HPP_

#include "assert.h"
//#include "cl_MSI_Equation_Object.hpp"               //FEM/INT/src

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
    class Equation_Object;
    enum class Dof_Type;
//------------------------------------------------------------------------------
    /**
     * \brief element block class that communicates with the mesh interface
     */
    class Equation_Block
    {
    protected:
        Cell< MSI::Equation_Object* > mElements;

        // map of the element active dof types
        moris::Cell< enum MSI::Dof_Type >         mEqnObjDofTypeList; // List of dof types of this equation obj

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------
        /**
         * constructor
         *
         * @param[ in ]     List of mtk::Cell pointer
         */
//        Equation_Block( moris::Cell< mtk::Cell* > & aCell,
//                       Element_Type                aElementType,
//                       Cell< IWG* >              & aIWGs,
//                       Cell< Node_Base* >        & aNodes){};

        /**
         * trivial constructor
         */
        Equation_Block( ){};

//------------------------------------------------------------------------------

        /**
         * trivial destructor
         */
        virtual ~Equation_Block(){};

//------------------------------------------------------------------------------

//        void delete_pointers();

//------------------------------------------------------------------------------

        virtual void finalize( MSI::Model_Solver_Interface * aModelSolverInterface )
        {
            MORIS_ERROR(false,"Equation_Block::finalize(), not implemented");
        };

//------------------------------------------------------------------------------

//        void create_dof_type_lists();

//------------------------------------------------------------------------------

//        void create_unique_dof_type_lists();

//------------------------------------------------------------------------------

//        void create_unique_list_of_first_dof_type_of_group();

//------------------------------------------------------------------------------

        Cell< MSI::Equation_Object * > & get_equation_object_list()
        {
            return mElements;
        };

//------------------------------------------------------------------------------

//        moris::Cell< Field_Interpolator* > & get_block_field_interpolator()
//        {
//            return mFieldInterpolators;
//        }

//------------------------------------------------------------------------------

//        Geometry_Interpolator * get_block_geometry_interpolator()
//        {
//            return mGeometryInterpolator;
//        }

//------------------------------------------------------------------------------

//        uint get_num_IWG()
//        {
//            return mIWGs.size();
//        }

//------------------------------------------------------------------------------

//        moris::Cell< IWG* > & get_IWGs()
//        {
//            return mIWGs;
//        }

//------------------------------------------------------------------------------

        moris::Cell< enum MSI::Dof_Type > & get_unique_dof_type_list()
        {
            return mEqnObjDofTypeList;
        }

//------------------------------------------------------------------------------

//        moris::Matrix< DDSMat > & get_interpolator_dof_type_map()
//        {
//            return mInterpDofTypeMap;
//        }

//------------------------------------------------------------------------------

//        moris::Cell< Cell< enum MSI::Dof_Type > > & get_interpolator_dof_type_list()
//        {
//            return mInterpDofTypeList;
//        }

//------------------------------------------------------------------------------

//        uint & get_num_interpolators()
//        {
//            return mNumOfInterp;
//        }

//------------------------------------------------------------------------------

        /**
         * auto detect full integration scheme
         */
        //FIXME: works for Lagrange only
//        mtk::Interpolation_Order get_auto_interpolation_order( const moris::uint aNumVertices,
//                                                               const mtk::Geometry_Type aGeometryType );

//------------------------------------------------------------------------------

//        fem::Interpolation_Type get_auto_time_interpolation_type( const moris::uint aNumVertices );

//------------------------------------------------------------------------------

//        void create_field_interpolators( MSI::Model_Solver_Interface * aModelSolverInterface );

    };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_MSI_EQUATION_BLOCK_HPP_ */
