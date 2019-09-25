/*
 * cl_GE_SDF.hpp
 *
 *  Created on: Dec 28, 2018
 *      Author: sonne
 */

#ifndef PROJECTS_MTK_GE_SRC_CL_GE_SDF_HPP_
#define PROJECTS_MTK_GE_SRC_CL_GE_SDF_HPP_

// GE includes
#include "cl_GE_Geometry.hpp"

// HMR includes
#include "cl_HMR_Field.hpp"

// SDF includes
#include "../projects/GEN/SDF/src/cl_SDF_Generator.hpp"

// other includes
#include "assert.hpp"

namespace moris
{
namespace ge
{
	class SDF : public Geometry
	{
	public:
	    /*
	     * @brief signed distance field geometry representation class
	     */
		SDF()
	    {
            mMySpaceInterpType  = fem::Interpolation_Type::LAGRANGE;
            mMySpaceInterpOrder = mtk::Interpolation_Order::LINEAR;
            mMyTimeInterpType   = fem::Interpolation_Type::CONSTANT;
            mMyTimeInterpOrder  = mtk::Interpolation_Order::CONSTANT;
	    };
		~SDF()
		{};

        //------------------------------------------------------------------------------
        /*
         * @brief function to report geometry representation type
         */
        enum GeomType
        get_geom_type() const
        {
            return GeomType::SDF;
        }
        //------------------------------------------------------------------------------
        void
        set_my_mesh(mtk::Mesh_Manager* aMyMesh)
        {
            //fixme set the hmr fields/add_hmr_field(), at this point as well?
            mMyMesh = aMyMesh;
        }
        //------------------------------------------------------------------------------
        void
        set_my_interpolation_rules( fem::Interpolation_Type  aSpaceInterpType,
                                    mtk::Interpolation_Order aSpaceInterpOrder,
                                    fem::Interpolation_Type  aTimeInterpType  = fem::Interpolation_Type::CONSTANT,
                                    mtk::Interpolation_Order aTimeInterpOrder = mtk::Interpolation_Order::CONSTANT )
        {
            mMySpaceInterpType  = aSpaceInterpType;
            mMySpaceInterpOrder = aSpaceInterpOrder;
            mMyTimeInterpType   = aTimeInterpType;
            mMyTimeInterpOrder  = aTimeInterpOrder;
        }
        //------------------------------------------------------------------------------
        fem::Interpolation_Type
        get_my_space_interpolation_type()
        {
            return mMySpaceInterpType;
        }
        //------------------------------------------------------------------------------
        mtk::Interpolation_Order
        get_my_space_interpolation_order()
        {
            return mMySpaceInterpOrder;
        }
        //------------------------------------------------------------------------------
        fem::Interpolation_Type
        get_my_time_interpolation_type()
        {
            return mMyTimeInterpType;
        }
        //------------------------------------------------------------------------------
        mtk::Interpolation_Order
        get_my_time_interpolation_order()
        {
            return mMyTimeInterpOrder;
        }
        //------------------------------------------------------------------------------
        mtk::Mesh_Manager*
        get_my_mesh()
        {
            MORIS_ASSERT( mMyMesh != nullptr, "ge::Geometry::get_my_mesh(): the associated mesh has not been set" );
            return mMyMesh;
        }
        //------------------------------------------------------------------------------
        moris_index
        add_hmr_field( std::shared_ptr< hmr::Field > &aField )
        {
            mHMRFields.push_back( aField );
            return mHMRFields.size()-1;
        }
        //------------------------------------------------------------------------------
        uint get_number_of_sub_types()
        {
            return mHMRFields.size();
        }
        //------------------------------------------------------------------------------
        /*
         * @brief performs all steps necessary to create the signed distance field
         */
        void
        initialize_sdf( const std::string & aObjectPath,
                        std::shared_ptr< mtk::Mesh > aMTKMesh,
                        const bool aVerboseFlag = true,
                        const uint aFieldIndex = 0 )
        {
            sdf::SDF_Generator tSdfGen( aObjectPath );
            mMySDFGenerator = & tSdfGen;

            Matrix< IndexMat > tSurfaceElements;
            //fixme need to know when this is not going to be the 0th interpolation matrix in the manager
//            mMySDFGenerator->raycast( mMyMesh->get_interpolation_mesh(0), tSurfaceElements );
            mMySDFGenerator->raycast( aMTKMesh, tSurfaceElements );
            // hmr mesh implementation only assert
//            MORIS_ASSERT( mMyMesh->get_interpolation_mesh(0)->get_mesh_type() == MeshType::HMR, "ge::SDF::initialize_sdf(): currently only have hmr implementation" );
            //------------------HMR MESH IMPLEMENTATION------------------
            // use the created field to calculate the sdf, an SDF field needs to have been created on the mesh
            //fixme implement a check to make sure there is a field at this index and that it is the right one
//            tSdfGen.calculate_sdf( mMyMesh->get_interpolation_mesh(0), mMyHMRFields(aFieldIndex)->get_node_values() );

            tSdfGen.calculate_sdf( aMTKMesh, mHMRFields(aFieldIndex)->get_node_values() );
            //-----------------------------------------------------------
            mMySDFVals = mHMRFields(aFieldIndex)->get_node_values();
        }
        //------------------------------------------------------------------------------
        /*
         * @brief return the matrix with all the computed sdf values
         */
        real
        get_sdf_vals( moris_index aIndex )
        {
            return mMySDFVals( aIndex );
        }

    private:
        sdf::SDF_Generator* mMySDFGenerator = nullptr;
        mtk::Mesh_Manager*  mMyMesh = nullptr;

        Matrix< DDRMat > mMySDFVals;

        moris::Cell< std::shared_ptr< hmr::Field > > mHMRFields;

        fem::Interpolation_Type  mMySpaceInterpType;
        mtk::Interpolation_Order mMySpaceInterpOrder;
        fem::Interpolation_Type  mMyTimeInterpType;
        mtk::Interpolation_Order mMyTimeInterpOrder;

    protected:

	};
} /* namespace gen */
} /* namespace moris */


#endif /* PROJECTS_MTK_GE_SRC_CL_GE_SDF_HPP_ */
