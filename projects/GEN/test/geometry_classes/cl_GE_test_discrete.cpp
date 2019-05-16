/*
 * cl_GE_test_discrete.cpp
 *
 *  Created on: Apr 22, 2019
 *      Author: sonne
 */

#include "cl_GE_Core.hpp"
#include "catch.hpp"

//------------------------------------------------------------------------------
// GE includes
#include "cl_GE_Element.hpp"
#include "cl_GE_Factory.hpp"
#include "cl_GE_Node.hpp"

// HMR includes
#include "cl_HMR.hpp"
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Mesh.hpp"
#include "cl_HMR_Parameters.hpp"
#include "fn_HMR_Exec_perform_mapping.hpp"

// linalg includes
#include "cl_Matrix.hpp"
#include "fn_all_true.hpp"
#include "fn_equal_to.hpp"
#include "linalg_typedefs.hpp"
#include "op_equal_equal.hpp"

// MTK includes
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Vertex.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"
//------------------------------------------------------------------------------

using namespace moris;
using namespace ge;

real
CircleFunction( const Matrix< DDRMat > & aPoint )
{
    return std::pow(aPoint(0),2) + std::pow(aPoint(1),2) - std::pow(0.6,2);
}

TEST_CASE("discrete_functionalities_test_01","[GE],[discrete_functionalities_hmr_mesh]")
{
    if(par_size()<=1)
    {
        /*
         * *create an hmr mesh to be used
         * *B-spline order = Lagrange Order = 2
         */
        hmr::ParameterList tParameters = hmr::create_hmr_parameter_list();
        tParameters.set( "number_of_elements_per_dimension", "2, 2" );
        tParameters.set( "bspline_orders", "2" );
        tParameters.set( "lagrange_orders", "2" );
        tParameters.set( "verbose", 0 );                    // displays debug info when turned on
        tParameters.set( "domain_dimensions", "2, 2" );     // if no dimensions set, hmr assumes 1x1
        tParameters.set( "domain_offset", "-1, -1" );       // offset so that the center node is at ( 0, 0 )

        hmr::HMR tHMR( tParameters );
        /*
         * *create a field and add to the hmr mesh
         */
        std::shared_ptr< hmr::Field > tField = tHMR.create_field( "circle" );

//        Matrix< DDRMat > tDiscreteNodalVals( 1,25 );    // values of circle LS at nodes, centered at origin with r=0.6
//        tDiscreteNodalVals(0) = 1.64;
//        tDiscreteNodalVals(9) = 1.64;
//        tDiscreteNodalVals(21) = 1.64;
//        tDiscreteNodalVals(16) = 1.64;
//
//        tDiscreteNodalVals(1) = 0.64;
//        tDiscreteNodalVals(10) = 0.64;
//        tDiscreteNodalVals(15) = 0.64;
//        tDiscreteNodalVals(3) = 0.64;
//
//        tDiscreteNodalVals(4) = 0.89;
//        tDiscreteNodalVals(11) = 0.89;
//        tDiscreteNodalVals(12) = 0.89;
//        tDiscreteNodalVals(22) = 0.89;
//        tDiscreteNodalVals(23) = 0.89;
//        tDiscreteNodalVals(18) = 0.89;
//        tDiscreteNodalVals(19) = 0.89;
//        tDiscreteNodalVals(7) = 0.89;
//
//        tDiscreteNodalVals(2) = -0.36;
//        tField->put_scalar_values_on_field( tDiscreteNodalVals );
        tField->evaluate_scalar_function(CircleFunction);

        mtk::Interpolation_Mesh* tInterpMesh1 = mtk::create_interpolation_mesh( MeshType::HMR, tMeshData );
        mtk::Integration_Mesh*   tIntegMesh1  = mtk::create_integration_mesh_from_interpolation_mesh(MeshType::HMR,tInterpMesh1);
        mtk::Mesh_Manager tMeshManager;
        uint tMeshIndex = tMeshManager.register_mesh_pair(tInterpMesh1,tIntegMesh1);


//        tHMR.save_to_exodus( "circle01.exo" );
    }
}
