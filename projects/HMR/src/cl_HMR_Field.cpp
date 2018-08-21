/*
 * cl_HMR_Data.cpp
 *
 *  Created on: Jun 25, 2018
 *      Author: messe
 */

#include "op_times.hpp" //LNA/src
#include "cl_HMR_Field.hpp" //HMR/src

namespace moris
{
    namespace hmr
        {
//------------------------------------------------------------------------------

        Field::Field(
                const Parameters     * aParameters,
                const std::string  & aLabel,
                Background_Mesh_Base * aBackgroundMesh,
                BSpline_Mesh_Base  * aBSplineMesh,
                Lagrange_Mesh_Base * aLagrangeMesh ) :
                    mParameters( aParameters ),
                    mLabel( aLabel ),
                    mBackgroundMesh( aBackgroundMesh ),
                    mBSplineMesh( aBSplineMesh ),
                    mLagrangeMesh( aLagrangeMesh ),
                    mTMatrix( aParameters, aBSplineMesh, aLagrangeMesh ),
                    mNumberOfBasis (  aBSplineMesh->get_number_of_active_basis_on_proc() ),
                    mNumberOfNodes ( aLagrangeMesh->get_number_of_nodes_on_proc() )
        {

            // initialize B-Spline values
            mBSplineValues.set_size(
                    mDimension,
                    mNumberOfBasis,
                    0 );

            // initialize Lagrange Values
            mLagrangeValues.set_size(
                    mDimension,
                    mNumberOfNodes );

        }

//------------------------------------------------------------------------------

         void
         Field::calculate_lagrange_values( )
         {
             // unflag all nodes
             for( uint k=0; k<mNumberOfNodes; ++k )
             {
                 mLagrangeMesh->get_basis_by_memory_index( k )->unflag();
             }

             // get number of elements
             luint tNumberOfElements
                 = mBackgroundMesh->get_number_of_active_elements_on_proc();

             // transposed T-Matrix
             Mat< real > tT;

             // lagrange conversion matrix
             Mat< real > tL = mTMatrix.get_lagrange_matrix();

             Cell< Basis* > tDOFs;

             // number of nodes per element
             uint tNumberOfNodesPerElement
                 = mLagrangeMesh->get_number_of_basis_per_element();

             // loop over all active elements
             for( uint e=0; e<tNumberOfElements; ++e )
             {
                 // get pointer to background element
                 Background_Element_Base * tBackElement = mBackgroundMesh->get_element( e );

                 mTMatrix.calculate_truncated_t_matrix(
                         tBackElement->get_memory_index(),
                         tT,
                         tDOFs );


                 // get number of dofs
                 uint tNumberOfDofs = tDOFs.size();

                 // fill matrix with values
                 Mat< real > tBValues ( tNumberOfDofs, 1 );
                 for( uint k=0; k<tNumberOfDofs; ++k )
                 {
                     tBValues( k ) = mBSplineValues ( tDOFs( k )->get_active_index() );
                 }

                 // calculate Lagrange values for this element
                 Mat< real > tLValues = tL * tT * tBValues;

                 /* std::cout << "Element " << tBackElement->get_domain_id() << std::endl;

                 tBValues.print("B");
                 tT.print("T");
                 tL.print("L");


                 tLValues.print("Values"); */

                 // get Lagrange Element
                 Element * tElement
                     = mLagrangeMesh->get_element_by_memory_index( tBackElement->get_memory_index() );

                 // write values into array
                 for( uint k=0; k<tNumberOfNodesPerElement; ++k )
                 {
                     // get node
                     Basis * tBasis = tElement->get_basis( k );

                     // test if basis has been processed
                     if ( ! tBasis->is_flagged() )
                     {
                         // write value
                         mLagrangeValues( tBasis->get_memory_index() )
                                 = tLValues( k );

                         // flag this basis
                         tBasis->flag();
                     }
                 }
             }
         }

//------------------------------------------------------------------------------

         void
         Field::append_to_mtk_object( MTK * aMTK )
         {
             aMTK->add_node_data( mLabel, mLagrangeValues );
         }

        } /* namespace hmr */
} /* namespace moris */
