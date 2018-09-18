/*
 * cl_HMR.hpp
 *
 *  Created on: May 5, 2018
 *      Author: messe
 */

#ifndef SRC_HMR_CL_HMR_HPP_
#define SRC_HMR_CL_HMR_HPP_

#include "cl_Cell.hpp"             //CON/src

#include "cl_HMR_Factory.hpp"        //HMR/src
#include "cl_HMR_Lagrange_Mesh.hpp"  //HMR/src
#include "cl_HMR_Mesh.hpp"      //HMR/src
#include "cl_HMR_Parameters.hpp"     //HMR/src
#include "cl_HMR_T_Matrix.hpp"       //HMR/src

namespace moris
{
    namespace hmr
    {
        class Field;
// -----------------------------------------------------------------------------
        /**
         * \brief the main class of HMR
         */
        class HMR
        {
        public :
            //! object containing user settings
            Parameters *          mParameters;

            //! flag telling if parameter pointer is suppposed to be deleted on destruction
            bool                        mDeleteParametersOnDestruction = false;

            //! pointer to background mesh
            Background_Mesh_Base*       mBackgroundMesh;

            //! cell of pointers to B-Spline meshes
            Cell< BSpline_Mesh_Base* >  mBSplineMeshes;

            //! cell of pointers to Lagrange meshes
            Cell< Lagrange_Mesh_Base* > mLagrangeMeshes;

            //! calculation object that calculates the T-Matrices
            Cell< T_Matrix* >           mTMatrix;

            //! communication table for this mesh. Created during finalize.
            Mat< uint >                 mCommunicationTable;

            //! container with field objects
            //Cell< Field* >              mFields;

// -----------------------------------------------------------------------------
        public :
// -----------------------------------------------------------------------------

            /**
             * the  HMR constructor
             *
             * @param[in] aParameters  ref to container of user defined settings
             */
            HMR ( Parameters * aParameters ) ;

// -----------------------------------------------------------------------------

            /**
             * alternative constructor using a ref
             *
             * @param[in] aParameters  ref to container of user defined settings
             */
            HMR ( Parameters & aParameters ) ;

// -----------------------------------------------------------------------------

            /**
             * alternative constructor using a parameter list
             *
             * @param[in] aParameters  ref to container of user defined settings
             */
            HMR ( ParameterList & aParameterList ) ;
// -----------------------------------------------------------------------------

            /**
             * alternative constructor which loads a mesh from a h5 file
             */
            HMR( const std::string & aPath );

// -----------------------------------------------------------------------------

            /**
             * default destructor of HMR
             */
            ~HMR ( ) ;

// -----------------------------------------------------------------------------

            /**
             * exposes the parameters pointer
             */
            Parameters *
            get_parameters()
            {
                return mParameters;
            }

// -----------------------------------------------------------------------------

            /**
             *  this function updates the meshes after an refinement step
             */
            void
            update_meshes();

// -----------------------------------------------------------------------------

            /**
             * returns the number of ( active ) elements on this proc
             */
            auto
            get_number_of_elements_on_proc()
                -> decltype( mBackgroundMesh->get_number_of_active_elements_on_proc() )
            {
                return mBackgroundMesh->get_number_of_active_elements_on_proc();
            }

// -----------------------------------------------------------------------------

             /**
              * returns the number of dimensions in space
              */
             auto
             get_number_of_dimensions() const
                 -> decltype( mParameters->get_number_of_dimensions() )
             {
                 return mParameters->get_number_of_dimensions();
             }

// -----------------------------------------------------------------------------

             /**
              * returns the number of Lagrange meshes
              */
             uint
             get_number_of_lagrange_meshes() const
             {
                 return mLagrangeMeshes.size();
             }

// -----------------------------------------------------------------------------

             /**
              * returns the number of connected fields
              */
             /*uint
             get_number_of_fields() const
             {
                 return mFields.size();
             } */

// -----------------------------------------------------------------------------

             /**
              * returns the pointer to a Lagrange mesh, needed by interface
              * constructor
              */
             Lagrange_Mesh_Base*
             get_lagrange_mesh_by_index( const uint& aIndex )
             {
                 return mLagrangeMeshes( aIndex );
             }

// -----------------------------------------------------------------------------

             /**
              * Creates an STK interface object. Per default,  the output
              * pattern is selected
              */
             Mesh *
             create_mtk_interface();
// -----------------------------------------------------------------------------

             /**
              * Creates an STK interface object with respect to a specified
              * output pattern. Used internally for L2 projection.
              */
             Mesh *
             create_mtk_interface(  const uint & aActivationPattern );

//-----------------------------------------------------------------------------

             /**
              * creates an MTK pointer to the input mesh
              */
             Mesh *
             create_input_mesh();

//-----------------------------------------------------------------------------

             /**
              * creates an MTK  pointer to the input mesh
              */
             Mesh *
             create_output_mesh();

//-----------------------------------------------------------------------------

             /**
              * popules the member variables of the relevant nodes
              * with their T-Matrices
              */
             void
             finalize();

//------------------------------------------------------------------------------
             /**
              * Node IDs are calculated with respect to used T-Matrices.
              * This function activates the T-Matrices of all active elements.
              * Needed for STK output if not connected to FEM module
              */
             void
             activate_all_t_matrices();

//------------------------------------------------------------------------------

             /**
              * provides a moris::Mat<uint> containing the IDs this mesh has
              * to communicate with
              */
             Mat< uint >
             get_communication_table() const
             {
                 return mCommunicationTable;
             }

// -----------------------------------------------------------------------------

             /**
              * Temporary function to add field data to the mesh object.
              * Needed for testing.
              */
             /* Field *
             create_field(
                     const std::string & aLabel,
                     const uint        & aLagrangeIndex ) */

// -----------------------------------------------------------------------------

             /**
              * experimental funciton
              */
             void
             flag_element( const uint & aIndex )
             {
                 mBackgroundMesh->get_element( aIndex )->put_on_refinement_queue();
             }

// -----------------------------------------------------------------------------

             /**
              * set active pattern of background mesh
              */
             void
             set_activation_pattern( const uint & aPattern )
             {
                 mBackgroundMesh->set_activation_pattern( aPattern );
             }

// -----------------------------------------------------------------------------

             /**
              * returns the active pattern
              */
             auto
             get_activation_pattern() const
                 -> decltype( mBackgroundMesh->get_activation_pattern() )
             {
                 return  mBackgroundMesh->get_activation_pattern();
             }

// -----------------------------------------------------------------------------

             /**
              * creates a union of two patterns
              */
             void
             unite_patterns(
                     const uint & aSourceA,
                     const uint & aSourceB,
                     const uint & aTarget );

// -----------------------------------------------------------------------------

             /**
              * copies a source pattern to a target pattern
              */
             void
             copy_pattern(
                     const uint & aSource,
                     const uint & aTarget );

// -----------------------------------------------------------------------------

             /**
              * runs the refinement scheme
              */
             void
             perform_refinement();

// -----------------------------------------------------------------------------

             // fixme: this function needs to be moved
             void
             save_to_exodus( const uint & aBlock, const std::string & aPath );

// -----------------------------------------------------------------------------

             // fixme: this function needs to be moved
             void
             save_to_exodus( const std::string & aPath );

// -----------------------------------------------------------------------------

             void
             save_to_hdf5( const std::string & aPath );

// -----------------------------------------------------------------------------

             // fixme: this function needs to be moved
             /**
              * aTarget must be a refined variant of aSource
              */
             void
             interpolate_field(
                     const uint       & aSourcePattern,
                     const mtk::Field * aSource,
                     const uint       & aTargetPattern,
                           mtk::Field * aTarget );

// -----------------------------------------------------------------------------

             /**
              * returns the pointer to a T-Matrix object. Needed by field
              */
             T_Matrix *
             get_t_matrix( const uint & aLagrangeMeshIndex );

// -----------------------------------------------------------------------------

             /**
              * This function checks if t-matrix flags on the neighbor procs
              * have been set. If so, it makes sure that basis owned by current
              * proc are created
              */
             void
             synchronize_t_matrix_flags();

// -----------------------------------------------------------------------------

             /**
              * Project an input field to the output mesh
              */
             /* Field *
             map_field_to_output_mesh(  Field * aSource ); */

// -----------------------------------------------------------------------------

             /**
              * Project an input field to the output mesh, added function
              * and error for testing purpose
              */
             /* Field *
             map_field_to_output_mesh(
                     Field * aSource,
                     real & aIntegrationError,
                     real (*aFunction)( const Mat< real > & aPoint ) ); */

// -----------------------------------------------------------------------------

             /**
              * needed for exodus output of cubic meshes, called by finalize
              */
             void
             add_extra_refinement_step_for_exodus();

// -----------------------------------------------------------------------------

             /**
              * add field pointer to internal list
              */
             /* void
             add_field( Field * aField ); */

// -----------------------------------------------------------------------------

             /**
              * Extract values from source and copy them to target.
              * Needed for testing
              * aSource must be a refined variant of aTarget
              */
             /* void
             extract_field( Field * aSource, Field* aTarget ); */

// -----------------------------------------------------------------------------

             /**
              * fixme: obsolete: remove this
              *
              * calls the refinement manager and refines against a given
              * nodal field
              *
              * @param[ in ] aNodalValues        Nodal field with data
              *
              */
             void
             flag_against_nodal_field(
                     const Mat< real > & aNodalValues );

// -----------------------------------------------------------------------------

             /**
              * flags active elements
              *
              * @param[ in ]   aElementIndices  indices of elements to be considered
              * @param[ in ]   aLevelLowPass    only elements below this level are considered
              * @param[ in ]   aPattern         choose activation for processing ( default: output )
              */
             void
             flag_elements(
                     const Mat< moris_index >   & aElementIndices,
                     const uint                 aLevelLowPass = MORIS_UINT_MAX,
                     const uint                 aPattern      = MORIS_UINT_MAX );


// -----------------------------------------------------------------------------

             mtk::Field *
             map_field_on_mesh( mtk::Field * aField, Mesh* aMesh );

// -----------------------------------------------------------------------------

             /**
              * creates a union mesh of the input and the output patterns
              */
             void
             create_union_pattern()
             {
                 this->unite_patterns(
                       mParameters->get_input_pattern(),
                       mParameters->get_output_pattern(),
                       mParameters->get_union_pattern() );
             }

// -----------------------------------------------------------------------------
        private:
// -----------------------------------------------------------------------------

            /**
             * this function initializes the Lagrange and B-Spline Meshes
             * is complete
             *
             * @return void
             */
            void
            create_meshes();

// -----------------------------------------------------------------------------
            /**
             * this function deletes the Lagrange and B-Spline meshes
             * the function is called before create_meshes
             */
            void
            delete_meshes();
// -----------------------------------------------------------------------------

            /**
             * initializes the T-Matrix objects
             */
            void
            init_t_matrices();

// -----------------------------------------------------------------------------

            /**
             * deletes the T-Matrix objects
             */
            void
            delete_t_matrices();

// -----------------------------------------------------------------------------

            /**
             * creates the communication table and writes it into
             * mCommunicationTable. Must be called after mesh has been finalized.
             */
            void
            create_communication_table();


        }; /* HMR */

    } /* namespace hmr */
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_HPP_ */
