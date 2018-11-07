/*
 * cl_HMR.hpp
 *
 *  Created on: May 5, 2018
 *      Author: messe
 */

#ifndef SRC_HMR_CL_HMR_HPP_
#define SRC_HMR_CL_HMR_HPP_


#include "cl_HMR_Parameters.hpp"     //HMR/src
#include "cl_HMR_Database.hpp"     //HMR/src
#include "cl_HMR_Mesh.hpp"
#include "cl_HMR_Element.hpp"
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

            std::shared_ptr< Database > mDatabase;

            //! flag telling if perform_refinement() has been called
            bool                        mPerformRefinementCalled = false;
            bool                        mUpdateRefinementCalled = false;

            //! mesh which points to input pattern
            Cell< std::shared_ptr< Mesh > > mInputMeshes;

            //! mesh which points to output pattern
            Cell< std::shared_ptr< Mesh > > mOutputMeshes;

            //! container with field objects
            Cell< std::shared_ptr< Field > > mFields;

            //! map for Lagrange orders
            Matrix< DDUMat > mLagrangeOrderToInputMeshIndexMap;

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
            ~HMR ( ){};

// -----------------------------------------------------------------------------

            void
            load_output_pattern_from_path( const std::string & aPath );

// -----------------------------------------------------------------------------
            /**
             * save the mesh to an exodus file
             */
            void
            save_to_exodus(
                    const uint        & aMeshIndex,
                    const std::string & aPath,
                    const double aTimeStep = 0.0 );

// -----------------------------------------------------------------------------

            /**
             * save the mesh to an exodus file
             */
            void
            save_to_exodus(
                    const std::string & aPath,
                    const double aTimeStep = 0.0 );

// -----------------------------------------------------------------------------

            /**
             * save the mesh to an exodus file
             */
            void
            save_last_step_to_exodus(
                    const std::string & aPath,
                    const double aTimeStep = 0.0 );

// -----------------------------------------------------------------------------

            /**
             * save the mesh to an hdf5 file
             */
            void
            save_to_hdf5( const std::string & aPath );


// -----------------------------------------------------------------------------

            /**
             * store the T-Matrices and B-Spline IDs into a file
             */
            void
            save_coeffs_to_binary_files(
                    const std::string & aFilePath );

// -----------------------------------------------------------------------------

            /**
             * store the T-Matrices and B-Spline IDs into a file
             */
            void
            save_coeffs_to_hdf5_file( const std::string & aFilePath );

// -----------------------------------------------------------------------------

            /**
             * loads a field from an HDF5 file and creates a smart pointer
             * to it
             */
            //std::shared_ptr< Field>
            std::shared_ptr< Field >
            load_field_from_hdf5_file(
                    const std::string & aFilePath,
                    const uint          aLagrangeOrder=0,
                    const uint          aBSpineOrder=0 );

// -----------------------------------------------------------------------------

            /**
             * flags active elements
             *
             * @param[ in ]   aElements            element pointers that are to be flagged
             * @param[ in ]   aMinRefinementLevel  if the level of the child is less than this value
             *                                     the child is automatically flagged in the next iteration
             */
            void
            flag_elements(
                    Cell< mtk::Cell* > & aElements,
                    const uint         aMinRefinementLevel = 0 );

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
             * runs the refinement scheme
             */
            void
            perform_refinement();

// -----------------------------------------------------------------------------

            /**
             * copy output pattern to input pattern
             */
            void
            update_refinement_pattern();

// -----------------------------------------------------------------------------

            std::shared_ptr< Field>
            map_field_on_mesh(
                    std::shared_ptr< Field > aField,
                    std::shared_ptr< Mesh >   aMesh );

// -----------------------------------------------------------------------------

            /**
             * Creates an STK interface object.
             * Default: Max Lagrange Order, Outpot pattern
             */
            std::shared_ptr< Mesh >
            create_mesh();

// -----------------------------------------------------------------------------

            /**
             * Creates an STK interface object.
             * Default: Output pattern
             */
            std::shared_ptr< Mesh >
            create_mesh( const uint & aLagrangeOrder );

// -----------------------------------------------------------------------------

            std::shared_ptr< Mesh >
            create_mesh( const uint & aLagrangeOrder, const uint & aPattern );

// -----------------------------------------------------------------------------

            std::shared_ptr< Field >
            create_field( const std::string & aLabel );

// -----------------------------------------------------------------------------


            std::shared_ptr< Field >
            create_field(
                    const std::string & aLabel,
                    const uint        & aLagrangeOrder,
                    const uint        & aBSplineOrder );

// -----------------------------------------------------------------------------

            /**
             * grab the pointer to the datavase
             */
            auto
            get_database() -> decltype ( mDatabase )
            {
                return mDatabase;
            }

// -----------------------------------------------------------------------------

            /**
             * flags elements on the surface and inside of a level set
             */
            uint
            flag_volume_and_surface_elements(
                    const std::shared_ptr<Field> aScalarField );

// -----------------------------------------------------------------------------

            /**
             * flags elements on the surface of a level set
             */
            uint
            flag_surface_elements(
                    const std::shared_ptr<Field> aScalarField);


// -----------------------------------------------------------------------------

            /**
             * Returns elements below the volume and surface refinement criterion.
             * These elements are passed to the Geometry Engine
             */
            void
            get_candidates_for_refinement(
                    Cell< mtk::Cell* > & aCandidates,
                    const uint           aMaxLevel=gMaxNumberOfLevels );

// -----------------------------------------------------------------------------

            /**
             * Flag an element for refinement. Needed for Testing.
             */
            void
            flag_element( const moris_index aElementIndex )
            {
                mDatabase->flag_element( aElementIndex );
            }

// -----------------------------------------------------------------------------

            /**
             * calculate T-Matrices, faces and edges
             */
            void
            finalize();

// -----------------------------------------------------------------------------
// Debug files
// -----------------------------------------------------------------------------

            /**
             * Dumps the background mesh into a VTK file
             *
             * @param[ in ] aFilePath  path of VTK file
             */
            void
            save_background_mesh_to_vtk( const std::string & aFilePath );

// -----------------------------------------------------------------------------

            /**
             * Dumps the Bsplines into a VTK file
             *
             * @param[ in ] aFilePath  path of VTK file
             */
            void
            save_bsplines_to_vtk( const std::string & aFilePath );

// -----------------------------------------------------------------------------

            /**
             * Dumps the faces into a VTK file
             *
             * @param[ in ] aFilePath  path of VTK file
             */
            void
            save_faces_to_vtk( const std::string & aFilePath );

// -----------------------------------------------------------------------------

            /**
             * Dumps the edges into a VTK file
             *
             * @param[ in ] aFilePath  path of VTK file
             */
            void
            save_edges_to_vtk( const std::string & aFilePath );

// -----------------------------------------------------------------------------

            /**
             * Dumps the lagrange mesh into a VTK file
             *
             * @param[ in ] aFilePath  path of VTK file
             */
            void
            save_mesh_to_vtk( const std::string & aFilePath );

// -----------------------------------------------------------------------------

            /**
             * calls the refinement procedure, calculates T-Matrices and
             * performs the mapping
             */
            void
            perform_refinement_and_map_fields();

// -----------------------------------------------------------------------------

            void
            flag_all_active_input_parents();

// -----------------------------------------------------------------------------

            void
            create_input_and_output_meshes();

// -----------------------------------------------------------------------------

        }; /* HMR */

    } /* namespace hmr */
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_HPP_ */
