/*
 * cl_HMR.hpp
 *
 *  Created on: May 5, 2018
 *      Author: messe
 */

#ifndef SRC_HMR_CL_HMR_HPP_
#define SRC_HMR_CL_HMR_HPP_


#include "cl_HMR_Database.hpp"     //HMR/src
#include "cl_HMR_Element.hpp"
#include "cl_HMR_Field_Param.hpp"
#include "cl_HMR_Mesh.hpp"
#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR_Mesh_Integration.hpp"
#include "cl_HMR_Parameters.hpp"     //HMR/src
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
            bool mPerformRefinementCalled = false;
            bool mUpdateRefinementCalled  = false;

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
             * alternative constructor which loads input and output patterns from path
             */
            HMR( const std::string & aInPath, const std::string & aOutPath );

// -----------------------------------------------------------------------------

            /**
             * default destructor of HMR
             */
            ~HMR ( ){};

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
                    const double aTimeStep = 0.0,
                    const uint aOutputOrder = 0 );

// -----------------------------------------------------------------------------

            /**
             * renumber nodes for femdoc and save the mesh to an exodus file. HACK with consent of Kurt
             */
            void
            renumber_and_save_to_exodus(
                    const std::string & aPath,
                    const double aTimeStep = 0.0,
                    const uint aOutputOrder = 0 );

// -----------------------------------------------------------------------------

            /**
             * save the mesh to an exodus file
             */
            void
            save_last_step_to_exodus(
                    const std::string & aPath,
                    const double aTimeStep = 0.0,
                    const uint aOutputOrder = 0 );

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
            save_coeffs_to_hdf5_file( const std::string & aFilePath );

            // -----------------------------------------------------------------------------

            void
            renumber_and_save_coeffs_to_hdf5_file_HACK( const std::string & aFilePath );

// -----------------------------------------------------------------------------

            /**
             * store the T-Matrices and B-Spline IDs into a file
             */
            void
            save_mesh_relations_to_hdf5_file ( const std::string & aFilePath );

// -----------------------------------------------------------------------------

            /**
             * loads a field from an HDF5 file and creates a smart pointer
             * to it
             */
            std::shared_ptr< Field >
            load_field_from_hdf5_file(
                    const std::string & aLabel,
                    const std::string & aFilePath,
                    const uint          aLagrangeOrder=0,
                    const uint          aBSpineOrder=0 );

// -----------------------------------------------------------------------------

            std::shared_ptr< Field >
            load_field_from_exo_file(
                    const std::string & aLabel,
                    const std::string & aFilePath,
                    const uint          aLagrangeOrder=0,
                    const uint          aBSpineOrder=0 );

// -----------------------------------------------------------------------------

            std::shared_ptr< Field >
            load_field_from_file(
                    const std::string & aLabel,
                    const std::string & aFilePath,
                    const uint          aLagrangeOrder=0,
                    const uint          aBSpineOrder=0 );

// -----------------------------------------------------------------------------

            /**
             * flags active elements
             *
             * @param[ in ]   aElements            element pointers that are to be flagged
             *  @param[ in ]  aLagrangeSwitch      off: flag b-spline elements, on: flag lagrange elements
             * @param[ in ]   aMinRefinementLevel  if the level of the child is less than this value
             *                                     the child is automatically flagged in the next iteration
             */
            void flag_elements(       Cell< mtk::Cell* > & aElements,
                                const uint                 aMinRefinementLevel = 0 );

// -----------------------------------------------------------------------------

            /**
             * exposes the parameters pointer
             */
            Parameters * get_parameters()
            {
                return mParameters;
            }

// -----------------------------------------------------------------------------

            /**
             * runs the refinement scheme
             */
            void perform_refinement( const enum RefinementMode aRefinementMode );

// -----------------------------------------------------------------------------

            /**
             * copy output pattern to input pattern
             */
            void update_refinement_pattern();

// -----------------------------------------------------------------------------

            /**
             * Creates an STK interface object.
             * Default: Max Lagrange Order, Outpot pattern
             */
            std::shared_ptr< Mesh > create_mesh();

// -----------------------------------------------------------------------------

            /**
             * Creates an STK interface object.
             * Default: Output pattern
             */
            std::shared_ptr< Mesh > create_mesh( const uint & aLagrangeOrder );

// -----------------------------------------------------------------------------

            std::shared_ptr< Mesh > create_mesh( const uint & aLagrangeOrder, const uint & aPattern );

// -----------------------------------------------------------------------------

            std::shared_ptr< Interpolation_Mesh_HMR >
            create_interpolation_mesh( const uint & aLagrangeOrder, const uint & aPattern );

// -----------------------------------------------------------------------------

            std::shared_ptr< Integration_Mesh_HMR >
            create_integration_mesh( const uint &             aLagrangeOrder,
                                     const uint &             aPattern,
                                     Interpolation_Mesh_HMR & aInterpolationMesh);

// -----------------------------------------------------------------------------
            std::shared_ptr< Field >
            create_field( const std::string & aLabel );

// -----------------------------------------------------------------------------

            std::shared_ptr< Field >
            create_field( const std::string & aLabel,
                          const uint        & aLagrangeOrder,
                          const uint        & aBSplineOrder );

 // -----------------------------------------------------------------------------

            // create field from parameter list
            std::shared_ptr< Field >
            create_field( const Field_Param & aParameters );

// -----------------------------------------------------------------------------

            /**
             * grab the pointer to the datavase
             */
            // std::shared_ptr< Database >
            auto
            get_database() -> decltype ( mDatabase )
            {
                return mDatabase;
            }

// -----------------------------------------------------------------------------

            /**
             * Flag an element for refinement. Needed for Testing.
             */
            void flag_element( const moris_index aElementIndex )
            {
                mDatabase->flag_element( aElementIndex );
            }

// -----------------------------------------------------------------------------

            /**
             * for flagging
             */
            void
            get_candidates_for_refinement(
                    Cell< mtk::Cell* > & aCandidates,
                    const uint           aMaxLevel=gMaxNumberOfLevels );


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
                    const std::shared_ptr<Field> aScalarField );

// -----------------------------------------------------------------------------

            /**
             * special function for tutorial
             */
            void
            perform_refinement_and_map_fields();

// -----------------------------------------------------------------------------

            /**
             * funciton for L2 test
             */
            void
            map_field_to_output( std::shared_ptr< Field > aField );

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

            void
            create_input_and_output_meshes();

// -----------------------------------------------------------------------------

            void
            perform_initial_refinement();

// -----------------------------------------------------------------------------

            void
            user_defined_flagging(
                    int (*aFunction)(
                                 Element                    * aElement,
                            const Cell< Matrix< DDRMat > >   & aElementLocalValues,
                                  ParameterList              & aParameters ),
                            Cell< std::shared_ptr< Field > > & aFields,
                                  ParameterList              & aParameters );

// -----------------------------------------------------------------------------

            uint
            get_mesh_index( const uint aOrder, const uint aPattern );

// -----------------------------------------------------------------------------
        }; /* HMR */

    } /* namespace hmr */
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_HPP_ */
