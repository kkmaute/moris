/*
 * cl_HMR.hpp
 *
 *  Created on: May 5, 2018
 *      Author: schmidt
 */

#ifndef SRC_HMR_CL_HMR_HPP_
#define SRC_HMR_CL_HMR_HPP_

//#include "cl_HMR_Database.hpp"     //HMR/src
#include "cl_HMR_Element.hpp"
#include "cl_HMR_Field_Param.hpp"
#include "cl_HMR_Parameters.hpp"     //HMR/src
#include "cl_Library_IO.hpp"

namespace moris
{
    namespace mtk
    {
        class Mesh_Manager;
    }
    namespace hmr
    {
        class Field;
        class Database;
        class Mesh;
        class Lagrange_Mesh_Base;
        class Interpolation_Mesh_HMR;
        class Integration_Mesh_HMR;
        // -----------------------------------------------------------------------------
        /**
         * \brief the main class of HMR
         */
        class HMR
        {
            public :
                //! object containing user settings
                Parameters * mParameters;

                std::shared_ptr< Database > mDatabase;

                //! flag telling if perform_refinement() has been called
                bool mPerformRefinementCalled = false;
                bool mUpdateRefinementCalled  = false;

                //! mesh which points to input pattern
                Cell< std::shared_ptr< Mesh > > mMeshes;
                Cell< std::shared_ptr< Mesh > > mInputMeshes;

                //! mesh which points to output pattern
                //Cell< std::shared_ptr< Mesh > > mOutputMeshes;

                //! container with field objects
                Cell< std::shared_ptr< Field > > mFields;

                //! map for Lagrange orders
                Matrix< DDUMat > mLagrangeOrderToInputMeshIndexMap;

                std::shared_ptr< mtk::Mesh_Manager > mMTKPerformer = nullptr;

                /*
                 * @brief determines elements (cells) intersected by the level set
                 *
                 * @param[in] aCells        - elements to be flagged for refinement
                 * @param[in] aCandidates   - candidates for refinement
                 * @param[in] aVertexValues - vertex values of scalar field
                 * @param[in] aLowerBound   - lower bound of LS
                 * @param[in] aUpperBound   - upper bound of LS
                 */
                void find_cells_intersected_by_levelset(
                        Cell< hmr::Element * >   & aCells,
                        Cell< hmr::Element * >   & aCandidates,
                        const  Matrix< DDRMat >  & aVertexValues,
                        const  real                aLowerBound = -0.0001,
                        const  real                aUpperBound =  0.0001);

                void find_low_level_cells_intersected_by_levelset(
                        Cell< hmr::Element * >   & aCells,
                        Cell< hmr::Element * >   & aCandidates,
                        const  Matrix< DDRMat >  & aVertexValues,
                        const  real                aLowerBound = -0.0001,
                        const  real                aUpperBound =  0.0001);

                /*
                 * @brief determines volume elements (cells)
                 *
                 * @param[in] aCells        - elements to be flagged for refinement
                 * @param[in] aCandidates   - candidates for refinement
                 * @param[in] aVertexValues - vertex values of scalar field
                 * @param[in] aUpperBound   - upper bound of LS
                 */
                void find_cells_within_levelset(
                        Cell< hmr::Element * >  & aCells,
                        Cell< hmr::Element * >  & aCandidates,
                        const  Matrix< DDRMat > & aVertexValues,
                        const  uint               aUpperBound = 0.0 );


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
                HMR ( 
                    ParameterList                       & aParameterList,
                    std::shared_ptr<moris::Library_IO>    aLibrary = nullptr ) ;

                // -----------------------------------------------------------------------------

                /**
                 * alternative constructor which loads a mesh from a h5 file
                 */
                HMR( const std::string & aPath );

                // -----------------------------------------------------------------------------

                /**
                 * alternative constructor which loads input and output patterns from path
                 */
                HMR(
                        const std::string & aInPath,
                        const std::string & aOutPath );

                // -----------------------------------------------------------------------------

                /**
                 * default destructor of HMR
                 */
                ~HMR ();

                // -----------------------------------------------------------------------------
                /**
                 * save the mesh to an exodus file
                 */
                void save_to_exodus(
                        const uint        & aMeshIndex,
                        const std::string & aPath,
                        const double        aTimeStep = 0.0 );

                // -----------------------------------------------------------------------------

                /**
                 * save the mesh to an exodus file
                 */
                void save_last_step_to_exodus(
                        const uint          aIndex,
                        const std::string & aPath,
                        const double        aTimeStep = 0.0 );

                // -----------------------------------------------------------------------------

                /**
                 * save the mesh to an hdf5 file
                 */
                void save_to_hdf5(
                        const std::string & aPath,
                        const uint          aLagrangeMeshIndex );

                // -----------------------------------------------------------------------------

                /**
                 * store the T-Matrices and B-Spline IDs into a file
                 */
                void save_coeffs_to_hdf5_file(
                        const std::string & aFilePath,
                        const uint        & aLagrangeMeshIndex );

                // -----------------------------------------------------------------------------

                /**
                 * store the T-Matrices and B-Spline IDs into a file
                 */
                void save_mesh_relations_to_hdf5_file (
                        const std::string & aFilePath,
                        const uint        & aLagrangeMeshIndex,
                        const uint        & aBsplineMeshIndex );

                // -----------------------------------------------------------------------------

                /**
                 * loads a field from an HDF5 file and creates a smart pointer
                 * to it
                 */
                std::shared_ptr< Field > load_field_from_hdf5_file(
                        const std::string & aLabel,
                        const std::string & aFilePath,
                        const uint          aLagrangeIndex = 0,
                        const uint          aBSpineIndex = 0 );

                // -----------------------------------------------------------------------------

                std::shared_ptr< Field > load_field_from_exo_file(
                        const std::string & aLabel,
                        const std::string & aFilePath,
                        const uint          aLagrangeIndex=0,
                        const uint          aBSpineIndex=0 );

                // -----------------------------------------------------------------------------

                std::shared_ptr< Field > load_field_from_file(
                        const std::string & aLabel,
                        const std::string & aFilePath,
                        const uint          aLagrangeIndex=0,
                        const uint          aBSpineIndex=0 );

                // -----------------------------------------------------------------------------

                // -----------------------------------------------------------------------------

                /**
                 * flags active elements on working pattern
                 *
                 * @param[ in ]   aElements            element pointers that are to be flagged
                 * @param[ in ]   aMinRefinementLevel  if the level of the child is less than this value
                 *                                     the child is automatically flagged in the next iteration
                 */
                void flag_elements_on_working_pattern(
                        Cell< hmr::Element* > & aElements,
                        const uint               aMinRefinementLevel = 0 );

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
                void perform_refinement_based_on_working_pattern(
                        const uint aPattern,
                        const bool aResetPattern = false );

                void perform_refinement( const uint aPattern );

                // -----------------------------------------------------------------------------

                void put_elements_on_refinment_queue( Cell< hmr::Element* > & aElements);

                // -----------------------------------------------------------------------------

                /**
                 * copy output pattern to input pattern
                 */
                void update_refinement_pattern( const uint aPattern) ;

                // -----------------------------------------------------------------------------

                /**
                 * set activation pattern
                 */
                void set_activation_pattern( const uint & aActivationPattern );

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

                std::shared_ptr< Mesh > create_mesh( const uint & aLagrangeOrder,
                        const uint & aPattern );

                std::shared_ptr< Mesh > create_mesh(
                        const uint & aLagrangeOrder,
                        const uint & aLagrangePattern,
                        const uint & aBsplinePattern );

                // -----------------------------------------------------------------------------

                Interpolation_Mesh_HMR * create_interpolation_mesh( const uint & aLagrangeMeshIndex );

                // -----------------------------------------------------------------------------

                Interpolation_Mesh_HMR * create_interpolation_mesh(
                        const std::string & aName,
                        bool                aTMatricesExist = false );

                // -----------------------------------------------------------------------------

                Interpolation_Mesh_HMR * create_interpolation_mesh(
                        const uint & aLagrangeOrder,
                        const uint & aPattern );

                Interpolation_Mesh_HMR * create_interpolation_mesh(
                        const uint & aOrder,
                        const uint & aLagrangePattern,
                        const uint & aBsplinePattern);

                // -----------------------------------------------------------------------------

                Integration_Mesh_HMR * create_integration_mesh(
                        const uint             & aLagrangeOrder,
                        const uint             & aPattern,
                        Interpolation_Mesh_HMR * aInterpolationMesh);

                // -----------------------------------------------------------------------------

                Integration_Mesh_HMR * create_integration_mesh(
                        const uint             & aLagrangeMeshIndex,
                        Interpolation_Mesh_HMR * aInterpolationMesh);

                // -----------------------------------------------------------------------------
                std::shared_ptr< Field > create_field( const std::string & aLabel );

                // -----------------------------------------------------------------------------

                std::shared_ptr< Field > create_field(
                        const std::string & aLabel,
                        const uint        & aLagrangeOrder,
                        const uint        & aBSplineOrder );

                // -----------------------------------------------------------------------------

                // create field from parameter list
                std::shared_ptr< Field > create_field( const Field_Param & aParameters );

                // -----------------------------------------------------------------------------

                /**
                 * grab the pointer to the datavase
                 */
                // std::shared_ptr< Database >
                auto get_database() -> decltype ( mDatabase )
                {
                    return mDatabase;
                 }

                // -----------------------------------------------------------------------------

                /**
                 * Flag an element for refinement. Needed for Testing.
                 */
                void flag_element( const moris_index aElementIndex );

                // -----------------------------------------------------------------------------

                /**
                 * for flagging
                 */
                void get_candidates_for_refinement(
                        Cell< hmr::Element* > & aCandidates,
                        const uint              aLagrangeMeshIndex);

                void get_candidates_for_refinement(
                        Cell< hmr::Element* > & aCandidates,
                        Lagrange_Mesh_Base    * aMesh );

                void get_active_candidates_for_refinement(
                        Cell< hmr::Element* > & aCandidates,
                        Lagrange_Mesh_Base    * aMesh );


                // -----------------------------------------------------------------------------

                /**
                 * flags elements on the surface and inside of a level set
                 */
                uint flag_volume_and_surface_elements_on_working_pattern( const std::shared_ptr<Field> aScalarField );

                // -----------------------------------------------------------------------------

                /**
                 * flags elements on the surface of a level set
                 */
                uint flag_surface_elements_on_working_pattern( const std::shared_ptr<Field> aScalarField );

                uint based_on_field_put_elements_on_queue(
                        const Matrix< DDRMat > & aFieldValues,
                        uint                     aLagrangeMeshIndex,
                        sint                     aFunctionIndex = -1);

                uint based_on_field_put_elements_on_queue(
                        const Matrix< DDRMat > & aFieldValues,
                        uint                     aPattern,
                        uint                     aOrder,
                        sint                     aFunctionIndex);

                uint based_on_field_put_low_level_elements_on_queue(
                        const Matrix< DDRMat > & aFieldValues,
                        uint                     aLagrangeMeshIndex,
                        sint                     aFunctionIndex);

                uint based_on_field_put_low_level_elements_on_queue(
                        const Matrix< DDRMat > & aFieldValues,
                        uint                     aPattern,
                        uint                     aOrder,
                        sint                     aFunctionIndex);

                // -----------------------------------------------------------------------------

                /**
                 * special function for tutorial
                 */
                //            void perform_refinement_and_map_fields( const uint aPattern );

                // -----------------------------------------------------------------------------

                /**
                 * funciton for L2 test
                 */
                void map_field_to_output(
                        std::shared_ptr< Field > aField,
                        const uint               aMesh_Index,
                        const uint               aBsplineMeshIndex);

                void map_field_to_output_union(
                        std::shared_ptr< Field > aField,
                        const uint               aUnionOrder );

                // -----------------------------------------------------------------------------

                /**
                 * calculate T-Matrices, faces and edges
                 */
                void finalize();

                // -----------------------------------------------------------------------------

                /**
                 * resets HMR. deletes T-Matrices, faces and edges
                 */

                void reset_HMR();

                // -----------------------------------------------------------------------------

                /**
                 * performs HMR specific task. Right now only finalize. more will be added
                 */
                void perform();

                // -----------------------------------------------------------------------------

                void set_performer( std::shared_ptr< mtk::Mesh_Manager > aMTKPerformer );

                // -----------------------------------------------------------------------------

                bool get_mesh_name_exists( const std::string & aName ) const;

                // -----------------------------------------------------------------------------
                // Debug files
                // -----------------------------------------------------------------------------

                /**
                 * Dumps the background mesh into a VTK file
                 *
                 * @param[ in ] aFilePath  path of VTK file
                 */
                void save_background_mesh_to_vtk( const std::string & aFilePath );

                // -----------------------------------------------------------------------------

                /**
                 * Dumps the Bsplines into a VTK file
                 *
                 * @param[ in ] aFilePath  path of VTK file
                 */
                void save_bsplines_to_vtk( const std::string & aFilePath,
                        const uint        & aLagrangeMeshIndex,
                        const uint        & aBsplineMeshIndex  );

                // -----------------------------------------------------------------------------

                void calculate_bspline_coordinates( const uint        & aLagrangeMeshIndex,
                        const uint        & aBsplineMeshIndex  );

                // -----------------------------------------------------------------------------

                /**
                 * Dumps the faces into a VTK file
                 *
                 * @param[ in ] aFilePath  path of VTK file
                 */
                void save_faces_to_vtk( const std::string & aFilePath,
                        const uint        & aLagrangeMeshIndex );

                // -----------------------------------------------------------------------------

                /**
                 * Dumps the edges into a VTK file
                 *
                 * @param[ in ] aFilePath  path of VTK file
                 */
                void save_edges_to_vtk( const std::string & aFilePath,
                        const uint        & aLagrangeMeshIndex );

                // -----------------------------------------------------------------------------

                /**
                 * Dumps the lagrange mesh into a VTK file
                 *
                 * @param[ in ] aFilePath  path of VTK file
                 */
                void save_mesh_to_vtk( const std::string & aFilePath,
                        const uint        & aLagrangeMeshIndex );

                // -----------------------------------------------------------------------------

                void create_input_and_output_meshes();

                // -----------------------------------------------------------------------------

                void perform_initial_refinement();

                // -----------------------------------------------------------------------------

        private:

            void user_defined_flagging(
                    Cell< hmr::Element * >   & aCells,
                    Cell< hmr::Element * >   & aCandidates,
                    const  Matrix< DDRMat >  & aVertexValues,
                    uint                       aFunctionIndex);

        }; /* HMR */

    } /* namespace hmr */
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_HPP_ */
