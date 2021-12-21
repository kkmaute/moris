#ifndef MORIS_FN_PRM_GEN_PARAMETERS_HPP_
#define MORIS_FN_PRM_GEN_PARAMETERS_HPP_

#include "cl_Param_List.hpp"

namespace moris
{
    namespace prm
    {
        // Local to this translation unit, can't be used outside
        namespace
        {
            /**
             * Creates a parameter list that can be used to construct a GEN field (either a geometry or property field).
             *
             * @return Field parameter list
             */
            inline ParameterList
            create_field_parameter_list()
            {
                ParameterList tParameterList;

                tParameterList.insert( "type", "" );// Type of field
                tParameterList.insert( "name", "" );// Name of field
                tParameterList.insert( "field_variable_indices", "" );// Indices of field variables to fill
                tParameterList.insert( "adv_indices", "" );// ADVs used to fill in variables
                tParameterList.insert( "constant_parameters", "" );// Remaining geometry parameters that are constant
                tParameterList.insert( "number_of_refinements", "" );// Number of refinement steps using HMR
                tParameterList.insert( "refinement_mesh_index", "" );// Refinement pattern
                tParameterList.insert( "refinement_function_index", -1 );// Index of user-defined refinement function (-1 = default)
                tParameterList.insert( "discretization_mesh_index", -2 );// Index of B-spline mesh to put this field on (-2 = none, -1 = store)
                tParameterList.insert( "discretization_lower_bound", -1.0 );// Lower bound of level set field (if bspline_mesh_index >= 0)
                tParameterList.insert( "discretization_upper_bound", 1.0 );// Upper bound of level set field (if bspline_mesh_index >= 0)

                return tParameterList;
            }
        }// namespace

        /**
         * Creates a parameter list used for the construction of the geometry engine. One of these parameter lists is
         * always required as the only parameter list in the first cell (index 0).
         *
         * @return GEN parameter list
         */
        inline ParameterList
        create_gen_parameter_list()
        {
            ParameterList tParameterList;

            // Level set parameters
            tParameterList.insert( "intersection_mode", "LEVEL_SET" );// Level-set isocontour level
            tParameterList.insert( "isocontour_threshold", 0.0 );// Level-set isocontour level
            tParameterList.insert( "isocontour_tolerance", 1e-12 );// Interface tolerance based on geometry value
            tParameterList.insert( "intersection_tolerance", 1e-12 );// Interface tolerance based on intersection distance
            tParameterList.insert( "output_mesh_file", "" );// File name for exodus mesh, if default no mesh is written
            tParameterList.insert( "geometry_field_file", "" );// Base file name (without extension) for saving geometry fields
            tParameterList.insert( "time_offset", 0.0 );// Time offset for writing files in optimization process

            // ADVs/IQIs
            tParameterList.insert( "initial_advs", "" );// Initial ADVs, string converted into vector
            tParameterList.insert( "advs_size", 0 );// Specify size and fill value for ADVs in addition to
            tParameterList.insert( "initial_advs_fill", 0.0 );// explicitly defined ADVs (above)
            tParameterList.insert( "lower_bounds", "" );// Lower bounds on advs, string converted into vector
            tParameterList.insert( "lower_bounds_fill", 0.0 );// Fill value for lower bounds up to size of ADV vector
            tParameterList.insert( "upper_bounds", "" );// Upper bounds on advs, string converted into vector
            tParameterList.insert( "upper_bounds_fill", 0.0 );// Fill value for upper bounds up to size of ADV vector
            tParameterList.insert( "IQI_types", "" );// Requested IQI types for sensitivity analysis
            tParameterList.insert( "PDV_types", "" );// Requested PDV types for sensitivity analysis

            // Phase table
            tParameterList.insert( "phase_table", "" );// Construct phase table directly
            tParameterList.insert( "phase_function_name", "" );// User-defined function for determining phase indices
            tParameterList.insert( "number_of_phases", 0 );// Number of phases for a user-defined phase function
            tParameterList.insert( "print_phase_table", false );// Whether or not to print the phase table

            // diagnostics
            tParameterList.insert( "diagnostics", false );
            tParameterList.insert( "diagnostics_id", "" );
            tParameterList.insert( "diagnostics_path", "" );

            return tParameterList;
        }

        /**
         * Creates a parameter list that can be used to construct a geometry field. Any number of these can be added to
         * the second cell (index 1) of then parameter lists for GEN.
         *
         * @return Geometry parameter list
         */
        inline ParameterList
        create_geometry_parameter_list()
        {
            ParameterList tParameterList = create_field_parameter_list();
            tParameterList.insert( "multilinear_intersections", false );

            return tParameterList;
        }

        /**
         * Same as a geometry parameter list, but forces the user to specify two more parameters which name the
         * user-defined functions used to evaluate a geometry field and to evaluate sensitivities.
         *
         * @return User-defined geometry parameter list
         */
        inline ParameterList
        create_user_defined_geometry_parameter_list()
        {
            ParameterList tParameterList = create_geometry_parameter_list();

            tParameterList.set( "type", "user_defined" );// User-defined geometry
            tParameterList.insert( "field_function_name", "" );// Function name for evaluating the geometry field
            tParameterList.insert( "sensitivity_function_name", "" );// Function name for evaluating the sensitivity of the field

            return tParameterList;
        }

        /**
         * Same as a geometry parameter list, but forces the user to specify the path to the voxel file
         *
         * @return Voxel geometrie parameterlist
         */
        inline ParameterList
        create_voxel_field_parameter_list()
        {
            ParameterList tParameterList = create_geometry_parameter_list();

            tParameterList.set( "type", "voxel" );// User-defined geometry
            tParameterList.insert( "voxel_field_file", "" );    // voxel file
            tParameterList.insert( "domain_dimensions", "" );   // domain size
            tParameterList.insert( "domain_offset", "" );

            return tParameterList;
        }

        /**
         * Same as a geometry parameter list, but forces the user to specify the path to the sdf file
         *
         * @return sdf geometrie parameterlist
         */
        inline
        ParameterList create_sdf_field_parameter_list()
        {
            ParameterList tParameterList = create_geometry_parameter_list();

            tParameterList.set("type", "sdf_field");             // User-defined geometry
            tParameterList.insert("sdf_object_path", "");        // obj file
            tParameterList.insert("sdf_object_offset", "");
            tParameterList.insert("sdf_shift", 0.0);

            return tParameterList;
        }

        /**
         * Same as a geometry parameter list, but forces the user to specify the path to the level set file
         *
         * @return sdf geometrie parameterlist
         */
        inline
        ParameterList create_nodal_field_from_file_parameter_list()
        {
            ParameterList tParameterList = create_geometry_parameter_list();

            tParameterList.set("type", "nodal_field_from_file");   // field defined on mesh file
            tParameterList.insert("file_name",   "");              // file name
            tParameterList.insert("field_name",  "");              // field name
            tParameterList.insert("file_format", "exodus");        // file format
            tParameterList.insert("offset", 0.0);                  // offset of field value

            return tParameterList;
        }

        /**
         * Creates a parameter list that will be used to create a swiss cheese slice. The parameters here are different
         * than those of a typical geometry.
         *
         * @return Swiss cheese slice parameter list
         */
        inline ParameterList
        create_swiss_cheese_slice_parameter_list()
        {
            ParameterList tParameterList;

            tParameterList.insert( "type", "swiss_cheese_slice" );// Type of geometry, do not change
            tParameterList.insert( "name", "" );// Name of geometry, can change

            // Must change
            tParameterList.insert( "left_bound", 0.0 );// Left-most hole center
            tParameterList.insert( "right_bound", 0.0 );// Right-most hole center
            tParameterList.insert( "bottom_bound", 0.0 );// Bottom-most hole center
            tParameterList.insert( "top_bound", 0.0 );// Top-most hole center
            tParameterList.insert( "hole_x_semidiameter", 0.0 );// Superellipse semi-diameter in the x direction
            tParameterList.insert( "hole_y_semidiameter", 0.0 );// Superellipse semi-diameter in the y direction

            // One of two options for hole spacing
            tParameterList.insert( "number_of_x_holes", 0 );// Number of holes in the x direction
            tParameterList.insert( "number_of_y_holes", 0 );// Number of holes in the y direction
            tParameterList.insert( "target_x_spacing", 0.0 );// Targeted spacing between hole centers in the x direction
            tParameterList.insert( "target_y_spacing", 0.0 );// Targeted spacing between hole centers in the y direction
            tParameterList.insert( "allow_less_than_target_spacing", true );

            // Optional
            tParameterList.insert( "superellipse_exponent", 2.0 );// Superellipse exponent
            tParameterList.insert( "superellipse_scaling", 1.0 );// Superellipse scaling
            tParameterList.insert( "superellipse_regularization", 1e-8 );// Superellipse regularization
            tParameterList.insert( "superellipse_shift", 1e-6 );// Superellipse shift

            tParameterList.insert( "row_offset", 0.0 );// Offset to be applied on subsequent rows
            tParameterList.insert( "number_of_refinements", "" );// Number of refinement steps using HMR
            tParameterList.insert( "refinement_mesh_index", "" );// Refinement pattern
            tParameterList.insert( "refinement_function_index", -1 );// Index of user-defined refinement function (-1 = none)
            tParameterList.insert( "discretization_mesh_index", -1 );// Index of B-spline mesh to create level set field on (-1 = none)
            tParameterList.insert( "discretization_lower_bound", -1.0 );// Lower bound of level set field (if bspline_mesh_index >= 0)
            tParameterList.insert( "discretization_upper_bound", 1.0 );// Upper bound of level set field (if bspline_mesh_index >= 0)
            tParameterList.insert( "multilinear_intersections", false );

            return tParameterList;
        }

        /**
         * Creates a parameter list that can be used to construct a property field. Any number of these can be added to
         * the third cell (index 2) of the parameter lists for GEN.
         *
         * @return GEN property parameter list
         */
        inline ParameterList
        create_gen_property_parameter_list()
        {
            ParameterList tParameterList = create_field_parameter_list();

            tParameterList.insert( "dependencies", "" );// Names of other fields that this property depends on
            tParameterList.insert( "pdv_type", "" );// The type of PDV that this property will be assigned to
            tParameterList.insert( "pdv_mesh_type", "interpolation" );// Mesh type for assigning PDVs
            tParameterList.insert( "pdv_mesh_set_names", "" );// Mesh set names for assigning PDVs
            tParameterList.insert( "pdv_mesh_set_indices", "" );// Mesh set indices for assigning PDVs

            return tParameterList;
        }

        /**
         * Same as a property parameter list, but forces the user to specify two more parameters which name the
         * user-defined functions used to evaluate a property field and to evaluate sensitivities.
         *
         * @return User-defined GEN property parameter list
         */
        inline ParameterList
        create_user_defined_property_parameter_list()
        {
            ParameterList tParameterList = create_gen_property_parameter_list();

            tParameterList.set( "type", "user_defined" );// User-defined property
            tParameterList.insert( "field_function_name", "" );// Function name for evaluating the property field
            tParameterList.insert( "sensitivity_function_name", "" );// Function name for evaluating the sensitivity of the field

            return tParameterList;
        }

    }// namespace prm
}// namespace moris

#endif /* MORIS_FN_PRM_GEN_PARAMETERS_HPP_ */
