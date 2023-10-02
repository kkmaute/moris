/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_PRM_GEN_Parameters.hpp
 *
 */

#pragma once

#include "cl_Param_List.hpp"

namespace moris::prm
{
    /**
     * Inserts all general field parameters to the given parameter list.
     *
     * @param aDesignParameterList Parameter list for a design (level set geometry or property)
     */
    static void insert_field_parameters( ParameterList& aDesignParameterList )
    {
        aDesignParameterList.insert( "field_type", "" );                   // Type of field
        aDesignParameterList.insert( "name", "" );                         // Name of field
        aDesignParameterList.insert( "dependencies", "" );                 // Names of other fields that this field depends on
        aDesignParameterList.insert( "field_variable_indices", "" );       // Indices of field variables to fill
        aDesignParameterList.insert( "adv_indices", "" );                  // ADVs used to fill in variables
        aDesignParameterList.insert( "constant_parameters", "" );          // Remaining geometry parameters that are constant
        aDesignParameterList.insert( "number_of_refinements", "" );        // Number of refinement steps using HMR
        aDesignParameterList.insert( "refinement_mesh_index", "" );        // Refinement pattern
        aDesignParameterList.insert( "refinement_function_index", -1 );    // Index of user-defined refinement function (-1 = default)
        aDesignParameterList.insert( "discretization_mesh_index", -2 );    // Index of B-spline mesh to put this field on (-2 = none, -1 = store)
        aDesignParameterList.insert( "discretization_lower_bound", -1.0 ); // Lower bound of level set field (if bspline_mesh_index >= 0)
        aDesignParameterList.insert( "discretization_upper_bound", 1.0 );  // Upper bound of level set field (if bspline_mesh_index >= 0)
    }

    /**
     * Inserts parameters to a field parameter list useful for getting user-defined functions from an input file.
     *
     * @param aDesignParameterList Parameter list for a design field with field parameters already inserted
     */
    static void insert_user_defined_field_parameters( ParameterList& aDesignParameterList )
    {
        aDesignParameterList.set( "field_type", "user_defined" );          // User-defined geometry
        aDesignParameterList.insert( "field_function_name", "" );          // Function name for evaluating the geometry field
        aDesignParameterList.insert( "sensitivity_function_name", "" );    // Function name for evaluating the sensitivity of the field
    }
    
    /**
     * Creates a parameter list for the construction of a geometry.
     * 
     * @return Geometry parameter list
     */
    static ParameterList create_geometry_parameter_list()
    {
        ParameterList tParameterList;
        tParameterList.insert( "design_type", "geometry" ); // Set the design type to a geometry

        return tParameterList;
    }

    //------------------------------------------------------------------------------

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
        tParameterList.insert( "evaluate_new_pts_as_linear", false );    // for the new vertices, should I use full background cell info or only the linear version
        tParameterList.insert( "output_mesh_file", "" );                 // File name for exodus mesh, if default no mesh is written
        tParameterList.insert( "geometry_field_file", "" );              // Base file name (without extension) for saving geometry fields
        tParameterList.insert( "time_offset", 0.0 );                     // Time offset for writing files in optimization process

        // ADVs/IQIs
        tParameterList.insert( "initial_advs", "" );          // Initial ADVs, string converted into vector
        tParameterList.insert( "advs_size", 0 );              // Specify size and fill value for ADVs in addition to
        tParameterList.insert( "initial_advs_fill", 0.0 );    // explicitly defined ADVs (above)
        tParameterList.insert( "lower_bounds", "" );          // Lower bounds on advs, string converted into vector
        tParameterList.insert( "lower_bounds_fill", 0.0 );    // Fill value for lower bounds up to size of ADV vector
        tParameterList.insert( "upper_bounds", "" );          // Upper bounds on advs, string converted into vector
        tParameterList.insert( "upper_bounds_fill", 0.0 );    // Fill value for upper bounds up to size of ADV vector
        tParameterList.insert( "IQI_types", "" );             // Requested IQI types for sensitivity analysis
        tParameterList.insert( "PDV_types", "" );             // Requested PDV types for sensitivity analysis

        // Phase table
        tParameterList.insert( "phase_table", "" );             // Construct phase table directly
        tParameterList.insert( "phase_function_name", "" );     // User-defined function for determining phase indices
        tParameterList.insert( "number_of_phases", 0 );         // Number of phases for a user-defined phase function
        tParameterList.insert( "print_phase_table", false );    // Whether to print the phase table

        // diagnostics
        tParameterList.insert( "diagnostics", false );
        tParameterList.insert( "diagnostics_id", "" );
        tParameterList.insert( "diagnostics_path", "" );

        return tParameterList;
    }

    //------------------------------------------------------------------------------

    /**
     * Creates a parameter list that can be used to construct a geometry field. Any number of these can be added to
     * the second cell (index 1) of then parameter lists for GEN.
     *
     * @return Geometry parameter list
     */
    inline ParameterList
    create_level_set_geometry_parameter_list()
    {
        ParameterList tParameterList = create_geometry_parameter_list(); // Inserts all geometry parameters
        tParameterList.insert( "geometry_type", "level_set" );           // Set the geometry type to level set
        insert_field_parameters( tParameterList );                       // Inserts all field parameters
        tParameterList.insert( "multilinear_intersections", false );     // Whether to use multilinear interoplation for calculating intersections
        tParameterList.insert( "intersection_mode", "LEVEL_SET" );       // Deprecated
        tParameterList.insert( "isocontour_threshold", 0.0 );            // Level set isocontour level
        tParameterList.insert( "isocontour_tolerance", 1e-12 );          // Interface tolerance based on geometry value
        tParameterList.insert( "intersection_tolerance", 1e-12 );        // Interface tolerance based on intersection distance

        return tParameterList;
    }

    //------------------------------------------------------------------------------

    /**
     * Same as a geometry parameter list, but forces the user to specify two more parameters which name the
     * user-defined functions used to evaluate a geometry field and to evaluate sensitivities.
     *
     * @return User-defined geometry parameter list
     */
    inline ParameterList
    create_user_defined_geometry_parameter_list()
    {
        ParameterList tParameterList = create_level_set_geometry_parameter_list(); // Level set geometry parameters
        insert_user_defined_field_parameters( tParameterList );                    // Parameters for reading user-defined field functions

        return tParameterList;
    }

    //------------------------------------------------------------------------------

    /**
     * Same as a geometry parameter list, but forces the user to specify the path to the voxel file
     *
     * @return Voxel geometry parameterlist
     */
    inline ParameterList
    create_voxel_field_parameter_list()
    {
        ParameterList tParameterList = create_level_set_geometry_parameter_list();

        tParameterList.set( "field_type", "voxel" );                // User-defined geometry
        tParameterList.insert( "voxel_field_file", "" );      // voxel file
        tParameterList.insert( "domain_dimensions", "" );     // domain size
        tParameterList.insert( "domain_offset", "" );         // domain offset
        tParameterList.insert( "grain_id_value_map", "" );    // grain id to value map

        return tParameterList;
    }

    //------------------------------------------------------------------------------

    /**
     * Same as a geometry parameter list, but forces the user to specify the path to the sdf file
     *
     * @return sdf geometry parameterlist
     */
    inline ParameterList
    create_sdf_field_parameter_list()
    {
        ParameterList tParameterList = create_level_set_geometry_parameter_list();

        tParameterList.set( "field_type", "sdf_field" );           // User-defined geometry
        tParameterList.insert( "sdf_object_path", "" );      // obj file
        tParameterList.insert( "sdf_object_offset", "" );    // offset of object
        tParameterList.insert( "sdf_shift", 0.0 );           // sdf shift

        return tParameterList;
    }

    //------------------------------------------------------------------------------

    /**
     * Same as a geometry parameter list, but forces the user to specify the path to the image-based sdf file
     *
     * @return sdf geometry parameterlist
     */
    inline ParameterList
    create_image_sdf_field_parameter_list()
    {
        ParameterList tParameterList = create_level_set_geometry_parameter_list();

        tParameterList.set( "field_type", "image_sdf" );                  // sdf field generated from image
        tParameterList.insert( "image_file", "" );                  // image file (hdf5 format)
        tParameterList.insert( "image_dimensions", "" );            // domain size
        tParameterList.insert( "image_offset", "" );                // domain offset
        tParameterList.insert( "image_sdf_scaling", 0.0 );          // sdf scaling factor (0: automatic scaling)
        tParameterList.insert( "image_sdf_shift", 0.0 );            // sdf shift value
        tParameterList.insert( "image_sdf_default", -1.0 );         // sdf value outside image
        tParameterList.insert( "image_sdf_interpolate", false );    // whether sdf value is interpolated

        return tParameterList;
    }

    //------------------------------------------------------------------------------

    /**
     * Same as a geometry parameter list, but forces the user to specify the path to the level set file
     *
     * @return sdf geometry parameterlist
     */
    inline ParameterList
    create_nodal_field_from_file_parameter_list()
    {
        ParameterList tParameterList = create_level_set_geometry_parameter_list();

        tParameterList.set( "field_type", "nodal_field_from_file" );    // field defined on mesh file
        tParameterList.insert( "file_name", "" );                 // file name
        tParameterList.insert( "field_name", "" );                // field name
        tParameterList.insert( "file_format", "exodus" );         // file format
        tParameterList.insert( "offset", 0.0 );                   // offset of field value

        return tParameterList;
    }

    //------------------------------------------------------------------------------

    /**
     * Creates a parameter list that will be used to create a swiss cheese slice. The parameters here are different
     * than those of a typical geometry.
     *
     * @return Swiss cheese slice parameter list
     */
    inline ParameterList
    create_swiss_cheese_slice_parameter_list()
    {
        ParameterList tParameterList = create_level_set_geometry_parameter_list();

        // Sets the type of the geometry
        tParameterList.set( "field_type", "swiss_cheese_slice" );    // Type of geometry, do not change

        // Must change
        tParameterList.insert( "left_bound", 0.0 );             // Left-most hole center
        tParameterList.insert( "right_bound", 0.0 );            // Right-most hole center
        tParameterList.insert( "bottom_bound", 0.0 );           // Bottom-most hole center
        tParameterList.insert( "top_bound", 0.0 );              // Top-most hole center
        tParameterList.insert( "hole_x_semidiameter", 0.0 );    // Superellipse semi-diameter in the x direction
        tParameterList.insert( "hole_y_semidiameter", 0.0 );    // Superellipse semi-diameter in the y direction

        // One of two options for hole spacing
        tParameterList.insert( "number_of_x_holes", 0 );     // Number of holes in the x direction
        tParameterList.insert( "number_of_y_holes", 0 );     // Number of holes in the y direction
        tParameterList.insert( "target_x_spacing", 0.0 );    // Targeted spacing between hole centers in the x direction
        tParameterList.insert( "target_y_spacing", 0.0 );    // Targeted spacing between hole centers in the y direction
        tParameterList.insert( "allow_less_than_target_spacing", true );

        // Optional
        tParameterList.insert( "superellipse_exponent", 2.0 );           // Superellipse exponent
        tParameterList.insert( "superellipse_scaling", 1.0 );            // Superellipse scaling
        tParameterList.insert( "superellipse_regularization", 1e-8 );    // Superellipse regularization
        tParameterList.insert( "superellipse_shift", 1e-6 );             // Superellipse shift
        tParameterList.insert( "row_offset", 0.0 );                     // Offset to be applied on subsequent rows
        tParameterList.set( "discretization_mesh_index", -1 );       // Index of B-spline mesh to create level set field on (-1 = none)

        return tParameterList;
    }

    //------------------------------------------------------------------------------

    /**
     * Creates a parameter list that can be used to construct a property field. Any number of these can be added to
     * the third cell (index 2) of the parameter lists for GEN.
     *
     * @return GEN property parameter list
     */
    inline ParameterList
    create_gen_property_parameter_list()
    {
        ParameterList tParameterList;                                 // For right now, all properties are fields
        tParameterList.insert( "design_type", "property" );           // Set the design type to a property
        insert_field_parameters( tParameterList );                    // Inserts all field parameters
        tParameterList.insert( "pdv_type", "" );                      // The type of PDV that this property will be assigned to
        tParameterList.insert( "pdv_mesh_type", "interpolation" );    // Mesh type for assigning PDVs
        tParameterList.insert( "pdv_mesh_set_names", "" );            // Mesh set names for assigning PDVs
        tParameterList.insert( "pdv_mesh_set_indices", "" );          // Mesh set indices for assigning PDVs

        return tParameterList;
    }

    //------------------------------------------------------------------------------

    /**
     * Same as a property parameter list, but forces the user to specify two more parameters which name the
     * user-defined functions used to evaluate a property field and to evaluate sensitivities.
     *
     * @return User-defined GEN property parameter list
     */
    inline ParameterList
    create_user_defined_property_parameter_list()
    {
        ParameterList tParameterList = create_gen_property_parameter_list(); // Create property parameter list
        insert_user_defined_field_parameters( tParameterList );              // Parameters for reading user-defined field functions

        return tParameterList;
    }

    //------------------------------------------------------------------------------

}
