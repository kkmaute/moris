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
     * Creates a parameter list used for the construction of the geometry engine. One of these parameter lists is
     * always required as the only parameter list in the first cell (index 0).
     *
     * @return GEN parameter list
     */
    inline ParameterList
    create_gen_parameter_list()
    {
        ParameterList tGENParameterList;

        // Level set parameters
        tGENParameterList.insert( "evaluate_new_pts_as_linear", false );    // for the new vertices, should I use full background cell info or only the linear version
        tGENParameterList.insert( "output_mesh_file", "" );                 // File name for exodus mesh, if default no mesh is written
        tGENParameterList.insert( "geometry_field_file", "" );              // Base file name (without extension) for saving geometry fields
        tGENParameterList.insert( "time_offset", 0.0 );                     // Time offset for writing files in optimization process

        // ADVs/IQIs
        tGENParameterList.insert( "initial_advs", "" );          // Initial ADVs, string converted into vector
        tGENParameterList.insert( "advs_size", 0 );              // Specify size and fill value for ADVs in addition to
        tGENParameterList.insert( "initial_advs_fill", 0.0 );    // explicitly defined ADVs (above)
        tGENParameterList.insert( "lower_bounds", "" );          // Lower bounds on advs, string converted into vector
        tGENParameterList.insert( "lower_bounds_fill", 0.0 );    // Fill value for lower bounds up to size of ADV vector
        tGENParameterList.insert( "upper_bounds", "" );          // Upper bounds on advs, string converted into vector
        tGENParameterList.insert( "upper_bounds_fill", 0.0 );    // Fill value for upper bounds up to size of ADV vector
        tGENParameterList.insert( "IQI_types", "" );             // Requested IQI types for sensitivity analysis
        tGENParameterList.insert( "PDV_types", "" );             // Requested PDV types for sensitivity analysis

        // Phase table
        tGENParameterList.insert( "phase_table", "" );             // Construct phase table directly
        tGENParameterList.insert( "phase_function_name", "" );     // User-defined function for determining phase indices
        tGENParameterList.insert( "number_of_phases", 0 );         // Number of phases for a user-defined phase function
        tGENParameterList.insert( "print_phase_table", false );    // Whether to print the phase table

        // diagnostics
        tGENParameterList.insert( "diagnostics", false );
        tGENParameterList.insert( "diagnostics_id", "" );
        tGENParameterList.insert( "diagnostics_path", "" );

        return tGENParameterList;
    }

    /**
     * Inserts all general field parameters to the given parameter list.
     *
     * @param aDesignParameterList Parameter list for a design field
     */
    static void insert_field_parameters( ParameterList& aDesignParameterList )
    {
        aDesignParameterList.insert( "name", "" );                      // Name of field
        aDesignParameterList.insert( "field_type", "" );                // Type of field
        aDesignParameterList.insert( "dependencies", "" );              // Names of other fields that this field depends on
        aDesignParameterList.insert( "field_variable_indices", "" );    // Indices of field variables to fill
        aDesignParameterList.insert( "adv_indices", "" );               // ADVs used to fill in variables
        aDesignParameterList.insert( "constant_parameters", "" );       // Remaining geometry parameters that are constant
    }

    /**
     * Creates a field parameter list for creation before being added to a level set geometry or a property
     *
     * @return Parameter list with only field parameters
     */
    inline ParameterList create_field_parameter_list()
    {
        ParameterList tParameterList;
        insert_field_parameters( tParameterList );
        return tParameterList;
    }

    /**
     * Inserts all parameters related to design fields to the given parameter list.
     *
     * @param aDesignParameterList Parameter list for a design (level set geometry or property)
     * @param aIncludeField Whether or not to include field parameters. If not, an existing field name must be assigned.
     */
    static void insert_design_field_parameters( ParameterList& aDesignParameterList, bool aIncludeField )
    {
        if ( aIncludeField )
        {
            insert_field_parameters( aDesignParameterList );
        }
        else
        {
            aDesignParameterList.insert( "assigned_field", "" );
        }
        aDesignParameterList.insert( "discretization_mesh_index", -2 );           // Index of B-spline mesh to put this field on (-2 = none, -1 = store)
        aDesignParameterList.insert( "discretization_lower_bound", -1.0 );        // Lower bound of level set field (if bspline_mesh_index >= 0)
        aDesignParameterList.insert( "discretization_upper_bound", 1.0 );         // Upper bound of level set field (if bspline_mesh_index >= 0)
        aDesignParameterList.insert( "use_multilinear_interpolation", false );    // Whether to use multilinear interpolation for derived node field values
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
     * Creates a parameter list with parameters shared by all designs
     *
     * @return Design parameter list
     */
    static ParameterList create_design_parameter_list()
    {
        ParameterList tDesignParameterList;

        tDesignParameterList.insert( "design_type", "" );                  // Insert the design type parameter
        tDesignParameterList.insert( "number_of_refinements", "" );        // Number of refinement steps using HMR
        tDesignParameterList.insert( "refinement_mesh_index", "" );        // Refinement pattern
        tDesignParameterList.insert( "refinement_function_index", -1 );    // Index of user-defined refinement function (-1 = default)

        return tDesignParameterList;
    }

    /**
     * Creates a parameter list for the construction of a geometry.
     *
     * @return Geometry parameter list
     */
    static ParameterList create_geometry_parameter_list()
    {
        ParameterList tGeometryParameterList = create_design_parameter_list();
        tGeometryParameterList.set( "design_type", "geometry" );             // Set the design type to a geometry
        tGeometryParameterList.insert( "geometry_type", "" );                // Insert the geometry type parameter
        tGeometryParameterList.insert( "intersection_tolerance", 1e-12 );    // Interface tolerance based on intersection distance

        return tGeometryParameterList;
    }

    //------------------------------------------------------------------------------

    /**
     * Creates a parameter list that can be used to construct a geometry field. Any number of these can be added to
     * the second cell (index 1) of then parameter lists for GEN.
     *
     * @param aIncludeField Whether or not to include field parameters. If not, an existing field name must be assigned.
     * @return Geometry parameter list
     */
    inline ParameterList
    create_level_set_geometry_parameter_list( bool aIncludeField = true )
    {
        ParameterList tLevelSetParameterList = create_geometry_parameter_list();    // Inserts all geometry parameters
        insert_design_field_parameters( tLevelSetParameterList, aIncludeField );    // Inserts all design parameters
        tLevelSetParameterList.set( "geometry_type", "level_set" );                 // Sets the geometry type to level set
        tLevelSetParameterList.insert( "intersection_mode", "LEVEL_SET" );          // Deprecated
        tLevelSetParameterList.insert( "isocontour_threshold", 0.0 );               // Level set isocontour level
        tLevelSetParameterList.insert( "isocontour_tolerance", 1e-12 );             // Interface tolerance based on geometry value

        return tLevelSetParameterList;
    }

    //------------------------------------------------------------------------------

    /**
     * Creates a parameter list that can be used to construct a surface mesh geometry. Any number of these can be added to
     * the second cell (index 1) of then parameter lists for GEN.
     *
     * @return Geometry parameter list
     */
    inline ParameterList
    create_surface_mesh_geometry_parameter_list( bool aIncludeField = false )
    {
        ParameterList tSurfaceMeshParameterList = create_geometry_parameter_list();    // Inserts all geometry parameters
        insert_design_field_parameters( tSurfaceMeshParameterList, aIncludeField );    // Inserts all design parameters
        tSurfaceMeshParameterList.insert( "offset", "0,0,0" );                         // offset all points in the geometry by this much
        tSurfaceMeshParameterList.insert( "scale", "1.0,1.0,1.0" );                    // scaling factor for all points in the geometry
        tSurfaceMeshParameterList.insert( "file_path", "" );                           // path to .obj file
        tSurfaceMeshParameterList.set( "geometry_type", "surface_mesh" );              // set the geometry type to surface mesh
        tSurfaceMeshParameterList.set( "intersection_tolerance", 1e-8 );             // interface tolerance for raycast computations
        tSurfaceMeshParameterList.erase( "field_type" );

        return tSurfaceMeshParameterList;
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
        ParameterList tGeometryParameterList = create_level_set_geometry_parameter_list();    // Level set geometry parameters
        insert_user_defined_field_parameters( tGeometryParameterList );                       // Parameters for reading user-defined field functions

        return tGeometryParameterList;
    }

    //------------------------------------------------------------------------------

    /**
     * Same as a geometry parameter list, but forces the user to specify the path to the voxel file
     *
     * @return Voxel geometry parameterlist
     */
    inline ParameterList
    create_voxel_geometry_parameter_list()
    {
        ParameterList tVoxelParameterList = create_geometry_parameter_list();

        tVoxelParameterList.set( "geometry_type", "voxel" );       // Set the geometry type to a voxel geometry
        tVoxelParameterList.insert( "voxel_field_file", "" );      // voxel file
        tVoxelParameterList.insert( "domain_dimensions", "" );     // domain size
        tVoxelParameterList.insert( "domain_offset", "" );         // domain offset
        tVoxelParameterList.insert( "grain_id_value_map", "" );    // grain id to value map

        return tVoxelParameterList;
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
        ParameterList tSDFParameterList = create_level_set_geometry_parameter_list();

        tSDFParameterList.set( "field_type", "sdf_field" );     // SDF field type
        tSDFParameterList.insert( "sdf_object_path", "" );      // obj file
        tSDFParameterList.insert( "sdf_object_offset", "" );    // offset of object
        tSDFParameterList.insert( "sdf_shift", 0.0 );           // sdf shift

        return tSDFParameterList;
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
        ParameterList tImageSDFParameterList = create_level_set_geometry_parameter_list();

        tImageSDFParameterList.set( "field_type", "image_sdf" );            // sdf field generated from image
        tImageSDFParameterList.insert( "image_file", "" );                  // image file (hdf5 format)
        tImageSDFParameterList.insert( "image_dimensions", "" );            // domain size
        tImageSDFParameterList.insert( "image_offset", "" );                // domain offset
        tImageSDFParameterList.insert( "image_sdf_scaling", 0.0 );          // sdf scaling factor (0: automatic scaling)
        tImageSDFParameterList.insert( "image_sdf_shift", 0.0 );            // sdf shift value
        tImageSDFParameterList.insert( "image_sdf_default", -1.0 );         // sdf value outside image
        tImageSDFParameterList.insert( "image_sdf_interpolate", false );    // whether sdf value is interpolated

        return tImageSDFParameterList;
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
        ParameterList tNodalFieldParameterList = create_level_set_geometry_parameter_list();

        tNodalFieldParameterList.set( "field_type", "nodal_field_from_file" );    // field defined on mesh file
        tNodalFieldParameterList.insert( "file_name", "" );                       // file name
        tNodalFieldParameterList.insert( "field_name", "" );                      // field name
        tNodalFieldParameterList.insert( "file_format", "exodus" );               // file format
        tNodalFieldParameterList.insert( "offset", 0.0 );                         // offset of field value

        return tNodalFieldParameterList;
    }

    //------------------------------------------------------------------------------

    /**
     * Creates a parameter list that will be used to construct a field array.
     *
     * @return Swiss cheese slice parameter list
     */
    inline ParameterList
    create_field_array_parameter_list()
    {
        ParameterList tFieldArrayParameterList = create_level_set_geometry_parameter_list();    // TODO not tied to geometry

        // Lower and upper bounds on the array with respect to the reference coordinates of the field
        tFieldArrayParameterList.insert( "lower_bound_x", 0.0 );
        tFieldArrayParameterList.insert( "upper_bound_x", 0.0 );
        tFieldArrayParameterList.insert( "lower_bound_y", 0.0 );
        tFieldArrayParameterList.insert( "upper_bound_y", 0.0 );
        tFieldArrayParameterList.insert( "lower_bound_z", 0.0 );
        tFieldArrayParameterList.insert( "upper_bound_z", 0.0 );

        // Number of fields to create in each dimension.
        // If any are equal to zero, GEN will create the maximum number of fields the minimum spacing will allow.
        tFieldArrayParameterList.insert( "number_of_fields_x", 0 );
        tFieldArrayParameterList.insert( "number_of_fields_y", 0 );
        tFieldArrayParameterList.insert( "number_of_fields_z", 0 );

        // Minimum spacing between reference coordinates of the field
        tFieldArrayParameterList.insert( "minimum_spacing_x", 0.0 );
        tFieldArrayParameterList.insert( "minimum_spacing_y", 0.0 );
        tFieldArrayParameterList.insert( "minimum_spacing_z", 0.0 );

        // Amount that the reference points will be shifted over for each subsequent row in the specified direction.
        tFieldArrayParameterList.insert( "offset_per_row_x", 0.0 );
        tFieldArrayParameterList.insert( "offset_per_row_y", 0.0 );
        tFieldArrayParameterList.insert( "offset_per_row_z", 0.0 );

        // Whether to use the minimum value of all fields in the array or maximum
        tFieldArrayParameterList.insert( "minimum", true );

        return tFieldArrayParameterList;
    }

    //------------------------------------------------------------------------------

    /**
     * Creates a parameter list that can be used to construct a property field. Any number of these can be added to
     * the third cell (index 2) of the parameter lists for GEN.
     *
     * @param aIncludeField Whether or not to include field parameters. If not, an existing field name must be assigned.
     * @return GEN property parameter list
     */
    inline ParameterList
    create_gen_property_parameter_list( bool aIncludeField = true )
    {
        ParameterList tPropertyParameterList = create_design_parameter_list();      // Create a design parameter list
        tPropertyParameterList.set( "design_type", "property" );                    // Set the design type to a property
        insert_design_field_parameters( tPropertyParameterList, aIncludeField );    // Inserts all design field parameters
        tPropertyParameterList.insert( "pdv_type", "" );                            // The type of PDV that this property will be assigned to
        tPropertyParameterList.insert( "pdv_mesh_type", "interpolation" );          // Mesh type for assigning PDVs
        tPropertyParameterList.insert( "pdv_mesh_set_names", "" );                  // Mesh set names for assigning PDVs
        tPropertyParameterList.insert( "pdv_mesh_set_indices", "" );                // Mesh set indices for assigning PDVs

        return tPropertyParameterList;
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
        ParameterList tPropertyParameterList = create_gen_property_parameter_list();    // Create property parameter list
        insert_user_defined_field_parameters( tPropertyParameterList );                 // Parameters for reading user-defined field functions

        return tPropertyParameterList;
    }

    //------------------------------------------------------------------------------

}    // namespace moris::prm
