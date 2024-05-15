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

#include "cl_Parameter_List.hpp"
#include "GEN_Data_Types.hpp"

namespace moris::prm
{

    /**
     * Creates a parameter list used for the construction of the geometry engine. One of these parameter lists is
     * always required as the only parameter list in the first cell (index 0).
     *
     * @return GEN parameter list
     */
    inline Parameter_List
    create_gen_parameter_list()
    {
        Parameter_List tGENParameterList;

        // Level set parameters
        tGENParameterList.insert( "output_mesh_file", "" );                 // File name for exodus mesh, if default no mesh is written
        tGENParameterList.insert( "geometry_field_file", "" );              // Base file name (without extension) for saving geometry fields
        tGENParameterList.insert( "time_offset", 0.0 );                     // Time offset for writing files in optimization process

        // ADVs/IQIs
        tGENParameterList.insert( "initial_advs", Vector< real >() );    // Initial ADV vector
        tGENParameterList.insert( "advs_size", 0 );                      // Specify size and fill value for ADVs in addition to
        tGENParameterList.insert( "initial_advs_fill", 0.0 );            // explicitly defined ADVs (above)
        tGENParameterList.insert( "lower_bounds", Vector< real >() );    // Lower bounds on advs, string converted into vector
        tGENParameterList.insert( "lower_bounds_fill", 0.0 );            // Fill value for lower bounds up to size of ADV vector
        tGENParameterList.insert( "upper_bounds", Vector< real >() );    // Upper bounds on advs, string converted into vector
        tGENParameterList.insert( "upper_bounds_fill", 0.0 );            // Fill value for upper bounds up to size of ADV vector
        tGENParameterList.insert( "IQI_types", "" );                     // Requested IQI types for sensitivity analysis
        tGENParameterList.insert( "PDV_types", "" );                     // Requested PDV types for sensitivity analysis

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
     * @param aParameterList Parameter list for a field
     * @param aFieldType Type of field to create
     */
    static void insert_field_parameters( Parameter_List& aParameterList, gen::Field_Type aFieldType )
    {
        if ( aFieldType != gen::Field_Type::NONE )
        {
            // Set field type
            aParameterList.insert( "field_type", gen::Field_Type::NONE, gen::Field_Type::CONSTANT, gen::Field_Type::USER_DEFINED );
            aParameterList.set( "field_type", aFieldType );

            // General field parameters
            aParameterList.insert( "name", "" );                                    // Name of field
            aParameterList.insert( "field_variable_indices", Vector< uint >() );    // Indices of field variables to fill
            aParameterList.insert( "adv_indices", Vector< uint >() );               // ADVs used to fill in variables
            aParameterList.insert( "constant_parameters", Vector< real >() );       // Remaining geometry parameters that are constant

            // Insert specific field parameters
            switch ( aFieldType )
            {
                case gen::Field_Type::NONE:
                    break;
                case gen::Field_Type::CONSTANT:
                {
                    aParameterList.insert< Design_Variable >( "constant", 0.0 );
                }
                case gen::Field_Type::LINE:
                {
                    aParameterList.insert< Design_Variable >( "center_x", 0.0 );
                    aParameterList.insert< Design_Variable >( "center_y", 0.0 );
                    aParameterList.insert< Design_Variable >( "normal_x", 1.0 );
                    aParameterList.insert< Design_Variable >( "normal_y", 0.0 );
                    break;
                }
                case gen::Field_Type::CIRCLE:
                {
                    aParameterList.insert< Design_Variable >( "center_x", 0.0 );
                    aParameterList.insert< Design_Variable >( "center_y", 0.0 );
                    aParameterList.insert< Design_Variable >( "radius", 1.0, 0.0, MORIS_REAL_MAX );
                    break;
                }
                case gen::Field_Type::SUPERELLIPSE:
                {
                    aParameterList.insert< Design_Variable >( "center_x", 0.0 );
                    aParameterList.insert< Design_Variable >( "center_y", 0.0 );
                    aParameterList.insert< Design_Variable >( "semidiameter_x", 1.0, 0.0, MORIS_REAL_MAX );
                    aParameterList.insert< Design_Variable >( "semidiameter_y", 1.0, 0.0, MORIS_REAL_MAX );
                    aParameterList.insert( "exponent", 2.0 );
                    break;
                }
                case gen::Field_Type::PLANE:
                {
                    aParameterList.insert< Design_Variable >( "center_x", 0.0 );
                    aParameterList.insert< Design_Variable >( "center_y", 0.0 );
                    aParameterList.insert< Design_Variable >( "center_z", 0.0 );
                    aParameterList.insert< Design_Variable >( "normal_x", 1.0 );
                    aParameterList.insert< Design_Variable >( "normal_y", 0.0 );
                    aParameterList.insert< Design_Variable >( "normal_z", 0.0 );
                    break;
                }
                case gen::Field_Type::SPHERE:
                {
                    aParameterList.insert< Design_Variable >( "center_x", 0.0 );
                    aParameterList.insert< Design_Variable >( "center_y", 0.0 );
                    aParameterList.insert< Design_Variable >( "center_z", 0.0 );
                    aParameterList.insert< Design_Variable >( "radius", 1.0, 0.0, MORIS_REAL_MAX );
                    break;
                }
                case gen::Field_Type::SUPERELLIPSOID:
                {
                    aParameterList.insert< Design_Variable >( "center_x", 0.0 );
                    aParameterList.insert< Design_Variable >( "center_y", 0.0 );
                    aParameterList.insert< Design_Variable >( "center_z", 0.0 );
                    aParameterList.insert< Design_Variable >( "semidiameter_x", 1.0, 0.0, MORIS_REAL_MAX );
                    aParameterList.insert< Design_Variable >( "semidiameter_y", 1.0, 0.0, MORIS_REAL_MAX );
                    aParameterList.insert< Design_Variable >( "semidiameter_z", 1.0, 0.0, MORIS_REAL_MAX );
                    aParameterList.insert( "exponent", 2.0 );
                    break;
                }
                case gen::Field_Type::SCALED_FIELD:
                {
                    aParameterList.insert( "dependencies", Vector< std::string >() );
                    aParameterList.insert< Design_Variable >( "scaling_factor", 1.0 );
                    break;
                }
                case gen::Field_Type::COMBINED_FIELDS:
                {
                    aParameterList.insert( "dependencies", Vector< std::string >() );
                    aParameterList.insert( "use_minimum", true );
                    break;
                }
                case gen::Field_Type::NODAL:
                    break;
                case gen::Field_Type::NODAL_FROM_FILE:
                {
                    aParameterList.insert( "file_name", "" );            // file name
                    aParameterList.insert( "field_name", "" );           // field name
                    aParameterList.insert( "file_format", "exodus" );    // file format
                    aParameterList.insert( "offset", 0.0 );              // offset of field value
                    break;
                }
                case gen::Field_Type::SIGNED_DISTANCE_OBJECT:
                {
                    aParameterList.insert( "sdf_object_path", "" );      // obj file
                    aParameterList.insert( "sdf_object_offset", "" );    // offset of object
                    aParameterList.insert( "sdf_shift", 0.0 );           // sdf shift
                    break;
                }
                case gen::Field_Type::SIGNED_DISTANCE_IMAGE:
                {
                    aParameterList.insert( "image_file", "" );                  // image file (hdf5 format)
                    aParameterList.insert( "image_dimensions", "" );            // domain size
                    aParameterList.insert( "image_offset", "" );                // domain offset
                    aParameterList.insert( "image_sdf_scaling", 0.0 );          // sdf scaling factor (0: automatic scaling)
                    aParameterList.insert( "image_sdf_shift", 0.0 );            // sdf shift value
                    aParameterList.insert( "image_sdf_default", -1.0 );         // sdf value outside image
                    aParameterList.insert( "image_sdf_interpolate", false );    // whether sdf value is interpolated
                    break;
                }
                case gen::Field_Type::USER_DEFINED:
                {
                    aParameterList.insert( "field_function_name", "" );          // Function name for evaluating the geometry field
                    aParameterList.insert( "sensitivity_function_name", "" );    // Function name for evaluating the sensitivity of the field
                    break;
                }
            }
        }
    }

    /**
     * Creates a field parameter list for creation before being added to a level set geometry or a property
     *
     * @param aFieldType Type of field to create
     * @return Parameter list with only field parameters
     */
    inline Parameter_List create_field_parameter_list( gen::Field_Type aFieldType )
    {
        Parameter_List tParameterList;
        insert_field_parameters( tParameterList, aFieldType );
        return tParameterList;
    }

    /**
     * Inserts all parameters related to design fields to the given parameter list.
     *
     * @param aDesignParameterList Parameter list for a design (level set geometry or property)
     * @param aFieldType Type of field to create
     */
    static void insert_design_field_parameters( Parameter_List& aDesignParameterList, gen::Field_Type aFieldType )
    {
        insert_field_parameters( aDesignParameterList, aFieldType );
        aDesignParameterList.insert( "discretization_mesh_index", -2 );           // Index of B-spline mesh to put this field on (-2 = none, -1 = store)
        aDesignParameterList.insert( "discretization_lower_bound", -1.0 );        // Lower bound of level set field (if bspline_mesh_index >= 0)
        aDesignParameterList.insert( "discretization_upper_bound", 1.0 );         // Upper bound of level set field (if bspline_mesh_index >= 0)
        aDesignParameterList.insert( "use_multilinear_interpolation", false );    // Whether to use multilinear interpolation for derived node field values
    }

    /**
     * Creates a parameter list with parameters shared by all designs
     *
     * @return Design parameter list
     */
    static Parameter_List create_design_parameter_list()
    {
        Parameter_List tDesignParameterList;

        tDesignParameterList.insert( "design_type", "" );               // Insert the design type parameter
        tDesignParameterList.insert( "number_of_refinements", "" );     // Number of refinement steps using HMR
        tDesignParameterList.insert( "refinement_mesh_index", "" );     // Refinement pattern
        tDesignParameterList.insert( "refinement_function_index", -1 ); // Index of user-defined refinement function (-1 = default)

        return tDesignParameterList;
    }

    /**
     * Creates a parameter list for the construction of a geometry.
     *
     * @return Geometry parameter list
     */
    static Parameter_List create_geometry_parameter_list()
    {
        Parameter_List tGeometryParameterList = create_design_parameter_list();
        tGeometryParameterList.set( "design_type", "geometry" ); // Set the design type to a geometry
        tGeometryParameterList.insert( "geometry_type", "" );    // Insert the geometry type parameter

        return tGeometryParameterList;
    }

    //------------------------------------------------------------------------------

    /**
     * Creates a parameter list that can be used to construct a geometry field. Any number of these can be added to
     * the second cell (index 1) of then parameter lists for GEN.
     *
     * @param aFieldType Type of field to create
     * @return Geometry parameter list
     */
    inline Parameter_List
    create_level_set_geometry_parameter_list( gen::Field_Type aFieldType )
    {
        Parameter_List tLevelSetParameterList = create_geometry_parameter_list();  // Inserts all geometry parameters
        insert_design_field_parameters( tLevelSetParameterList, aFieldType );      // Inserts all design parameters
        tLevelSetParameterList.set( "geometry_type", "level_set" );                // Sets the geometry type to level set
        tLevelSetParameterList.insert( "isocontour_threshold", 0.0 );              // Level set isocontour level
        tLevelSetParameterList.insert( "isocontour_tolerance", 1e-12 );            // Interface tolerance based on geometry value
        tLevelSetParameterList.insert( "intersection_tolerance", 1e-12 );          // Interface tolerance based on intersection distance

        return tLevelSetParameterList;
    }

    //------------------------------------------------------------------------------

    /**
     * Creates a parameter list that can be used to construct a surface mesh geometry. Any number of these can be added to
     * the second cell (index 1) of then parameter lists for GEN.
     *
     * @return Geometry parameter list
     */
    inline Parameter_List
    create_surface_mesh_geometry_parameter_list()
    {
        Parameter_List tSurfaceMeshParameterList = create_geometry_parameter_list();           // Inserts all geometry parameters
        insert_design_field_parameters( tSurfaceMeshParameterList, gen::Field_Type::NONE );    // Inserts all design parameters
        tSurfaceMeshParameterList.insert( "offset", "0,0,0" );                                 // offset all points in the geometry by this much
        tSurfaceMeshParameterList.insert( "scale", "1.0,1.0,1.0" );                            // scaling factor for all points in the geometry
        tSurfaceMeshParameterList.insert( "file_path", "" );                                   // path to .obj file
        tSurfaceMeshParameterList.set( "geometry_type", "surface_mesh" );                      // set the geometry type to surface mesh
        tSurfaceMeshParameterList.erase( "field_type" );

        return tSurfaceMeshParameterList;
    }

    //------------------------------------------------------------------------------

    /**
     * Same as a geometry parameter list, but forces the user to specify the path to the voxel file
     *
     * @return Voxel geometry parameterlist
     */
    inline Parameter_List
    create_voxel_geometry_parameter_list()
    {
        Parameter_List tVoxelParameterList = create_geometry_parameter_list();

        tVoxelParameterList.set( "geometry_type", "voxel" );                    // Set the geometry type to a voxel geometry
        tVoxelParameterList.insert( "voxel_field_file", "" );                   // voxel file
        tVoxelParameterList.insert( "domain_dimensions", Vector< real >() );    // domain size
        tVoxelParameterList.insert( "domain_offset", Vector< real >() );        // domain offset

        return tVoxelParameterList;
    }

    //------------------------------------------------------------------------------

    /**
     * Creates a parameter list that will be used to construct a field array.
     *
     * @return Swiss cheese slice parameter list
     */
    inline Parameter_List
    create_field_array_parameter_list( gen::Field_Type aFieldType )
    {
        Parameter_List tFieldArrayParameterList = create_level_set_geometry_parameter_list( aFieldType );    // TODO not tied to geometry

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
    inline Parameter_List
    create_gen_property_parameter_list( gen::Field_Type aFieldType )
    {
        Parameter_List tPropertyParameterList = create_design_parameter_list();  // Create a design parameter list
        tPropertyParameterList.set( "design_type", "property" );                 // Set the design type to a property
        insert_design_field_parameters( tPropertyParameterList, aFieldType );    // Inserts all design field parameters
        tPropertyParameterList.insert( "pdv_type", "" );                         // The type of PDV that this property will be assigned to
        tPropertyParameterList.insert( "pdv_mesh_type", "interpolation" );       // Mesh type for assigning PDVs
        tPropertyParameterList.insert( "pdv_mesh_set_names", "" );               // Mesh set names for assigning PDVs
        tPropertyParameterList.insert( "pdv_mesh_set_indices", "" );             // Mesh set indices for assigning PDVs

        return tPropertyParameterList;
    }

    //------------------------------------------------------------------------------

}    // namespace moris::prm
