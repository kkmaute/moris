/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Library_IO_Standard.cpp
 *
 */

#include "cl_Library_IO_Standard.hpp"

namespace moris
{
    //------------------------------------------------------------------------------------------------------------------

    Library_IO_Standard::Library_IO_Standard()
            : Library_IO() // initialize base class data as usual
    {
        // set the type of this library
        mLibraryType = Library_Type::STANDARD;

        // list of supported parameter list types
        mSupportedParamListTypes = {
            Parameter_List_Type::OPT,
            Parameter_List_Type::HMR,
            Parameter_List_Type::STK,
            Parameter_List_Type::XTK,
            Parameter_List_Type::GEN,
            Parameter_List_Type::FEM,
            Parameter_List_Type::SOL,
            Parameter_List_Type::MSI,
            Parameter_List_Type::VIS,
            Parameter_List_Type::MIG,
            Parameter_List_Type::WRK,
            Parameter_List_Type::MORISGENERAL };
    }

    //------------------------------------------------------------------------------------------------------------------

    Library_IO_Standard::~Library_IO_Standard()
    {
        // do nothing extra
    }

    //------------------------------------------------------------------------------------------------------------------

    void
    Library_IO_Standard::finalize()
    {
        // check that an .xml input file has been specified
        MORIS_ERROR( mSoLibIsInitialized || mXmlParserIsInitialized,
                "Library_IO_Standard::finalize() - Neither an .xml nor a .so input file has been specified. "
                "At least one input file is required." );

        // load the standard parameters into the member variables
        this->load_all_standard_parameters();

        // if an .so file has been parsed, first use its parameters (if any were defined in it) to overwrite or add to the standard parameters
        if( mSoLibIsInitialized )
        {
            this->load_parameters_from_shared_object_library();
        }

        // load parameters from xml, overwrites parameters specified in either the standard parameters or an .so file if parsed
        if( mXmlParserIsInitialized )
        {
            this->load_parameters_from_xml();
        }

        // check the parameters for validity
        // TODO: this->check_parameters();

        // mark this library as finalized and lock it from modification
        mLibraryIsFinalized = true;

        // print receipt of the finalized library
        this->print_parameter_receipt( "./Parameter_Receipt.xml" ); // TODO: the file name and location should be user defineable
    }

    //------------------------------------------------------------------------------------------------------------------
    // STANDARD PARAMETER LIST FUNCTIONS
    //------------------------------------------------------------------------------------------------------------------

    void
    Library_IO_Standard::load_all_standard_parameters()
    {
        // FIXME: use PRM functions once everything is moved
        // do nothing for now
    }

    //------------------------------------------------------------------------------------------------------------------

} // namespace moris
