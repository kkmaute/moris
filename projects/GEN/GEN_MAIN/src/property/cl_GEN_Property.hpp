/*
 * cl_GEN_Property.hpp
 *
 *  Created on: Dec 4, 2019
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_SRC_PROPERTY_CL_GEN_PROPERTY_HPP_
#define PROJECTS_GEN_SRC_PROPERTY_CL_GEN_PROPERTY_HPP_

// FEM includes
#include "cl_FEM_Field_Interpolator.hpp"

// GE includes
//#include "cl_GEN_Enums.hpp"
#include "../projects/GEN/GEN_CORE/src/cl_GEN_Dv_Enums.hpp"

namespace moris
{
namespace ge
{
//typedef std::function< Matrix< DDRMat > ( moris::Cell< Matrix< DDRMat > >         & aCoeff,
//                                          moris::Cell< fem::Field_Interpolator* > & aDofFI,
//                                          moris::Cell< fem::Field_Interpolator* > & aDvFI,
//                                          fem::Geometry_Interpolator              * aGeometryInterpolator ) > PropertyFunc;

//typedef std::function< Matrix< DDRMat > ( moris::Cell< Matrix< DDRMat > >         & aCoeff,
//                                          moris::Cell< fem::Field_Interpolator* > & aDvFI,
//                                          fem::Geometry_Interpolator              * aGeometryInterpolator ) > PropertyFunc;

typedef std::function< Matrix< DDRMat > ( moris::Cell< Matrix< DDRMat > >         & aCoeff,
                                          moris::Cell< fem::Field_Interpolator* > & aDvFI ) > PropertyFunc;
//------------------------------------------------------------------------------
class GEN_Property
{
protected:
//    // active dof types
//    moris::Cell< moris::Cell< MSI::Dof_Type > > mDofTypes;

//    // single dof type
//    MSI::Dof_Type mDof = MSI::Dof_Type::END_ENUM;

//    // active dof type map
//    Matrix< DDSMat > mDofTypeMap;

//    // dof field interpolators
//    moris::Cell< moris::fem::Field_Interpolator* > mDofFI;

    // active dv types
    moris::Cell< moris::Cell< GEN_DV > > mDvTypes;

    // active dv type map
    Matrix< DDSMat > mDvTypeMap;

    // dv field interpolators
    moris::Cell< moris::fem::Field_Interpolator* > mDvFI;

    // parameters
    moris::Cell< Matrix< DDRMat> > mParameters;

    // value function
    PropertyFunc mValFunction = nullptr;

    // dof and dv derivative functions
//    moris::Cell< PropertyFunc > mDofDerFunctions;
    moris::Cell< PropertyFunc > mDvDerFunctions;

    // geometry interpolator
    moris::fem::Geometry_Interpolator* mGeometryInterpolator = nullptr;

    // flag for evaluation
    bool mPropEval = true;
//    moris::Cell< bool > mPropDofDerEval;
    moris::Cell< bool > mPropDvDerEval;

    // pdv type
    enum GEN_DV mPdvType = GEN_DV::END_ENUM;

    // storage
    Matrix< DDRMat > mProp;
//    moris::Cell< Matrix< DDRMat > > mPropDofDer;
    moris::Cell< Matrix< DDRMat > > mPropDvDer;

public:

    //------------------------------------------------------------------------------
    GEN_Property(  ){};

    //------------------------------------------------------------------------------
    GEN_Property( enum GEN_DV aPdvType ): mPdvType( aPdvType )
    {  };

    //------------------------------------------------------------------------------
    virtual ~GEN_Property( ){};

    //------------------------------------------------------------------------------
    void set_pdv_type( enum GEN_DV aPdvType )      //FIXME: use the 'set_dv_type()' function?
    {
        mPdvType = aPdvType;
    }

    enum GEN_DV get_pdv_type( )
    {
        MORIS_ASSERT( mPdvType == GEN_DV::END_ENUM, "cl_GEN_Property::get_pdv_type() - pdv type not set" );
        return mPdvType;
    }
    //------------------------------------------------------------------------------
    /**
     * set parameters
     * @param[ in ] aParameters list of parameters for property evaluation
     */
    void set_parameters( moris::Cell< moris::Matrix< DDRMat > > aParameters )
    {
        mParameters = aParameters;
    };

    //------------------------------------------------------------------------------
    /**
     * get parameters
     * @param[ out ] mParameters list of parameters for property evaluation
     */
    const moris::Cell< moris::Matrix< DDRMat > > & get_parameters( ) const
    {
        return mParameters;
    };

    //------------------------------------------------------------------------------
    /**
     * set val function
     * @param[ in ] aValFunction function for property evaluation
     */
    void set_val_function( PropertyFunc aValFunction )
    {
        mValFunction = aValFunction;
    };

    //------------------------------------------------------------------------------
    /**
     * get val function
     * @param[ out ] mValFunction function for property evaluation
     */
    const PropertyFunc & get_val_function( ) const
    {
        return mValFunction;
    };
    //------------------------------------------------------------------------------
//    /**
//     * set dof derivative functions
//     * @param[ in ] aDofDerFunctions list function for property derivatives wrt dof
//     */
//    void set_dof_derivative_functions( moris::Cell< PropertyFunc > aDofDerFunctions )
//    {
//        mDofDerFunctions = aDofDerFunctions;
//    };

    //------------------------------------------------------------------------------
//    /**
//     * get dof derivative functions
//     * @param[ out ] mDofDerFunctions list function for property derivatives wrt dof
//     */
//    const moris::Cell< PropertyFunc > & get_dof_derivative_functions() const
//    {
//        return mDofDerFunctions;
//    };

    //------------------------------------------------------------------------------
    /**
     * set dv derivative functions
     * @param[ in ] aDvDerFunctions list function for property derivatives wrt dv
     */
    void set_dv_derivative_functions( moris::Cell< PropertyFunc > aDvDerFunctions )
    {
        mDvDerFunctions = aDvDerFunctions;
    };

    //------------------------------------------------------------------------------
    /**
     * get dv derivative functions
     * @param[ out ] aDvDerFunctions list function for property derivatives wrt dv
     */
    const moris::Cell< PropertyFunc > & get_dv_derivative_functions( ) const
    {
        return mDvDerFunctions;
    };

    //------------------------------------------------------------------------------
    /**
     * set geometry interpolator
     * @param[ in ] aGeometryInterpolator geometry interpolator pointer
     */
    void set_geometry_interpolator( fem::Geometry_Interpolator * aGeometryInterpolator )
    {
        mGeometryInterpolator = aGeometryInterpolator;
    };

    //------------------------------------------------------------------------------
    /**
     * get geometry interpolator
     * @param[ out ] mGeometryInterpolator geometry interpolator pointer
     */
    fem::Geometry_Interpolator * get_geometry_interpolator( ) const
    {
        return mGeometryInterpolator;
    };
    //------------------------------------------------------------------------------
//    /**
//     * set a list of active dof types
//     * @param[ in ] aDofTypes list of dof types
//     */
//    void set_dof_type_list( moris::Cell< moris::Cell< MSI::Dof_Type > > aDofTypes )
//    {
//        // set dof type list
//        mDofTypes = aDofTypes;
//
//        // build dof type map
//        this->build_dof_type_map();
//
//        // number of dof types
//        uint tNumDofTypes = mDofTypes.size();
//
//        // set mDofDerFunctions size
//        mDofDerFunctions.assign( tNumDofTypes, nullptr );
//
//        // set mPropDofDerEval size
//        mPropDofDerEval.assign( tNumDofTypes, true );
//
//        // set mPropDofDer size
//        mPropDofDer.resize( tNumDofTypes );
//    };

    //------------------------------------------------------------------------------
//    /*
//     * set dof type for the single dof type case
//     * @param[ in ] aDof a dof type
//     */
//    void set_dof_type( MSI::Dof_Type aDof )
//    {
//        mDof = aDof;
//    };

    //------------------------------------------------------------------------------
//    /*
//     * get dof type for the single dof type case
//     * @param[ out ] mDof a dof type
//     */
//    MSI::Dof_Type get_dof_type(  )
//    {
//        return mDof;
//    };

    //------------------------------------------------------------------------------
//     /**
//      * return a list of active dof types
//      * @param[ out ] mDofTypes list of dof types
//      */
//     const moris::Cell< moris::Cell< MSI::Dof_Type > > & get_dof_type_list() const
//     {
//         return mDofTypes;
//     };

    //------------------------------------------------------------------------------
//    /**
//     * build a dof type map
//     */
//    void build_dof_type_map();

    //------------------------------------------------------------------------------
//    /**
//     * get a dof type map
//     * @param[ out ] mDofTypeMap map of the dof types
//     */
//    const moris::Matrix< DDSMat > & get_dof_type_map() const
//    {
//        return mDofTypeMap;
//    }

    //------------------------------------------------------------------------------
//     /**
//      * check if the property depends on a particular group of dof type
//      * @param[ in ]  aDofType cell of dof type
//      * @param[ out ] aBool    boolean, true if dependency on the dof type
//      */
//     bool check_dof_dependency( const moris::Cell< MSI::Dof_Type > aDofType );

    //------------------------------------------------------------------------------
//    /**
//     * set dof field interpolators
//     * @param[ in ] aFieldInterpolators cell of dof field interpolator pointers
//     */
//    void set_dof_field_interpolators( moris::Cell< moris::fem::Field_Interpolator* > & aFieldInterpolators );

    //------------------------------------------------------------------------------
//    /**
//     * get dof field interpolators
//     * @param[ out ] mDofFI cell of dof field interpolator pointers
//     */
//    const moris::Cell< fem::Field_Interpolator* > & get_dof_field_interpolators() const
//    {
//        return mDofFI;
//    }

    //------------------------------------------------------------------------------
//    /**
//     * check that dof field interpolators are assigned
//     */
//    void check_dof_field_interpolators();

    //------------------------------------------------------------------------------
    /**
     * set a list of dv types
     * @param[ in ] aDvTypes list of dv type
     */
    void set_dv_type_list( const moris::Cell< moris::Cell< GEN_DV > > & aDvTypes )
    {
        // set dv type list
        mDvTypes = aDvTypes;

        // build a dv type map
        this->build_dv_type_map();

        // number of dv types
        uint tNumDvTypes = mDvTypes.size();

        // set mDvDerFunctions size
        mDvDerFunctions.resize( tNumDvTypes, nullptr );

        // set mPropDvDerEval size
        mPropDvDerEval.assign( tNumDvTypes, true );

        // set mPropdvDer size
        mPropDvDer.resize( tNumDvTypes );
    };

    //------------------------------------------------------------------------------
    /**
     * return a list of dv types
     * @param[ out ] mDvTypes list of dv type
     */
    const moris::Cell< moris::Cell< GEN_DV > > & get_dv_type_list( ) const
    {
        return mDvTypes;
    };

    //------------------------------------------------------------------------------
    /**
     * build a dv type map
     */
    void build_dv_type_map( );

    //------------------------------------------------------------------------------
    /**
     * get a dv type map
     * @param[ out ] mDvTypeMap map of the dv types
     */
    const moris::Matrix< DDSMat > & get_dv_type_map( ) const
    {
        return mDvTypeMap;
    }

    //------------------------------------------------------------------------------
    /**
     * check if the property depends on a particular group of dv types
     * @param[ in ]  aDvType cell of dv type
     * @param[ out ] aBool   boolean, true if dependency on the dv type
     */
    bool check_dv_dependency( const moris::Cell< GEN_DV > aDvType );

    //------------------------------------------------------------------------------
    /**
     * set dv field interpolators
     * @param[ in ] aFieldInterpolators cell of dv field interpolator pointers
     */
    void set_dv_field_interpolators( moris::Cell< moris::fem::Field_Interpolator* > & aFieldInterpolators );

    //------------------------------------------------------------------------------
    /**
     * get dv field interpolators
     * @param[ out ] mDvFI cell of dv field interpolator pointers
     */
    const moris::Cell< fem::Field_Interpolator* > & get_dv_field_interpolators( ) const
    {
        return mDvFI;
    }

    //------------------------------------------------------------------------------
    /**
     * check that dv field interpolators are assigned
     */
    void check_dv_field_interpolators( );

    //------------------------------------------------------------------------------
    /**
     * reset evaluation flags
     */
    void reset_eval_flags( )
    {
        // reset the property value
        mPropEval = true;

//        // reset the property derivatives wrt dof type
//        mPropDofDerEval.assign( mDofTypes.size(), true );

        // reset the property derivatives wrt dv type
        mPropDvDerEval.assign( mDvTypes.size(), true );
    }

    //------------------------------------------------------------------------------
    /**
     * get the property value
     * @param[ out ] aVal matrix with property value
     */
    const Matrix< DDRMat > & val( );

    //------------------------------------------------------------------------------
    /**
     * evaluate property in terms of coefficients and variables
     */
    void eval_Prop( );

    //------------------------------------------------------------------------------
//    /**
//     * get property derivatives wrt a dof type
//     * @param[ in ]  aDofType   cell of dof type
//     * @param[ out ] adPropdDOF matrix with derivative wrt to the dof type
//     */
//    const Matrix< DDRMat > & dPropdDOF( const moris::Cell< MSI::Dof_Type > aDofType );

    //------------------------------------------------------------------------------
//    /**
//     * evaluate property derivatives wrt a dof type
//     * @param[ in ] aDofType cell of dof type
//     */
//    void eval_dPropdDOF( const moris::Cell< MSI::Dof_Type > aDofType );

    //------------------------------------------------------------------------------
    /**
     * get the property derivatives wrt a design variable
     * @param[ in ]  aDvType   cell of dv type
     * @param[ out ] adPropdDV matrix with derivative wrt to the dv type
     */
    const Matrix< DDRMat > & dPropdDV( const moris::Cell< GEN_DV > aDvType );

    //------------------------------------------------------------------------------
    /**
     * evaluate property derivatives wrt a design variable
     * @param[ in ] aDvType cell of dv type
     */
    void eval_dPropdDV( const moris::Cell< GEN_DV > aDvType );

    //------------------------------------------------------------------------------
//    /**
//     * get non unique dof type list
//     * @param[ in ] aDofType cell of dof type
//     */
//    void get_non_unique_dof_types( moris::Cell< MSI::Dof_Type > & aDofTypes )
//    {
//        // init counter
//        uint tCounter = 0;
//
//        // loop over dof types
//        for ( uint iDOF = 0; iDOF < mDofTypes.size(); iDOF++ )
//        {
//            // update counter
//            tCounter += mDofTypes( iDOF ).size();
//        }
//
//        // reserve memory for dof type list
//        aDofTypes.reserve( tCounter );
//
//        // loop over dof types
//        for ( uint iDOF = 0; iDOF < mDofTypes.size(); iDOF++ )
//        {
//            // populate the dof type list
//            aDofTypes.append( mDofTypes( iDOF ) );
//        }
//    }
    //------------------------------------------------------------------------------

};

}   // end ge namespace
}   // end moris namespace



#endif /* PROJECTS_GEN_SRC_PROPERTY_CL_GEN_PROPERTY_HPP_ */
