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
#include "cl_GEN_Pdv_Enums.hpp"

namespace moris
{
namespace ge
{

//typedef std::function< Matrix< DDRMat > ( moris::Cell< Matrix< DDRMat > >         & aCoeff,
//                                          moris::Cell< fem::Field_Interpolator* > & aDvFI ) > GENPropertyFunc;

typedef std::function< Matrix< DDRMat > ( moris::Cell< Matrix< DDRMat > >  & aCoeff ) > GENPropertyFunc;
//------------------------------------------------------------------------------
class GEN_Property
{

protected:
    // active dv types
    moris::Cell< moris::Cell< PDV > > mDvTypes;

    // active dv type map
    Matrix< DDSMat > mDvTypeMap;

    // dv field interpolators
    moris::Cell< moris::fem::Field_Interpolator* > mDvFI;

    // parameters
    moris::Cell< Matrix< DDRMat> > mParameters;

    // value function
    GENPropertyFunc mValFunction = nullptr;

    // dof and dv derivative functions
    moris::Cell< GENPropertyFunc > mDvDerFunctions;

    // geometry interpolator
    moris::fem::Geometry_Interpolator* mGeometryInterpolator = nullptr;

    // flag for evaluation
    bool mPropEval = true;
    moris::Cell< bool > mPropDvDerEval;

    // pdv type
    enum PDV mPdvType = PDV::UNDEFINED;

    // storage
    Matrix< DDRMat > mProp;
    moris::Cell< Matrix< DDRMat > > mPropDvDer;

public:

    //------------------------------------------------------------------------------
    GEN_Property(  ){};

    //------------------------------------------------------------------------------
    GEN_Property( enum PDV aPdvType ): mPdvType( aPdvType )
    {  };

    //------------------------------------------------------------------------------
    virtual ~GEN_Property( ){};

    //------------------------------------------------------------------------------
    void set_pdv_type( enum PDV aPdvType )
    {
        mPdvType = aPdvType;
    }
    //------------------------------------------------------------------------------
    enum PDV get_pdv_type( )
    {
        MORIS_ASSERT( mPdvType == PDV::UNDEFINED, "cl_GEN_Property::get_pdv_type() - pdv type not set" );
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
    void set_val_function( GENPropertyFunc aValFunction )
    {
        mValFunction = aValFunction;
    };
    //------------------------------------------------------------------------------
    /**
     * get val function
     * @param[ out ] mValFunction function for property evaluation
     */
    const GENPropertyFunc & get_val_function( ) const
    {
        return mValFunction;
    };
    //------------------------------------------------------------------------------
    /**
     * set dv derivative functions
     * @param[ in ] aDvDerFunctions list function for property derivatives wrt dv
     */
    void set_dv_derivative_functions( moris::Cell< GENPropertyFunc > aDvDerFunctions )
    {
        mDvDerFunctions = aDvDerFunctions;
    };
    //------------------------------------------------------------------------------
    /**
     * get dv derivative functions
     * @param[ out ] aDvDerFunctions list function for property derivatives wrt dv
     */
    const moris::Cell< GENPropertyFunc > & get_dv_derivative_functions( ) const
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
    /**
     * set a list of dv types
     * @param[ in ] aDvTypes list of dv type
     */
    void set_dv_type_list( const moris::Cell< moris::Cell< PDV > > & aDvTypes )
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
    const moris::Cell< moris::Cell< PDV > > & get_dv_type_list( ) const
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
    bool check_dv_dependency( const moris::Cell< PDV > aDvType );
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
    /**
     * get the property derivatives wrt a design variable
     * @param[ in ]  aDvType   cell of dv type
     * @param[ out ] adPropdDV matrix with derivative wrt to the dv type
     */
    const Matrix< DDRMat > & dPropdDV( const moris::Cell< PDV > aDvType );
    //------------------------------------------------------------------------------
    /**
     * evaluate property derivatives wrt a design variable
     * @param[ in ] aDvType cell of dv type
     */
    void eval_dPropdDV( const moris::Cell< PDV > aDvType );
    //------------------------------------------------------------------------------

};

}   // end ge namespace
}   // end moris namespace



#endif /* PROJECTS_GEN_SRC_PROPERTY_CL_GEN_PROPERTY_HPP_ */
