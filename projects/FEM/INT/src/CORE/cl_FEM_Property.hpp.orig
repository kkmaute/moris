/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Property.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_PROPERTY_HPP_
#define SRC_FEM_CL_FEM_PROPERTY_HPP_

#include "moris_typedefs.hpp"                     //MRS/COR/src
#include "cl_Matrix.hpp"                    //LNA/src
#include "cl_Vector.hpp"                    //MRS/CNT/src
#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_MSI_Dof_Type_Enums.hpp"        //FEM/MSI/src
#include "cl_FEM_Enums.hpp"                 //FEM/MSI/src
#include "fn_equal_to.hpp"

#include "GEN_Data_Types.hpp"

#include <functional>
#include <utility>

namespace moris::fem
{
    class Set;
    class Field_Interpolator_Manager;

    // definition of a property value/derivative function
    typedef std::function< void(
            moris::Matrix< moris::DDRMat >&         aPropMatrix,
            Vector< Matrix< DDRMat > >&             aCoeff,
            moris::fem::Field_Interpolator_Manager* aFIManager ) >
            PropertyFunc;

    // default constant property value function
    inline void
    Property_Function_Constant(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    //------------------------------------------------------------------------------
    /**
     * Property
     */
    class Property
    {
      protected:
        // field interpolator manager pointer
        Field_Interpolator_Manager* mFIManager = nullptr;

        // FEM set pointer
        Set* mSet = nullptr;

        // active dof types
        Vector< Vector< MSI::Dof_Type > > mDofTypes;

        // active dof type map
        Matrix< DDSMat > mDofTypeMap;

        // active dv types
        Vector< Vector< gen::PDV_Type > > mDvTypes;

        // active dv type map
        Matrix< DDSMat > mDvTypeMap;

        // active dv types
        Vector< Vector< mtk::Field_Type > > mFieldTypes;

        // active dv type map
        Matrix< DDSMat > mFieldTypeMap;

        // parameters
        Vector< Matrix< DDRMat > > mParameters;

        // value function
        PropertyFunc mValFunction = Property_Function_Constant;

        // space derivative functions
        Vector< PropertyFunc > mSpaceDerFunctions;

        // dof and dv derivative functions
        Vector< PropertyFunc > mDofDerFunctions;
        Vector< PropertyFunc > mDvDerFunctions;

        // storage
        Matrix< DDRMat >           mProp;
        Vector< Matrix< DDRMat > > mPropDofDer;
        Vector< Matrix< DDRMat > > mPropDvDer;
        Vector< Matrix< DDRMat > > mdPropdx;

        // property name
        std::string mName;

      private:
        // flag for evaluation
        bool                    mPropEval = true;
        moris::Matrix< DDBMat > mdPropdxEval;
        moris::Matrix< DDBMat > mPropDofDerEval;
        moris::Matrix< DDBMat > mPropDvDerEval;

        // flag for setting mValFunction, mSpaceDerFunction, mDofDerFunctions, mDvDerFunctions
        bool                    mSetValFunction = false;
        moris::Matrix< DDBMat > mSetSpaceDerFunctions;
        bool                    mSetDofDerFunctions = false;
        bool                    mSetDvDerFunctions  = false;

        //------------------------------------------------------------------------------

      public:
        //------------------------------------------------------------------------------
        /**
         * constructor
         */
        Property();

        //------------------------------------------------------------------------------
        /**
         * virtual destructor
         */
        virtual ~Property(){};

        //------------------------------------------------------------------------------
        /**
         * set name
         * param[ in ] aName a string for property name
         */
        void
        set_name( std::string aName )
        {
            mName = std::move( aName );
        }

        //------------------------------------------------------------------------------
        /**
         * get name
         * param[ out ] mName a string for property name
         */
        std::string
        get_name()
        {
            return mName;
        }

        //------------------------------------------------------------------------------
        /**
         * set field interpolator manager
         * @param[ in ] aFieldInterpolatorManager a field interpolator manager pointer
         */
        void set_field_interpolator_manager(
                Field_Interpolator_Manager* aFieldInterpolatorManager );

        //------------------------------------------------------------------------------
        /*
         * set member set pointer
         * @param[ in ] aSetPointer a FEM set pointer
         */
        void
        set_set_pointer( Set* aSetPointer )
        {
            mSet = aSetPointer;
        }

        //------------------------------------------------------------------------------
        /**
         * set parameters
         * @param[ in ] aParameters list of parameters for property evaluation
         */
        void
        set_parameters(
                const Vector< moris::Matrix< DDRMat > >& aParameters )
        {
            mParameters = aParameters;
        }

        //------------------------------------------------------------------------------
        /**
         * get parameters
         * @param[ out ] mParameters list of parameters for property evaluation
         */
        const Vector< moris::Matrix< DDRMat > >&
        get_parameters() const
        {
            return mParameters;
        }

        //------------------------------------------------------------------------------
        /**
         * set val function
         * @param[ in ] aValFunction function for property evaluation
         */
        void
        set_val_function( PropertyFunc aValFunction )
        {
            // set the value function
            mValFunction = std::move( aValFunction );

            // set setting flag
            mSetValFunction = true;
        }

        //------------------------------------------------------------------------------
        /**
         * set space derivative function
         * @param[ in ] aSpaceDerFunctions cell of functions for dnPropdxn evaluation
         */
        void set_space_der_functions(
                const Vector< PropertyFunc >& aSpaceDerFunctions );

        //------------------------------------------------------------------------------
        /**
         * get val function
         * @param[ out ] mValFunction function for property evaluation
         */
        const PropertyFunc&
        get_val_function() const
        {
            return mValFunction;
        }

        //------------------------------------------------------------------------------
        /**
         * set dof derivative functions
         * @param[ in ] aDofDerFunctions list function for property derivatives wrt dof
         */
        void
        set_dof_derivative_functions(
                const Vector< PropertyFunc >& aDofDerFunctions )
        {
            // set functions for derivatives wrt dof
            mDofDerFunctions = aDofDerFunctions;

            // set setting flag
            mSetDofDerFunctions = true;
        }

        //------------------------------------------------------------------------------
        /**
         * get dof derivative functions
         * @param[ out ] mDofDerFunctions list function for property derivatives wrt dof
         */
        const Vector< PropertyFunc >&
        get_dof_derivative_functions() const
        {
            return mDofDerFunctions;
        }

        //------------------------------------------------------------------------------
        /**
         * set dv derivative functions
         * @param[ in ] aDofDerFunctions list function for property derivatives wrt dv
         */
        void
        set_dv_derivative_functions(
                const Vector< PropertyFunc >& aDvDerFunctions )
        {
            // set functions for derivatives wrt dv
            mDvDerFunctions = aDvDerFunctions;

            // set setting flag
            mSetDvDerFunctions = true;
        }

        //------------------------------------------------------------------------------
        /**
         * get dv derivative functions
         * @param[ out ] mDofDerFunctions list function for property derivatives wrt dv
         */
        const Vector< PropertyFunc >&
        get_dv_derivative_functions() const
        {
            return mDvDerFunctions;
        }

        //------------------------------------------------------------------------------
        /**
         * set a list of active dof types
         * @param[ in ] aDofTypes list of dof types
         */
        void set_dof_type_list(
                const Vector< Vector< MSI::Dof_Type > >& aDofTypes );

        //------------------------------------------------------------------------------
        /**
         * return a list of active dof types
         * @param[ out ] mDofTypes list of dof types
         */
        const Vector< Vector< MSI::Dof_Type > >&
        get_dof_type_list() const
        {
            return mDofTypes;
        }

        //------------------------------------------------------------------------------
        /**
         * build a dof type map
         */
        void build_dof_type_map();

        //------------------------------------------------------------------------------
        /**
         * get a dof type map
         * @param[ out ] mDofTypeMap map of the dof types
         */
        const moris::Matrix< DDSMat >&
        get_dof_type_map() const
        {
            return mDofTypeMap;
        }

        //------------------------------------------------------------------------------
        /**
         * check if the property depends on a particular group of dof type
         * @param[ in ]  aDofType cell of dof type
         * @param[ out ] aBool    boolean, true if dependency on the dof type
         */
        bool check_dof_dependency( const Vector< MSI::Dof_Type >& aDofType );

        //------------------------------------------------------------------------------
        /**
         * set a list of dv types
         * @param[ in ] aDvTypes list of dv type
         */
        void set_dv_type_list( const Vector< Vector< gen::PDV_Type > >& aDvTypes );

        //------------------------------------------------------------------------------
        /**
         * return a list of dv types
         * @param[ out ] mDvTypes list of dv type
         */
        const Vector< Vector< gen::PDV_Type > >&
        get_dv_type_list() const
        {
            return mDvTypes;
        }

        //------------------------------------------------------------------------------
        /**
         * build a dv type map
         */
        void build_dv_type_map();

        //------------------------------------------------------------------------------
        /**
         * get a dv type map
         * @param[ out ] mDvTypeMap map of the dv types
         */
        const moris::Matrix< DDSMat >&
        get_dv_type_map() const
        {
            return mDvTypeMap;
        }

        //------------------------------------------------------------------------------
        /**
         * check if the property depends on a particular group of dv types
         * @param[ in ]  aDvType cell of dv type
         * @param[ out ] aBool   boolean, true if dependency on the dv type
         */
        bool check_dv_dependency( const Vector< gen::PDV_Type >& aDvType );

        //------------------------------------------------------------------------------
        /**
         * set a list of field types
         * @param[ in ] aFieldTypes list of dv type
         */
        void set_field_type_list( const Vector< Vector< mtk::Field_Type > >& aFieldTypes );

        //------------------------------------------------------------------------------
        /**
         * return a list of field types
         * @param[ out ] mFieldTypes list of field type
         */
        const Vector< Vector< mtk::Field_Type > >&
        get_field_type_list() const
        {
            return mFieldTypes;
        }

        //------------------------------------------------------------------------------
        /**
         * build a field type map
         */
        void build_field_type_map();

        //------------------------------------------------------------------------------
        /**
         * get a field type map
         * @param[ out ] mFieldTypeMap map of the field types
         */
        const moris::Matrix< DDSMat >&
        get_field_type_map() const
        {
            return mFieldTypeMap;
        }

        //------------------------------------------------------------------------------
        /**
         * check if the property depends on a particular group of field types
         * @param[ in ]  aFieldType cell of dv type
         * @param[ out ] aBool   boolean, true if dependency on the field type
         */
        bool check_field_dependency( const Vector< mtk::Field_Type >& aFieldType );

        //------------------------------------------------------------------------------
        /**
         * reset evaluation flags
         */
        void reset_eval_flags();

        //------------------------------------------------------------------------------
        /**
         * get the property value
         * @param[ out ] aVal matrix with property value
         */
        const Matrix< DDRMat >& val();

        //------------------------------------------------------------------------------
        /**
         * evaluate property in terms of coefficients and variables
         */
        void eval_Prop();

        //------------------------------------------------------------------------------
        /**
         * get property derivatives wrt a dof type
         * @param[ in ]  aDofType   cell of dof type
         * @param[ out ] adPropdDOF matrix with derivative wrt to the dof type
         */
        const Matrix< DDRMat >& dPropdDOF( const Vector< MSI::Dof_Type >& aDofType );

        //------------------------------------------------------------------------------
        /**
         * evaluate property derivatives wrt a dof type
         * @param[ in ] aDofType cell of dof type
         */
        void eval_dPropdDOF( const Vector< MSI::Dof_Type >& aDofType );

        //------------------------------------------------------------------------------
        /**
         * get the property derivatives wrt a design variable
         * @param[ in ]  aDvType   cell of dv type
         * @param[ out ] adPropdDV matrix with derivative wrt to the dv type
         */
        const Matrix< DDRMat >& dPropdDV( const Vector< gen::PDV_Type >& aDvType );

        //------------------------------------------------------------------------------
        /**
         * evaluate property derivatives wrt a design variable
         * @param[ in ] aDvType cell of dv type
         */
        void eval_dPropdDV( const Vector< gen::PDV_Type >& aDvType );

        //------------------------------------------------------------------------------
        /**
         * check if the property depends on space
         * @param[ in ]  aOrder   order of derivative wrt x
         * @param[ out ] aBool    boolean, true if dependency on the space order
         */
        bool check_space_dependency( const uint& aOrder );

        //------------------------------------------------------------------------------
        /**
         * get property derivative wrt x
         * @param[ in ]  aOrder   order of derivative wrt x
         * @param[ out ] mdPropdx matrix with derivative wrt x
         */
        const Matrix< DDRMat >& dnPropdxn( const uint& aOrder );

        //------------------------------------------------------------------------------
        /**
         * evaluate property derivative wrt x
         * @param[ in ]  aOrder   order of derivative wrt x
         */
        void eval_dnPropdxn( const uint& aOrder );

        //------------------------------------------------------------------------------
        /**
         * get non unique dof type list
         * @param[ in ] aDofType cell of dof type
         */
        void get_non_unique_dof_types( Vector< MSI::Dof_Type >& aDofTypes );

        //------------------------------------------------------------------------------
        /**
         * get non unique dof type list
         * @param[ in ] aDofType    cell of dof types
         * @param[ in ] aDvType     cell of dv types
         * @param[ in ] aFioeldType cell of field types
         */
        void get_non_unique_dof_dv_and_field_types(
                Vector< MSI::Dof_Type >&   aDofTypes,
                Vector< gen::PDV_Type >&   aDvTypes,
                Vector< mtk::Field_Type >& aFieldTypes );

        //------------------------------------------------------------------------------
    };

    //------------------------------------------------------------------------------
}    // namespace moris::fem

#endif /* SRC_FEM_CL_FEM_PROPERTY_HPP_ */

