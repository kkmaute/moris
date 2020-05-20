/*
 * cl_GEN_Field.hpp
 *
 *  Created on: Feb 26, 2020
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_GEN_MAIN_SRC_GEOMENG_CL_GEN_FIELD_HPP_
#define PROJECTS_GEN_GEN_MAIN_SRC_GEOMENG_CL_GEN_FIELD_HPP_

#include "cl_Matrix.hpp"

#include "cl_GEN_Property.hpp"

#include "cl_MTK_Mesh_Manager.hpp"

namespace moris
{
namespace ge
{

class GEN_Field
{
private :
    moris::Matrix< DDRMat >                 mFieldData;

    std::shared_ptr< Property >         mProperty = nullptr;

    moris::Cell<moris::mtk::Vertex const *> mAllVertices;

//------------------------------------------------------------------------------
public :
    // condtuctor to set the property and initialize right away
    GEN_Field( std::shared_ptr< Property > aProperty,
               mtk::Mesh_Manager*              aMeshManager )   :
        mProperty(aProperty)
    {
        this->initialize( aMeshManager );
    }
//------------------------------------------------------------------------------
    // constructor to only set the property pointer ( if using this, will need to call initialize() separately )
    GEN_Field( std::shared_ptr< Property > aProperty )   :
        mProperty(aProperty)
    {

    }
//------------------------------------------------------------------------------
    // empty constructor ( if using this, will need to call set_field_data() and/or initialize() separately )
    GEN_Field(){}
//------------------------------------------------------------------------------
    // destructor
    ~GEN_Field(){}
//------------------------------------------------------------------------------
    void initialize( mtk::Mesh_Manager* aMeshManager,
                     const uint         aBSplineMeshIndex = 0,
                     const moris_index  aWhichPair = 0 )
    {
        // get all the vertices from the interpolation mesh (including those on the aura)
        mAllVertices = aMeshManager->get_interpolation_mesh( aWhichPair )->get_all_vertices();

        uint tNumVertices = mAllVertices.size();

        moris::Matrix< DDRMat > tCoordinates;     // the property needs to know all the coordinates to perform the evaluation

        Matrix< DDRMat > tAllVals(tNumVertices, 1);
        for(uint i=0; i<tNumVertices; i++)  // fill the parameter cell with all the node coordinates
        {
            tAllVals(i) = mProperty->evaluate_field_value(mAllVertices(i)->get_id(), mAllVertices(i)->get_coords());
        }

        // check to see that the provided property is returning a row vector with values at all the nodes
        MORIS_ASSERT( tAllVals.n_rows() == tNumVertices, "GEN_Field::initialize() - the property function provided to the field needs to return a row vector with data at all nodes" );

        // set the size of the member data vector accordingly
        mFieldData.set_size( tNumVertices,1, std::numeric_limits<moris::real>::max() );    // fill with REAL_MAX values to use as a check in ASSERT cases

        for( uint iVert=0; iVert<tNumVertices; iVert++ )
        {
            if( mAllVertices( iVert )->get_interpolation( aBSplineMeshIndex )->get_weights()->numel() > 0 )    // check if the node is on the aura
            {
                MORIS_ASSERT( mProperty != nullptr, "GEN_Geom_Field::initialize() - property function not set" );    // check that the property is set

                mFieldData(iVert) = tAllVals(iVert);    // only add field values for nodes which are NOT on the aura
            }
        }
    }
//------------------------------------------------------------------------------
    bool is_analytic()
    {
        return true;
    }
//------------------------------------------------------------------------------
    void set_field_property( std::shared_ptr< Property > aProperty )
    {
        mProperty = aProperty;
    }
//------------------------------------------------------------------------------
    void set_field_data( moris::Matrix< DDRMat > aFieldData )
    {
        mFieldData = aFieldData;
    }
//------------------------------------------------------------------------------
    void set_field_data_value( const uint        aValueIndex,
                               const moris::real aValue )
    {
        mFieldData( aValueIndex ) = aValue;
    }
//------------------------------------------------------------------------------
    moris::Matrix< DDRMat > get_all_field_data()
    {
        MORIS_ASSERT( mFieldData.numel() > 0, "GEN_Field::get_field_data() - field data not set" );
        return mFieldData;
    }
//------------------------------------------------------------------------------
    moris::real get_field_val_at_vertex( const moris_index aIndex )
    {
        MORIS_ASSERT( mFieldData(aIndex) != std::numeric_limits<moris::real>::max(),"GEN_Field::get_field_val_at_vertex() - requested vertex is on the aura" );

        return mFieldData( aIndex );
    }
//------------------------------------------------------------------------------
    /**
     * returns the T-Matrix of a specific vertex
     */
    const Matrix< DDRMat > * get_t_matrix( const moris_index aVertexIndex,
                                           const uint aBSplineMeshIndex = 0 ) const
    {
        return mAllVertices( aVertexIndex )->get_interpolation( aBSplineMeshIndex )->get_weights();
    }
//------------------------------------------------------------------------------

};

}   /* end ge namespace */
}   /* end moris namespace */



#endif /* PROJECTS_GEN_GEN_MAIN_SRC_GEOMENG_CL_GEN_FIELD_HPP_ */
