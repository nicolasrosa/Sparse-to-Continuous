#ifndef AUX_HILBERT_H
#define AUX_HILBERT_H

#include <cvpp/interfaces/cpplot.h>

#include <cvpp/algorithms/hilbert_maps/real/hm_real.h>
#include <cvpp/algorithms/hilbert_maps/real/feats/hm_real_feat_sqexp.h>
#include <cvpp/algorithms/hilbert_maps/real/regrs/hm_real_regr_logistic.h>
#include <cvpp/algorithms/marching_cubes/marching3D.h>

#include "class_map.h"
#include "class_group.h"

using namespace cvpp;

/// Create Hilbert Map

void createHilbertMap( Map& main , Matd& grd , CPPlot& draw ,
                       int& buf_srf , int& buf_clr )
{
    HMreal hm( new HMfeatSqExp() , new HMregrLogistic() );

    SeqMatd P , M , S; Matd W;
    main.getData( P , M , S , W ); hm.add( P , M , S ); hm.weights() = W;
    Matd prb = hm.query( grd ); Matd srf = marching3D( grd , prb > 0 , prb , 0.7 );

    Matd MM = M; KDtreed tree( MM ); Matd clr( srf.r() );
    SSeqi idx; SSeqd dst; tree.knnSearch( srf , 1 , idx , dst );
    for( unsigned i = 0 ; i < idx.size() ; i++ )
        if( idx[i].size() > 0 ) clr(i) = W( idx[i][0] );

    draw.updBuffer3D( buf_srf , srf );
    draw.updBufferRGBjet( buf_clr , clr , 0.0 , 1.0 );
}

#endif
