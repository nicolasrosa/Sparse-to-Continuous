
#include <cvpp/interfaces/cpplot.h>

#include <cvpp/algorithms/hilbert_maps/real/hm_real.h>
#include <cvpp/algorithms/hilbert_maps/real/feats/hm_real_feat_sqexp.h>
#include <cvpp/algorithms/hilbert_maps/real/regrs/hm_real_regr_logistic.h>

#include <cvpp/algorithms/marching_cubes/marching3D.h>

using namespace cvpp;

int main( int argc , char* argv[] )
{

/////// GRID

//    Matd ftr( "datasets/features_trn" );
//    Matd lbl( "datasets/labels_trn" );
    Matd ftr( "datasets/features_tst" );
    Matd lbl( "datasets/labels_tst" );

    ftr.info();

    float gres = 0.05;
    Matd Xgr = MatGrid3d( -0.6 , +0.6 , gres );
    Matd Ygr = ftr.r(0).t();

/////// DRAW

    CPPlot draw("Window");
    draw[0].set3Dworld().setViewer(  -1.363180 , -1.211810 , 0.7920310 ,
                                     -0.680736 , -0.594176 , 0.401144 );

    unsigned buf_grd0 = draw.addBuffer3D( Xgr.sr( Ygr == 0 ) );
    unsigned buf_grd1 = draw.addBuffer3D( Xgr.sr( Ygr == 1 ) );

    bool change = false;

    int show = 0 , cnt = 0;
    disp( cnt , lbl.r(cnt).maxColsIDX() );

    while( draw.input() )
    {
        draw[0].clear();

        switch( show )
        {
        case 0:
            draw.psc(3,RED).pts3D( buf_grd1 );
            break;
        case 1:
            draw.psc(1,BLU).pts3D( buf_grd0 );
            draw.psc(3,RED).pts3D( buf_grd1 );
            break;
        }

        draw.updateWindow(30);

        if( draw.keys.enter )
            show = ++show % 2 , halt(100);

        if( draw.keys.right && cnt < ftr.r() - 1 ) { cnt++; change = true; }
        if( draw.keys.left  && cnt > 0           ) { cnt--; change = true; }

        if( change )
        {
            disp( cnt , lbl.r(cnt).maxColsIDX() );
            Ygr = ftr.r( cnt ).t();

            buf_grd0 = draw.addBuffer3D( Xgr.sr( Ygr == 0 ) );
            buf_grd1 = draw.addBuffer3D( Xgr.sr( Ygr == 1 ) );

            change = false;
            halt(100);
        }
    }

    return 0;
}






