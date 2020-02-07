
#include <cvpp/interfaces/cpplot.h>
#include <cvpp/algorithms/clustering/quick_means.h>

#include <cvpp/algorithms/hilbert_maps/real/hm_real.h>
#include <cvpp/algorithms/hilbert_maps/real/feats/hm_real_feat_sqexp.h>
#include <cvpp/algorithms/hilbert_maps/real/regrs/hm_real_regr_logistic.h>

#include <cvpp/algorithms/clustering/kmeans.h>
#include <cvpp/algorithms/clustering/kmeans/inits/kmeans_init_ss.h>
#include <cvpp/algorithms/clustering/kmeans/opts/kmeans_opt_full.h>

#include <cvpp/algorithms/marching_cubes/marching3D.h>

#include "hm_incremental/hm_incremental.h"

using namespace cvpp;

Matd calcSurf( HMincr& hm , const double& res , const double& thr )
{
    int n = hm.means.size();
    SeqMatd srf(n) , grd(n) , prb(n);

    KDtreed tree( hm.ctrs ); SSeqi idxs; SSeqd dsts;
    tree.knnSearch( hm.ctrs , 10 , idxs , dsts );

    Timer t;

    int ngrid = 0;

    #pragma omp parallel for
    forLOOPi( n )
    {
        double px = 3.5 * std::sqrt( hm.vars[i](0,0) );
        double py = 3.5 * std::sqrt( hm.vars[i](1,1) );
        double pz = 3.5 * std::sqrt( hm.vars[i](2,2) );

        double xmin = std::round( ( hm.means[i](0) - px ) / res ) * res;
        double xmax = std::round( ( hm.means[i](0) + px ) / res ) * res;
        double ymin = std::round( ( hm.means[i](1) - py ) / res ) * res;
        double ymax = std::round( ( hm.means[i](1) + py ) / res ) * res;
        double zmin = std::round( ( hm.means[i](2) - pz ) / res ) * res;
        double zmax = std::round( ( hm.means[i](2) + pz ) / res ) * res;

        grd[i] = MatGrid3d( xmin , xmax , ymin , ymax , zmin , zmax , res );
        prb[i] = hm.query3( grd[i] , idxs[i] );

        ngrid += grd[i].r();
    }
    t.ptick( "block_querying" );

    disp( ngrid );

    forLOOPi( n )
        srf[i] = marching3D( grd[i] , prb[i] , thr );
    t.ptick( "block_marching" );

    return Matd( srf );
}

int main( int argc , char* argv[] )
{
    Timer t;

    int thr = 5;
    double res = 0.05 , rad = 0.02 ;
    double reach = 25.0;

    Matd all( "../data/riegl1" );
    Matd Xtr = all.cl(3).clone().AddRand( 0.05 );
    Matd clr = all.cr(1).clone();

    Mati idx;

    idx = Xtr.c(0) > -reach; Xtr.SampleRows( idx ); clr.SampleRows( idx );
    idx = Xtr.c(0) < +reach; Xtr.SampleRows( idx ); clr.SampleRows( idx );
    idx = Xtr.c(1) > -reach; Xtr.SampleRows( idx ); clr.SampleRows( idx );
    idx = Xtr.c(1) < +reach; Xtr.SampleRows( idx ); clr.SampleRows( idx );
    idx = Xtr.c(2) > -reach; Xtr.SampleRows( idx ); clr.SampleRows( idx );
    idx = Xtr.c(2) < +reach; Xtr.SampleRows( idx ); clr.SampleRows( idx );

    Xtr.SampleRows( 5 );
    clr.SampleRows( 5 );

    t.ptick( "loading" );

    SeqMatd P,M,S;
//    quick_means( Xtr , P , M , S , rad , rad , thr );
    quick_means( Xtr , P , M , S , rad , rad , 10 , thr );
    t.ptick( "clustering" );

    HMreal hm( new HMfeatSqExp() , new HMregrLogistic() );
    hm.add( P , M , S ); hm.weights().setVal(1);
    t.ptick( "training" );

    t.tick();
    Matd Xte = MatGrid3d( Xtr.limRows( 0.5 ) , res ); Xte.info( "grid1" );
    Matd Yte = hm.query( Xte );
    t.ptick( "querying" );

    t.tick();
    Matd srf = marching3D( Xte , Yte > 0 , Yte , 0.5 );
    t.ptick( "marching" );


//    HMincr hm2;
//    hm2.add( M , S , MatONESd( M.size() ) );
//    Matd blk = calcSurf( hm2 , res , 0.5 );
//    t.ptick( "block" );

//    KDtreed tree( Xtr ); SSeqi idxs; SSeqd dsts;
//    tree.knnSearch( srf , 1 , idxs , dsts );
//    Matd srf_sclr = MatZEROSd( srf.r() );
//    forLOOPi( srf.r() ) if( idxs[i].size() > 0 ) srf_sclr(i) = clr( idxs[i][0] );
//    tree.knnSearch( blk , 1 , idxs , dsts );
//    Matd srf_bclr = MatZEROSd( blk.r() );
//    forLOOPi( blk.r() ) if( idxs[i].size() > 0 ) srf_bclr(i) = clr( idxs[i][0] );

    CPPlot draw("Window", ULHW( 960 , 1280 ));
    draw.useScreen(0).set3Dworld().setViewer( -2.92440 , 6.27890 , 6.54793 ,
                                              -2.43251 , 5.77697 , 5.83651 ).setBackground( WHI );

    unsigned buf_pts = draw.addBuffer3D( Xtr );
    unsigned buf_srf = draw.addBuffer3D( srf );
//    unsigned buf_blk = draw.addBuffer3D( blk );
    unsigned buf_ctr = draw.addBuffer3D( hm.ctrs() );
    unsigned buf_sclr = draw.addBufferRGBjet( srf.c(2).clone() );//srf_sclr , 0.0 , 1.0 );
//    unsigned buf_bclr = draw.addBufferRGBjet( blk.c(2).clone() );//srf_bclr , 0.0 , 1.0 );
    unsigned buf_pts_clr = draw.addBufferRGBjet( clr , 0.0 , 1.0 );

    // DRAW

    int show = 0;
    while( draw.input() )
    {
        if( draw.keys.enter ) { show = ++show % 2; halt(100); }

        draw[0].clear();

        switch( show )
        {
        case 0:
            draw.color(BLU).wsurf3D( buf_srf , buf_sclr );
            break;
        case 1:
//            draw.color(BLU).wsurf3D( buf_blk , buf_bclr );
            break;
        }

        draw.updateWindow(30);
    }

    return 0;
}
