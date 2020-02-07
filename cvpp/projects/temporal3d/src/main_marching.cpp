
#include <cvpp/interfaces/cpplot.h>
#include <cvpp/algorithms/marching_cubes/marching3D.h>

#include "hm_incremental/hm_incremental.h"
#include "hm_incremental/cluster.h"

using namespace cvpp;

Matd calcSurf( HMincr& hm , const double& res , const double& thr )
{
    SeqMatd srf( hm.means.size() );

//    #pragma omp parallel for
    forLOOPi( hm.vars.size() )
    {
        double px = 5.0 * std::sqrt( hm.vars[i](0,0) );
        double py = 5.0 * std::sqrt( hm.vars[i](1,1) );
        double pz = 5.0 * std::sqrt( hm.vars[i](2,2) );

        double xmin = std::round( ( hm.means[i](0) - px ) / res ) * res;
        double xmax = std::round( ( hm.means[i](0) + px ) / res ) * res;
        double ymin = std::round( ( hm.means[i](1) - py ) / res ) * res;
        double ymax = std::round( ( hm.means[i](1) + py ) / res ) * res;
        double zmin = std::round( ( hm.means[i](2) - pz ) / res ) * res;
        double zmax = std::round( ( hm.means[i](2) + pz ) / res ) * res;

        Matd grid = MatGrid3d( xmin , xmax , ymin , ymax , zmin , zmax , res );
        srf[i] = marching3D( grid , hm.query3( grid ) , thr );
    }

    return Matd( srf );
}

int main()
{
    Timer t;

//    Matd pts( "../data/virtual" );
    Matd pts( "../data/riegl1" ); pts = pts.cl(3).clone();


    SeqMatd P,M,S;
    cluster( pts , P , M , S , 0.2 , 0.2 , 5 );
    Matd W = MatONESd( M.size() );

    HMincr hm; hm.add( M , S , W );

    double res = 0.5;
    double thr = 0.5;

    t.ptick( "pre" );
    Matd grid = MatGrid3d( pts.limRows( 0.5 ) , res );
    Matd prob = hm.query3( grid );
    Matd srf = marching3D( grid , prob , thr );
    t.ptick( "surf" );
    Matd blk = calcSurf( hm , res , thr );
    t.ptick( "block" );

    CPPlot draw( "Window" );
    draw[0].set3Dworld().setViewer( -2.72474 , -3.86654 , 8.28832 ,
                                    -2.35697 , -3.33153 , 7.52772 );

    int buf_srf = draw.addBuffer3D( srf );
    int buf_sclr = draw.addBufferRGBjet( srf.c(2).clone() );

    int buf_blk = draw.addBuffer3D( blk );
    int buf_bclr = draw.addBufferRGBjet( blk.c(2).clone() );

    int show = 0;
    while( draw.input() )
    {
        if( draw.keys.enter ) { show = ++show % 2; halt(100); }

        draw[0].clear();

        switch( show )
        {
        case 0:
//            draw.psc(2,WHI).pts3D( pts );
//            draw.psc(5,RED).pts3D( Matd( M ) );
//            draw.psc(5,GRE).pts3D( blk );
            draw.wsurf3D( buf_srf , buf_sclr );
            break;
        case 1:
//            draw.psc(2,WHI).pts3D( pts );
            draw.wsurf3D( buf_blk , buf_bclr );
//            forLOOPi( M.size() ) draw.clr(RED).ellipse3D( M[i] , S[i] );
            break;
        }

        draw.updateWindow(30);
    }

    return 0;
}
