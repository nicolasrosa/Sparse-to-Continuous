
#include <cvpp/interfaces/cpplot.h>

#include "hilbert/hm_map.h"
#include "auxiliar/aux_load.h"

using namespace cvpp;

void ptsColor( const Matd& pts , Map& main , CPPlot& draw , int& buf_clr )
{
    KDtreed kd( Matd( main.occs.M ) ); Matd clr( pts.r() );
    SSeqi idx; SSeqd dst; kd.knnSearch( pts , 1 , idx , dst );

    for( unsigned i = 0 ; i < idx.size() ; i++ )
        if( idx[i].size() > 0 )
            clr(i) = main.T[ main.O[ idx[i][0] ] - 1 ].v.rsqsum();
    buf_clr = draw.addBufferRGBjet( clr , 0.0 , 1.0 );
}

int main()
{
    int t = 0 , n = 100;
    double rad = 0.2 , reach = 18.0;
    String file = "../proc/kitti_drive_0027";

    SeqMatd occs , free , grnd;
    loadKITTIproc( file , n , occs , free , grnd , reach );

    Matd circle( 360 * 10 , 3 );
    forLOOPi( circle.r() )
    {
        double ang = double(i) / ( 2 * PI * 10 );
        circle.row(i) << reach * std::cos(ang) , reach * std::sin(ang) , -2.0 ;
    }

    Map main( occs , free , t , rad ) , prev , next;
    main.detectObjects( rad ); prev = main.copy();

    CPPlot draw( "Window" );
    draw[0].set3Dworld().setViewer(-26.3897 , -3.17565 , 9.92102 ,
                                   -25.5010 , -3.06472 , 9.47620 ).setBackground(WHI);

    int buf_occs[ occs.size() ] , buf_clr;
    int buf_free[ free.size() ];

    forLOOPi( occs.size() )
    {
        buf_occs[i] = draw.addBuffer3D( occs[i] );
        buf_free[i] = draw.addBuffer3D( free[i] );
    }

    Matd L1,L2;
    bool has = true , change = true;
    while( draw.input() )
    {
        if( draw.keys.enter )
        {
            prev = next.copy();
            change = true;
        }

        if( change )
        {
            disp( "###################################" , t++ , n );

            if( t < n )
            {
                next = Map( occs , free , t , rad );
                next.detectObjects( rad );
                alignObjects( main , next );
            }
            else
            {
                has = false;
                next = Map();
            }

            propagateObjects( prev , next , main , 0.9 );

            int n_max = std::max( main.O.max() , next.O.max() );

            L1.reset( n_max , 3 ); L1.setVal(-999);
            L2.reset( n_max , 3 ); L2.setVal(-999);

            forLOOPi( L1.r() )
            {
                    Mati I1 = main.O == i + 1; Mati I2 = next.O == i + 1;
                    if( I1.filled() ) L1.r(i) = Matd( main.occs.M ).sr( I1 ).meanRows();
                    if( I2.filled() ) L2.r(i) = Matd( next.occs.M ).sr( I2 ).meanRows();
            }

            change = false;
            halt(100);
        }

        draw[0].clear();

        ptsColor( occs[t] , main , draw , buf_clr );
        draw.ps(2).pts3D( buf_occs[t] , buf_clr );
        draw.psc(3,BLA).pts3D( circle );
        main.drawObjects( draw );

        draw.updateWindow(30);
    }

    return 0;
}
