
#include <cvpp/interfaces/cpplot.h>

#include "hilbert/hm_map.h"
#include "auxiliar/aux_load.h"

using namespace cvpp;

int main()
{
//    int t = 0 , n = 200;
//    String file = "../proc/kitti_drive_0015";

    int t = 0 , n = 40;
    String file = "../proc/kitti_drive_0027";

//    int t = 0 , n = 100;
//    String file = "../proc/kitti_drive_0018";

    double rad = 0.15 , reach = 15.0;

    SeqMatd occs , free , grnd; Matd circ;
    loadKITTIproc( file , n , occs , free , grnd , circ , reach );
    disp( "DATASET LOADED" );

    Map main( occs , free , t , rad ) , prev , next;
    disp( "MAIN STARTED" );

    main.detectObjects( rad ); prev = main.copy();
    disp( "FIRST OBJECTS DETECTED" );

    Matd grd = MatGrid3d( -35 , +35 , -35 , +35 , -5 , +8 , 0.25 );

    CPPlot draw( "Window" , ULHW( 600 , 2 * 800 ) , 1 , 2 );
    draw[0].set3Dworld().setViewer(-26.3897 , -3.17565 , 9.92102 ,
                                   -25.5010 , -3.06472 , 9.47620 ).setBackground(WHI);
    draw[1].set3Dworld().setViewer( draw.screen(0).viewer ).setBackground(WHI);

    int buf_occs[ occs.size() ];
    int buf_free[ free.size() ];
    int buf_srf,buf_clr;

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

            disp( "PRE PROPAGATE" );
            propagateObjects( prev , next , main , 0.9 );
            disp( "POST PROPAGATE" );

            int n_max = std::max( main.O.max() , next.O.max() );

            L1.reset( n_max , 3 ); L1.setVal(-999);
            L2.reset( n_max , 3 ); L2.setVal(-999);

            forLOOPi( L1.r() )
            {
                    Mati I1 = main.O == i + 1; Mati I2 = next.O == i + 1;
                    if( I1.filled() ) L1.r(i) = Matd( main.occs.M ).sr( I1 ).meanRows();
                    if( I2.filled() ) L2.r(i) = Matd( next.occs.M ).sr( I2 ).meanRows();
            }

            create_hm( main , grd , draw , buf_srf , buf_clr );

            change = false;
            halt(100);
        }

        draw[1].clear();

        draw.psc(1,BLU).pts3D( buf_occs[t-1] );
        draw.psc(1,BLA).pts3D( buf_occs[t  ] );
        draw.psc(3,BLA).pts3D( circ );

        forLOOPi( L1.r() )
        {
            if( L1(i,0) != -999 && L2(i,0) != -999 )
                draw.lwc(5,BLA).line3D( L1.r(i) | L2.r(i) );
        }

        draw[0].clear();

        draw.psc(8,GRE).pts3D( L1 );
        draw.psc(10,RED).pts3D( L2 );

//        draw.psc(5,YEL).pts3D( Matd( prev.occs.M ) );
//        draw.psc(5,MAG).pts3D( Matd( next.occs.M ) );

//        forLOOPi( prev.O.max() )
//            draw.psc(5,RGBint(i)).pts3D( prev.getObject(i) );
//        forLOOPi( next.O.max() )
//            draw.psc(5,RGBint(i)).pts3D( next.getObject(i) );
//        forLOOPi( main.O.max() )
//            draw.psc(9,RGBint(i)).pts3D( main.getObject(i) );

        draw.wsurf3D( buf_srf , buf_clr );

        main.drawObjects( draw );
//        main.drawVelocities( draw );

        draw.psc(3,BLA).pts3D( circ );

        draw.updateWindow(30);
    }

    return 0;
}
