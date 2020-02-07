
#include <cvpp/interfaces/cpplot.h>

#include "code/aux_load.h"
#include "code/aux_align.h"
#include "code/aux_propagate.h"
#include "code/aux_hilbert.h"
#include "code/class_map.h"

using namespace cvpp;

int main()
{
        int t = 0 , n = 100;
        String file = "../proc/kitti_drive_0027";
        double rad = 0.15 , reach = 15.0;

//        int t = 0 , n = 100;
//        String file = "../proc/kitti_drive_0018";
//        double rad = 0.15 , reach = 15.0;

        SeqMatd occs , free , grnd ; Matd circ ;
        loadKITTIproc( file , n , occs , free , grnd , circ , reach );
        disp( "DATASET LOADED" );

        Matd grd = MatGrid3d( -35 , +35 , -35 , +35 , -5 , +8 , 0.25 );

        Map *main , *prev , *next;
        main = new Map( occs[t] , free[t] , rad );
        prev = new Map( occs[t] , free[t] , rad );
        next = NULL;

        CPPlot draw( "Window" , ULHW( 600 , 2 * 800 ) , 1 , 2 );
        draw[0].set3Dworld().setViewer( -26.3897 , -3.17565 , 9.92102 ,
                                        -25.5010 , -3.06472 , 9.47620 ).setBackground( WHI );
        draw[1].set3Dworld().setViewer( draw.screen(0).viewer ).setBackground( WHI );

        int buf_circ = draw.addBuffer3D( circ );
        int buf_srf = draw.addBuffer3Dd( grd.r() );
        int buf_clr = draw.addBuffer3Df( grd.r() );

        bool change = true;
        while( draw.input() )
        {
            if( draw.keys.enter )
            {
                delete prev; prev = next;
                disp( "STORING OLD MAP" );

                change = true;
            }

            if( change )
            {
                Timer tt;

                disp( "###################################" , t++ , n );

                next = new Map( occs[t] , free[t] , rad );
                tt.ptick( "INITIATING NEXT MAP" );

                alignObjects( *main , *next , 2.0 );
                tt.ptick( "ALIGNING OBJECTS" );

                propagateObjects( *prev , *next , *main , 8 , 0.9 );
                tt.ptick( "PROPAGATING OBJECTS" );

                createHilbertMap( *main , grd , draw , buf_srf , buf_clr );
                tt.ptick( "CREATING HILBERT MAP" );

                change = false;
            }

            draw[0].clear();

//            draw.psc(1,BLU).pts3D( occs[t-1] );
            draw.psc(1,BLA).pts3D( occs[t] );
            draw.psc(3,BLA).pts3D( buf_circ );
//            main->drawObjects( draw );
//            next->drawObjects( draw );

            draw[1].clear();

            draw.wsurf3D( buf_srf , buf_clr );
            draw.psc(3,BLA).pts3D( buf_circ );

            main->drawBoxes( draw );
            main->drawHistory( draw );

            draw.updateWindow(30);
        }


    return 0;
}
