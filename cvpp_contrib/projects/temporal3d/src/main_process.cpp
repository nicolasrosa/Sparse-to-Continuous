
#include <cvpp/interfaces/cpplot.h>
#include "auxiliar/aux_process.h"

using namespace cvpp;

int main()
{
//    int n = 113;
//    float ground = -1.5;
//    String file_input = "../data/kitti_drive_0027/data";
//    String file_output = "../proc/kitti_drive_0027";

//    int n = 220; float ground = -1.7;
//    String file_input = "../data/kitti_drive_0015/data";
//    String file_output = "../proc/kitti_drive_0015";

    int n = 113; float ground = -1.5;
    String file_input = "../data/kitti_drive_0017/data";
    String file_output = "../proc/kitti_drive_0017";

//    int n = 269; float ground = -1.5;
//    String file_input = "../data/kitti_drive_0018/data";
//    String file_output = "../proc/kitti_drive_0018";

    Matd pos = loadKITTIpos( file_input , n );
    SeqMatd vel = loadKITTIvel( file_input , n );

    SeqPosed pose( pos.r() );
    SeqMatd xyz( pos.r() ) , grd( pos.r() );

    forLOOPi( pose.size() )
    {
        pose[i].setPose( pos.r(i) );
        vel[i] = pose[i].o2w( vel[i].cl(3) );
        xyz[i] = vel[i].sr( vel[i].c(2)  > ground );
        grd[i] = vel[i].sr( vel[i].c(2) <= ground );
    }

    forLOOPi( pose.size() )
        xyz[i] *= pose[0].R().t();
    pos.cl(3) *= pose[0].R().t();

    CPPlot draw( "Window" );
    draw[0].set3Dworld().setBackground(WHI);

    int buf_xyz[ vel.size() ];
    forLOOPi( vel.size() )
        buf_xyz[i] = draw.addBuffer3D( xyz[i] );

    Pts3d base;
    base.push( Pt3d( 0 , 0 , 0 ) );
    base.push( Pt3d( 5 , 0 , 0 ) );
    base.push( Pt3d( 0 , 5 , 0 ) );
    base.push( Pt3d( 0 , 0 , 5 ) );

    int t = 0;
    while( draw.input() )
    {
        if( draw.keys.up   ) { t++; if( t == vel.size() ) t--; disp(t); halt(100); }
        if( draw.keys.down ) { t--; if( t == -1         ) t++; disp(t); halt(100); }

        draw[0].clear();

        draw.psc(2,BLA).pts3D( buf_xyz[t] );

        Pts3d axis = pose[t].o2w( base );
        axis.mat() *= pose[0].R().t();

        draw.psc(8,RED).line3D( axis[0] , axis[1] );
        draw.psc(8,YEL).line3D( axis[0] , axis[2] );
        draw.psc(8,GRE).line3D( axis[0] , axis[3] );

        draw.psc(8,RED).pts3D( pos.cl(3).clone() );

        draw.updateWindow(30);
    }

    forLOOPi( xyz.size() )
    {
        xyz[i].saveBIN( file_output + "/pts/" + toString( i , 10 ) + ".bin" );
        grd[i].saveBIN( file_output + "/grd/" + toString( i , 10 ) + ".bin" );
    }
    pos.saveBIN( file_output + "/poses.bin" );

    return 0;
}
