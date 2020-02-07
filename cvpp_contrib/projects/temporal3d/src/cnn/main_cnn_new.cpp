
#include "aux_cnn.h"
#include <cvpp/interfaces/cpplot.h>
#include <cvpp/objects/object_camera.h>

using namespace cvpp;

int main()
{
    int n = 50;

    String dir = "../data/kitti_drive_0015";
    Seq<String> str_img = get_files( dir + "/data/image_02/data/" , n );
    Seq<String> str_vel = get_files( dir + "/data/velodyne_points/data/" , n );
    Seq<String> str_gps = get_files( dir + "/data/oxts/data/" , n );

    SeqMatd vels = load_vel( str_vel );
    Seq<Img3c> imgs = load_img( str_img );
    Seq<Posef> poses = load_pos( str_gps );

    int r = imgs[0].r() , c = imgs[0].c();

    CPPlot draw( "Window" , ULHW( 2 * r , c ) , 2 , 1 );
    draw[0].set3Dworld();
    draw[1].set2Dimage().setResolution( imgs[0] );

    glPixelStorei( GL_UNPACK_ALIGNMENT , 1 );

    Matf K, D, R, P;
    load_cam2cam( dir , K , D , R , P );
    Matf T = load_vel2cam( dir );

    Matf iP = MatIDTYf( 4 );
    iP.cl(3) = P.clone(); iP = iP.inv();
    Matf iR = R.t() , iT = T.inv();

    Matf TRP = T * R * P;
    Matf iPRT = iP * iR * iT;

    int tex_img = draw.addTexture3U( imgs[0] );

    int buf_uv = draw.addBuffer2Df( 200000 );
    int buf_dep = draw.addBuffer3Df( 200000 );

    int buf_vel = draw.addBuffer3Df( 200000 );
    int buf_clr = draw.addBuffer3Df( 2000000 );
    int buf_xyz = draw.addBuffer3Df( 2000000 );

    int t = 0 , tm = imgs.size();
    bool change = true;

    while( draw.input() )
    {
        if( draw.keys.up   ) { if( t < tm - 1 ) t++; change = true; }
        if( draw.keys.down ) { if( t >      0 ) t--; change = true; }

        if( change )
        {
            Matf vel = poses[t].o2w( vels[t].cl(3).toFloat() );
            Matf proj = poses[t].w2o( vel ).appR1() * TRP;

            filter_img( proj , r , c );

            Matf dep = proj.c(2) , uv = proj.cl(2) / dep;
            Matf xyz = poses[t].o2w( ( uv.appR1() % dep ).appR1() * iPRT );

            disp( uv.sum() , dep.sum() , xyz.sum() );

            draw.updTexture( tex_img , imgs[t] );

            draw.updBuffer3D( buf_vel , vel );
            draw.updBuffer3D( buf_clr , color_pts( uv , imgs[t] ) );
            draw.updBuffer3D( buf_xyz , xyz.cl(3).clone() );

            draw.updBuffer2D( buf_uv , uv );
            draw.updBufferRGBjet( buf_dep , dep );

            change = false;
        }

        draw[0].clear();
        draw.psc(1,WHI).pts3D( buf_vel );
        draw.ps(4).pts3D( buf_xyz , buf_clr );

        forLOOPi( poses.size() )
            draw.psc(3,BLU).pt3D( poses[i].getPosPt() );
        draw.psc(5,RED).pt3D( poses[t].getPosPt() );

        draw[1].clear();
        draw.useTexture( tex_img );
        draw.ps(2).pts2D( buf_uv , buf_dep );


        draw.updateWindow(30);
    }


    return 0;
}
