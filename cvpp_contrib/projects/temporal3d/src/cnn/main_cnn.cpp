
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
    n = str_img.size();

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

    Matf G = load_imu2vel( dir );
    Posef gps; gps.setPose( G );
    forLOOPi( poses.size() ) poses[i] -= gps;

    Matf iP = MatIDTYf( 4 );
    iP.cl(3) = P.clone(); iP = iP.inv();
    Matf iR = R.t() , iT = T.inv();

    Matf TRP = T * R * P;
    Matf iPRT = iP * iR * iT;

    SeqMatf vel(n),xyz(n),proj(n),uv(n),dep(n);
    Veci tex_img(n),tex_dsp(n),buf_vel(n),buf_xyz(n),buf_clr(n),buf_uv(n),buf_dep(n);

    #pragma omp parallel for
    forLOOPi( n )
    {
        process_frame( poses , vels , i , TRP , iPRT , r , c ,
                       vel[i] , xyz[i] , proj[i] , uv[i] , dep[i] );
    }

    forLOOPi( n )
    {
        tex_img[i] = draw.addTexture( imgs[i] );

        Img1f dspf( imgs[i].dims() ); dspf.setVal( 0.0 );
        forLOOPj( uv[i].r() )
            dspf( std::floor( uv[i](j,1) + 0.5 ) ,
                  std::floor( uv[i](j,0) + 0.5 ) ) = dep[i](j);
        Img1c dspc = ( dspf ).toUChar();
        tex_dsp[i] = draw.addTexture( dspc );

        buf_vel[i] = draw.addBuffer3D( vel[i] );
        buf_xyz[i] = draw.addBuffer3D( xyz[i] );
        buf_clr[i] = draw.addBuffer3D( color_pts( uv[i] , imgs[i] ) );
        buf_uv[i]  = draw.addBuffer2D( uv[i] );
        buf_dep[i] = draw.addBufferRGBjet( dep[i] );
    }

    int t = 0 , tm = imgs.size();
    bool change = true;

    int show = 0;
    while( draw.input() )
    {
        if( draw.keys.enter ) { show = ++show % 2; halt(100); }

        if( draw.keys.up   ) { if( t < tm - 1 ) t++; change = true; }
        if( draw.keys.down ) { if( t >      0 ) t--; change = true; }

        if( change )
        {
            change = false;
            halt(10);
        }

        draw[0].clear();
//        draw.psc(1,WHI).pts3D( buf_vel );

        forLOOPi( std::min( t + 1 , n ) )
            draw.ps(4).pts3D( buf_xyz[i] , buf_clr[i] );

//        forLOOPi( std::min( t + 1 , n ) )
//            draw.psc(2,WHI).pts3D( buf_vel[i] );

        forLOOPi( poses.size() )
            draw.psc(3,BLU).pt3D( poses[i].getPosPt() );
        draw.psc(5,RED).pt3D( poses[t].getPosPt() );

        draw[1].clear();

        switch( show )
        {
            case 0:
                draw.useTexture( tex_img[t] );
                draw.ps(2).pts2D( buf_uv[t] , buf_dep[t] ); break;
            case 1:
                draw.useTexture( tex_dsp[t] ); break;
        }

        draw.updateWindow(30);
    }


    return 0;
}
