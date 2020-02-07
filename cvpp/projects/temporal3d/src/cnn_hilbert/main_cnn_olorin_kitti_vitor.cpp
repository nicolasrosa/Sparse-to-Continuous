
#include "aux_loads.h"
#include "aux_process.h"

#include <cvpp/interfaces/cpplot.h>
#include <cvpp/objects/object_camera.h>

#include <cvpp/algorithms/clustering/quick_means.h>
#include <cvpp/algorithms/hilbert_maps/real/hm_real.h>
#include <cvpp/algorithms/hilbert_maps/real/feats/hm_real_feat_sqexp.h>
#include <cvpp/algorithms/hilbert_maps/real/regrs/hm_real_regr_logistic.h>
#include <cvpp/algorithms/marching_cubes/marching3D.h>

using namespace cvpp;

int main( int argc , char* argv[] )
{
    cvpp::randomise( 10 );

    int n = 0;

//    String path = "/media/vguizilini/ETERNIA/Datasets/kitti/raw_data";
//    String path = "../data/kitti";
//    String type = "residential";
//    String date = "2011_09_26";
//    String fldr = "drive_0019_sync";

    String path = argv[1];
    String type = argv[2];
    String date = argv[3];
    String fldr = argv[4];

    String suf = type + "_" + date + "_" + fldr + "_";
    String dir_data = path + "/data/" + type + "/" + date + "_" + fldr + "/";
    String dir_calib = path + "/calib_files/" + date + "/";

    Seq<String> str_img = get_files( dir_data + "image_02/data/" , n );
    Seq<String> str_vel = get_files( dir_data + "velodyne_points/data/" , n );
    Seq<String> str_gps = get_files( dir_data + "oxts/data/" , n );
    dir_data += "proc_kitti";

    n = str_img.size();
    prep_dirs( dir_data );

    Img3c img0 = load_img( str_img[0] );
    Seq<Posef> poses = load_pos( str_gps );
    int r = img0.r() , c = img0.c();

    Matf K, D, R, P;
    load_cam2cam( dir_calib , K , D , R , P );
    Matf T = load_vel2cam( dir_calib );

    Matf G = load_imu2vel( dir_calib );
    Posef gps; gps.setPose( G );
    forLOOPi( poses.size() ) poses[i] -= gps;

    Matf iP = MatIDTYf( 4 );
    iP.cl(3) = P.clone(); iP = iP.inv();
    Matf iR = R.t() , iT = T.inv();

    Matf TRP = T * R * P;
    Matf iPRT = iP * iR * iT;

    int tex_img, tex_dsp1, tex_dsp2, buf_vel;
    int buf_xyz, buf_xyzclr, buf_all, buf_allclr;
    int buf_srf, buf_srfclr, buf_uv, buf_dep;

    Img1f dspf1( img0.dims() );
    Img1f dspf2( img0.dims() );
    Img1c dspc1, dspc2;

//    CPPlot draw( "Window" , ULHW( 2 * r , c ) , 2 , 1 );
//    draw[0].set3Dworld().setViewer(  35.7973 , -13.4488 , 19.9121 ,
//                                     35.2661 , -12.6142 , 19.7661 );
//    draw[1].set2Dimage().setResolution( img0 );
//    glPixelStorei( GL_UNPACK_ALIGNMENT , 1 );

//    draw[0].calibrate( P.blu(3) , img0.r() , img0.c() , 0.1 , 10000.0 );

//    tex_img  = draw.addTexture3U( img0 );
//    tex_dsp1 = draw.addTexture3U( img0 );
//    tex_dsp2 = draw.addTexture3U( img0 );

//    buf_vel    = draw.addBuffer3Df(  1000000 );
//    buf_xyz    = draw.addBuffer3Df(  1000000 );
//    buf_xyzclr = draw.addBuffer3Df(  1000000 );
//    buf_all    = draw.addBuffer3Df(  5000000 );
//    buf_allclr = draw.addBuffer3Df(  5000000 );
//    buf_srf    = draw.addBuffer3Df(  5000000 );
//    buf_srfclr = draw.addBuffer3Df(  5000000 );
//    buf_uv     = draw.addBuffer3Df(  1000000 );
//    buf_dep    = draw.addBuffer3Df(  1000000 );

    int t = 0 , show = 0 , surf = 0;
    bool change = true , with_hilbert = false;

    while( t < n )
//    while( draw.input() )
    {
//        if( draw.keys.enter ) { show = ++show % 2; halt(100); }
//        if( draw.keys.space ) { surf = ++surf % 3; halt(100); }
//        if( draw.keys.h )     { with_hilbert = !with_hilbert; halt(100); }

//        if( draw.keys.up   ) { if( t < n - 1 ) t++; change = true; }
//        if( draw.keys.down ) { if( t >     0 ) t--; change = true; }

        if( change )
        {
            Timer timer;
            Matf vel, xyz, all, proj, uv1, uv2, dep1, dep2;
            Matf xyzclr, allclr, uvxyz, depxyz;

            Matd velst = load_vel( str_vel[t] );
            Img3c imgt = load_img( str_img[t] );

            process_frame( poses[t] , velst , TRP , iPRT , r , c ,
                           vel , xyz , all , proj , uv1 , dep1 );
            timer.ptick( "process" );

            create_images( uv1 , dep1 , dspf1 , dspc1 );
            timer.ptick( "create" );

            Matd pts  = all.toDouble(); SeqMatd P , M , S;
            quick_means( pts , P , M , S , 0.10 , 0.05 , 20 , 10 );
            HMreal hm( new HMfeatSqExp() , new HMregrLogistic() );
            hm.add( P , M , S ); hm.weights().setVal( 1 );
            Matd tst = MatGrid3d( pts.limRows( 0.5 ) , 0.1 );
            Matd occ = hm.query( tst );
            timer.ptick( "hilbert" );

            Matd srf = marching3D( tst , occ > 0 , occ , 0.8 );
//            draw.updBuffer3D( buf_srf , srf.toFloat() );
            timer.ptick( "marching" );

//            xyzclr = color_pts( uv1 , imgt );
//            allclr = xyzclr | xyzclr | xyzclr | xyzclr | xyzclr | xyzclr | xyzclr;
//            buf_srfclr = draw.addBufferRGBknn( srf.toFloat() , all , allclr );
//            timer.ptick( "color" );

///////////////////////

            Pts3d extra;
            extrapolate( srf , extra , 0.01 );
            timer.ptick( "extrapolate" );
            srf |= extra.mat();
            uv2 = poses[t].w2o( srf.toFloat() ).appR1() * TRP;
            timer.ptick( "project" );
            filter_img( uv2 , r , c );
            timer.ptick( "filter" );
            dep2 = uv2.c(2) , uv2 = uv2.cl(2) / dep2;
            create_images( uv2 , dep2 , dspf2 , dspc2 );
            timer.ptick( "images" );

            disp( dspf2.mat().mean() );

///////////////////////

//            draw.updTexture( tex_img  , imgt  );
//            draw.updTexture( tex_dsp1 , dspc1 );
//            draw.updTexture( tex_dsp2 , dspc2 );

//            draw.updBuffer3D( buf_vel     , vel    );
//            draw.updBuffer3D( buf_xyz     , xyz    );
//            draw.updBuffer3D( buf_xyzclr  , xyzclr );
//            draw.updBuffer3D( buf_all     , all    );
//            draw.updBuffer3D( buf_allclr  , allclr );
//            draw.updBuffer2D( buf_uv      , uv1    );
//            draw.updBufferRGBjet( buf_dep , dep1   );

//            timer.ptick( "8" );

            save_data( suf , dir_data , str_img[t] ,
                       imgt , dspc1 , dspc2 , xyz , uvxyz , all );

            change = false;
//            halt(100);
        }

//        draw[0].clear();

//        switch( surf )
//        {
//        case 0:
//            draw.ps(4).pts3D( buf_xyz , buf_xyzclr );
//            break;
//        case 1:
//            draw.ps(4).pts3D( buf_all , buf_allclr );
//            break;
//        case 2:
//            draw.surf3D( buf_srf , buf_srfclr );
//            break;
//        }

//        forLOOPi( poses.size() )
//                draw.psc(3,BLU).pt3D( poses[i].getPosPt() );
//        draw.psc(5,RED).pt3D( poses[t].getPosPt() );

//        draw[1].clear();

//        switch( show )
//        {
//        case 0:
//            draw.useTexture( tex_img );
//            draw.ps(4).pts2D( buf_uv , buf_dep );
//            break;
//        case 1:
//            draw.useTexture( with_hilbert ? tex_dsp2 : tex_dsp1 );
//            break;
//        }

//        draw.updateWindow(30);

        t++; change = true;
    }

    return 0;
}
