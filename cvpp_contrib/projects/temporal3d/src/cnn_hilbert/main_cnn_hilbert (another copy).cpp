
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

int main()
{
    int n = 0;

//    String path = "/media/vguizilini/ETERNIA/Datasets/kitti/raw_data";
    String path = "../data/kitti";

    String type = "residential";
    String date = "2011_09_26";
    String fldr = "drive_0019_sync";

//    String path = argv[1];
//    String type = argv[2];
//    String date = argv[3];
//    String fldr = argv[4];

    String suf = type + "_" + date + "_" + fldr + "_";
    String dir_data = path + "/data/" + type + "/" + date + "_" + fldr + "/";
    String dir_calib = path + "/calib_files/" + date + "/";

    Seq<String> str_img = get_files( dir_data + "image_02/data/" , n );
    Seq<String> str_vel = get_files( dir_data + "velodyne_points/data/" , n );
    Seq<String> str_gps = get_files( dir_data + "oxts/data/" , n );

    n = str_img.size();

    prep_dirs( dir_data );

    SeqMatd vels = load_vel( str_vel );
    Seq<Img3c> imgs = load_img( str_img );
    Seq<Posef> poses = load_pos( str_gps );

    int r = imgs[0].r() , c = imgs[0].c();

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

    Matf vel,xyz,proj,uv,dep;
    int tex_img,tex_dsp1,tex_dsp2,buf_vel,buf_xyz,buf_clr,buf_uv,buf_dep;
    int buf_srf , buf_tex;

    int t = 0 , tm = imgs.size();
    bool change = true;

    Img1f dspf1( imgs[t].dims() );
    Img1c dspc1;

    Img1f dspf2( imgs[t].dims() );
    Img1c dspc2;

    Matd all_xyzfar;

    CPPlot draw( "Window" , ULHW( 2 * r , c ) , 2 , 1 );
    draw[0].set3Dworld().setViewer(  35.7973 , -13.4488 , 19.9121 ,
                                     35.2661 , -12.6142 , 19.7661 );
    draw[1].set2Dimage().setResolution( imgs[0] );
    glPixelStorei( GL_UNPACK_ALIGNMENT , 1 );

    draw[0].calibrate( P.blu(3) , imgs[0].r() , imgs[0].c() , 0.1 , 10000.0 );

    tex_img = draw.addTexture3U( imgs[0] );
    tex_dsp1 = draw.addTexture3U( imgs[0] );
    tex_dsp2 = draw.addTexture3U( imgs[0] );

    buf_vel = draw.addBuffer3Df( 2000000 );
    buf_xyz = draw.addBuffer3Df( 2000000 );
    buf_clr = draw.addBuffer3Df( 2000000 );
    buf_uv  = draw.addBuffer2Df( 5000000 );
    buf_dep = draw.addBuffer3Df( 5000000 );
    buf_srf = draw.addBuffer3Df( 5000000 );
    buf_tex = draw.addBuffer3Df( 5000000 );

    bool with_hilbert = false;
    int show = 0 , surf = 0;
    while( draw.input() )
    {
        if( draw.keys.enter ) { show = ++show % 2; halt(100); }
        if( draw.keys.space ) { surf = ++surf % 2; halt(100); }
        if( draw.keys.h )     { with_hilbert = !with_hilbert; halt(100); }

        if( draw.keys.up   ) { if( t < tm - 1 ) t++; change = true; }
        if( draw.keys.down ) { if( t >      0 ) t--; change = true; }

        if( change )
        {
            process_frame( poses , vels , t , TRP , iPRT , r , c ,
                           vel , xyz , proj , uv , dep );
            create_images( uv , dep , dspf1 , dspc1 );

            Matd pts  = xyz.toDouble(); SeqMatd P , M , S;
            quick_means( pts , P , M , S , 0.10 , 0.10 , 20 , 30 );
            HMreal hm( new HMfeatSqExp() , new HMregrLogistic() );
            hm.add( P , M , S ); hm.weights().setVal( 1 );
            Matd tst = MatGrid3d( pts.limRows( 0.5 ) , 0.2 );
            Matd occ = hm.query( tst );

            Matd srf = marching3D( tst , occ > 0 , occ , 0.8 );
            draw.updBuffer3D( buf_srf , srf.toFloat() );

            Matf clrtmp = color_pts( uv , imgs[t] );
            Matf clr = clrtmp | clrtmp | clrtmp | clrtmp | clrtmp | clrtmp | clrtmp;
            buf_tex = draw.addBufferRGBknn( srf.toFloat() , xyz , clr );

///////////////////////

            uv = poses[t].w2o( srf.toFloat() ).appR1() * TRP;
            filter_img( uv , r , c );
            dep = uv.c(2) , uv = uv.cl(2) / dep;
            create_images( uv , dep , dspf2 , dspc2 );

/////////////////////

//            dspf2.setVal( 0 );
//            Matd xyzfar , ctr = poses[t].getPos().t().toDouble();
//            Matf tmp_uv( 1 , 4 );

//            int w = dspf2.c() - 1;
//            int h = dspf2.r() - 1;
//            int wdw = 2;

//            forLOOPk( uv.r() )
//            {
//                disp( k , uv.r() );

//                int mini = std::max( 0 , int( uv(k,0) - wdw ) );
//                int maxi = std::min( w , int( uv(k,0) + wdw ) );
//                int minj = std::max( 0 , int( uv(k,1) - wdw ) );
//                int maxj = std::min( h , int( uv(k,1) + wdw ) );

//                Matd dst, occ, diff;
//                Matf tmp;

//                for( int i = mini ; i <= maxi ; i++ )
//                {   for( int j = minj ; j <= maxj ; j++ )
//                    {
//                        if( dspf2(j,i) == 0 )
//                        {
//                            tmp_uv.eig() << i , j , 1.0 , 1.0 ;
//                            xyzfar = poses[t].o2w( ( tmp_uv * iPRT ).cl(3) ).toDouble();
//                            diff = xyzfar - ctr;

//                            dspf2(j,i) = -1;
//                            for( double u = 1.0 ; u < 85.0 ; u += 1.0 )
//                            {
//                                dst = ctr + u * diff;
//                                occ = hm.query3( dst );

//                                if( occ(0) > 0.8 )
//                                {
//                                    dspf2(j,i) = 3.0 * dst.rsqsum();
//                                    break;
//                                }
//                            }
//                        }
//                    }
//                }
//            }

//            forLOOPij( dspf2.r() , dspf2.c() )
//            {
//                    if( dspf2(i,j) == -1 ) dspf2(i,j) = 0;
//                    if( dspf2(i,j) > 255 ) dspf2(i,j) = 0;
//            }
//            dspc2 = dspf2.toUChar();

///////////////////////

            draw.updTexture( tex_img , imgs[t] );
            draw.updTexture( tex_dsp1 , dspc1 );
            draw.updTexture( tex_dsp2 , dspc2 );

            draw.updBuffer3D( buf_vel , vel );
            draw.updBuffer3D( buf_xyz , xyz );
            draw.updBuffer3D( buf_clr , clr );
            draw.updBuffer2D( buf_uv , uv );
            draw.updBufferRGBjet( buf_dep , dep );

//            save_data( dir , str_img[t] ,
//                       dspf , dspc , vel , xyz , uv , dep );

            change = false;
            halt(10);
        }

        draw[0].clear();

        if( surf == 0 )
            draw.ps(4).pts3D( buf_xyz , buf_clr );
        if( surf == 1 )
            draw.surf3D( buf_srf , buf_tex );

        forLOOPi( poses.size() )
            draw.psc(3,BLU).pt3D( poses[i].getPosPt() );
        draw.psc(5,RED).pt3D( poses[t].getPosPt() );

        draw[1].clear();

        switch( show )
        {
            case 0:
                draw.useTexture( tex_img );
                draw.ps(4).pts2D( buf_uv , buf_dep );
                break;
            case 1:
                draw.useTexture( with_hilbert ? tex_dsp2 : tex_dsp1 );
                break;
        }

        draw.updateWindow(30);
    }


    return 0;
}
