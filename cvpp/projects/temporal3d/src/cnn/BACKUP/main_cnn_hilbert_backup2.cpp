
#include "aux_cnn.h"
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
    int tex_img,tex_dsp,buf_vel,buf_xyz,buf_clr,buf_uv,buf_dep;
    int buf_srf , buf_tex;

    int t = 0 , tm = imgs.size();
    bool change = true;

    Img1f dspf( imgs[t].dims() );
    Img1c dspc;

    CPPlot draw( "Window" , ULHW( 2 * r , c ) , 2 , 1 );
    draw[0].set3Dworld().setViewer(  35.7973 , -13.4488 , 19.9121 ,
                                     35.2661 , -12.6142 , 19.7661 );
    draw[1].set2Dimage().setResolution( imgs[0] );
    glPixelStorei( GL_UNPACK_ALIGNMENT , 1 );

    draw[0].calibrate( P.blu(3) , imgs[0].r() , imgs[0].c() , 0.1 , 10000.0 );

    tex_img = draw.addTexture3U( imgs[0] );
    tex_dsp = draw.addTexture3U( imgs[0] );

    buf_vel = draw.addBuffer3Df(  200000 );
    buf_xyz = draw.addBuffer3Df(  200000 );
    buf_clr = draw.addBuffer3Df(  200000 );
    buf_uv  = draw.addBuffer2Df(  200000 );
    buf_dep = draw.addBuffer3Df(  200000 );
    buf_srf = draw.addBuffer3Df( 5000000 );
    buf_tex = draw.addBuffer3Df( 5000000 );

//    Pts2f tmp_uv;
//    forLOOPij( imgs[t].r() , imgs[t].c() ) tmp_uv.insert( Pt2f( j , i ) );
//    tmp_uv.update();

//    Matf tmp_xyz = poses[t].o2w( ( ( tmp_uv.mat().appR1() * 50.0 ).appR1() * iPRT ).cl(3) );
//    Img1c dsph( imgs[t].dims() );
//    Pts3f tmp_pos;

    int show = 0 , surf = 0;
    while( draw.input() )
    {
        if( draw.keys.enter ) { show = ++show % 2; halt(100); }
        if( draw.keys.space ) { surf = ++surf % 2; halt(100); }

        if( draw.keys.up   ) { if( t < tm - 1 ) t++; change = true; }
        if( draw.keys.down ) { if( t >      0 ) t--; change = true; }

        if( change )
        {
            process_frame( poses , vels , t , TRP , iPRT , r , c ,
                           vel , xyz , proj , uv , dep );
            create_images( uv , dep , dspf , dspc );

            Matd pts  = xyz.toDouble(); SeqMatd P , M , S;
            quick_means( pts , P , M , S , 0.25 , 0.10 , 20 , 50 );
            HMreal hm( new HMfeatSqExp() , new HMregrLogistic() );
            hm.add( P , M , S ); hm.weights().setVal( 1 );
            Matd tst = MatGrid3d( pts.limRows( 0.5 ) , 0.2 );
            Matd occ = hm.query( tst );
            Matd srf = marching3D( tst , occ > 0 , occ , 0.2 );
            draw.updBuffer3D( buf_srf , srf.toFloat() );

            Matf clrtmp = color_pts( uv , imgs[t] );
            Matf clr = clrtmp | clrtmp | clrtmp | clrtmp | clrtmp | clrtmp | clrtmp;

            buf_tex = draw.addBufferRGBknn( srf.toFloat() , xyz , clr );

///////////////////////

//            {
//                dsph.setVal( 0 );
//                Matf ctr = poses[t].getPos().t();
//                Matf tmp( 841 , 3 );
//                tmp_pos.reset();

//                forLOOPi( tmp_xyz.r() )
//                {
//                    disp( i , tmp_xyz.r() );

//                    int cnt = 0;
//                    for( double u = 0.05 ; u < 1.0 ; u += 0.02 )
//                        tmp.row( cnt++ ) = ctr.eig() + u * ( tmp_xyz.eig() - ctr.eig() );
//                    Matd occ = hm.query3( tmp.toDouble() );

//                    forLOOPj( occ.r() )
//                    {
//                        if( occ(j) > 0.5 )
//                        {
//                            tmp_pos.insert( Pt3f( tmp(i,0) , tmp(i,1) , tmp(i,2) ) );
//                            break;
//                        }
//                    }
//                }
//            }

//                    draw.lwc(3,YEL).line3D( poses[t].getPos().t() | tmp_xyz.r(i) );

//                    disp( i , j );

//                    tmp_uv.eig() << j , i ;
//                    for( double u = 1.0 ; u < 85.0 ; u += 0.1 )
//                    {
////                        tmp_occ.print();

////                        if( tmp_occ(0) > 0.8 )
////                        {
////                            dsph(j,i) = int( 3.0 * tmp_occ(0) );
////                            break;
////                        }
//                    }
//                }
//            }

///////////////////////

            draw.updTexture( tex_img , imgs[t] );
            draw.updTexture( tex_dsp , dspc );

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

//        draw.psc(1,WHI).pts3D( buf_vel );
        if( surf == 0 )
            draw.ps(4).pts3D( buf_xyz , buf_clr );
        if( surf == 1 )
            draw.surf3D( buf_srf , buf_tex );

//        forLOOPi( tmp_xyz.r() )
//            draw.lwc(3,YEL).line3D( poses[t].getPos().t() | tmp_xyz.r(i) );
//        draw.psc(5,YEL).pts3D( tmp_pos );

        forLOOPi( poses.size() )
            draw.psc(3,BLU).pt3D( poses[i].getPosPt() );
        draw.psc(5,RED).pt3D( poses[t].getPosPt() );

        draw[1].clear();

        switch( show )
        {
            case 0:
                draw.useTexture( tex_img );
                draw.ps(2).pts2D( buf_uv , buf_dep ); break;
            case 1:
                draw.useTexture( tex_dsp ); break;
        }

        draw.updateWindow(30);
    }


    return 0;
}
