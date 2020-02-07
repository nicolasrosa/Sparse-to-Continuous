
#include <cvpp/interfaces/cpplot.h>

#include <cvpp/algorithms/hilbert_maps/real/hm_real.h>
#include <cvpp/algorithms/hilbert_maps/real/feats/hm_real_feat_sqexp.h>
#include <cvpp/algorithms/hilbert_maps/real/regrs/hm_real_regr_logistic.h>

#include <cvpp/algorithms/marching_cubes/marching3D.h>

using namespace cvpp;

void new_points( const Pt3d& p1 , const Pt3d& p2 ,
                 const double& t , Pts3d& pts )
{
    Pt3d d12 = p2 - p1;
    double n12 = len( d12 );

    int q12 = std::ceil( n12 / t );
    double r12 = 1.0 / (double)q12;

    forLOOPj( q12 + 1 )
        pts.insert( p1 + d12 * (double)j * r12 );
    pts.update();
}

void load( const String& str , Pts3d& pts1 , Pts3d& pts2 ,
           const int& tt = 100 , const double& ss = 1.0 )
{
    IFile file( str );
    if( file.is_open() )
    {
        String off ; file >> off;
        int n , m , e ; file >> n >> m >> e ;
        double x , y , z ;

        forLOOPi( n )
        {
            file >> x >> y >> z;
            pts1.insert( Pt3d( x , y , z ) );
        }
        pts1.update();

        Matd lim = pts1.mat().limRows();
        double s = ( lim.r(1) - lim.r(0) ).max();
        double t = s / (double)tt;

        forLOOPi( m )
        {
            int n,a,b,c;
            file >> n >> a >> b >> c ;

            Pt3d p1 = pts1[a] , p2 = pts1[b] , p3 = pts1[c];
            Pt3d d12 = p2 - p1 , d23 = p3 - p2 , d31 = p1 - p3;
            double len12 = len( d12 ) , len23 = len( d23 ) , len31 = len( d31 );

            int q12 = std::ceil( len12 / t );
            int q23 = std::ceil( len23 / t );

            double q = std::max( q12 , q23 );
            double r = 1.0 / (double)q;

            forLOOPj( q )
            {
                Pt3d t1 = p1 + (double)j * r * ( p3 - p1 );
                Pt3d t2 = p2 + (double)j * r * ( p3 - p2 );

                new_points( t1 , t2 , t , pts2 );
            }
        }

        Matd mean = pts1.mat().meanRows();
        pts1.mat() -= mean; pts2.mat() -= mean;

        if( ss > 0 )
        {
            pts1.mat() /= s; pts2.mat() /= s;
            pts1.mat() *= ss; pts2.mat() *= ss;
        }
    }
}

int main( int argc , char* argv[] )
{

/////// DATA

    String prf = "../data/chair/train/chair_";
    String num , suf = ".off";

    int i = atoi( argv[1] );
    num = toString( i , 4 );

    Pts3d pts1,pts2;
    load( prf + num + suf , pts1 , pts2 );

    pts1.mat().RotateZ( atoi( argv[2] ) );
    pts2.mat().RotateZ( atoi( argv[2] ) );

/////// GRID

    Matd Xtr = pts2.mat();

    float gres = 0.05;
    Matd Xgr = MatGrid3d( -0.6 , +0.6 , gres );
    Matd Ygr( Xgr.r() ); Ygr.setVal(0);

    KDtreed kd;
    kd.add( Xgr );

    SSTDi idx; SSTDd dst;
    kd.knnSearch( Xtr , 1 , idx , dst );
    forLOOPi( idx.size() ) Ygr[ idx[i][0] ] = 1;

/////// DRAW

    Xtr.info();
    Xgr.info();

    CPPlot draw("Window");
    draw[0].set3Dworld().setViewer(  -1.363180 , -1.211810 , 0.7920310 ,
                                     -0.680736 , -0.594176 , 0.401144 );

    unsigned buf_grd0 = draw.addBuffer3D( Xgr.sr( Ygr == 0 ) );
    unsigned buf_grd1 = draw.addBuffer3D( Xgr.sr( Ygr == 1 ) );

    int show = 0;
    while( draw.input() )
    {
        draw[0].clear();

        switch( show )
        {
        case 0:
            draw.psc(4,RED).pts3D( pts1 );
            draw.psc(2,WHI).pts3D( pts2 );
            break;
        case 1:
            draw.psc(3,RED).pts3D( buf_grd1 );
            break;
        case 2:
            draw.psc(1,BLU).pts3D( buf_grd0 );
            draw.psc(3,RED).pts3D( buf_grd1 );
            break;
        }

        draw.updateWindow(30);

        if( draw.keys.enter )
            show = ++show % 3 , halt(100);
    }

    return 0;
}






