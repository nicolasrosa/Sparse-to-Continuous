
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

        double minZ = pts2.mat().c(2).min() + 0.5;
        pts1.mat().c(2) -= minZ; pts2.mat().c(2) -= minZ;
    }
}

int main( int argc , char* argv[] )
{
    STD<String> prf;
//    prf.push_back( "../data/chair/train/chair_" );
//    prf.push_back( "../data/table/train/table_" );
    prf.push_back( "../data/chair/test/chair_" );
    prf.push_back( "../data/table/test/table_" );

    String num , suf = ".off";

    STDi st,fn;
//    st.push_back( 1 ); fn.push_back( 500 );
//    st.push_back( 1 ); fn.push_back( 390 );
    st.push_back( 890 ); fn.push_back( 989 );
    st.push_back( 393 ); fn.push_back( 492 );

    int nn = 0; forLOOPi( st.size() )
        nn += fn[i] - st[i] + 1;

    disp( nn );

    float gres = 0.05;
    Matd grd = MatGrid3d( -0.6 , +0.6 , gres );
    Matd bin( grd.r() );

    KDtreed kd; kd.add( grd );

    Matd ftr( nn , bin.s() );
    Matd lbl( ftr.r() , prf.size() ); lbl.setVal(0);

    int cnt = 0;
    forLOOPj( prf.size() )
    {
        disp( prf[j] );
        forLOOPii( st[j] , fn[j] + 1 )
        {
            disp( i );
            num = toString( i , 4 );

            Pts3d pts1,pts2;
            load( prf[j] + num + suf , pts1 , pts2 );

            SSTDi idx; SSTDd dst;
            kd.knnSearch( pts2.mat() , 1 , idx , dst );

            bin.setVal(0);
            forLOOPj( idx.size() )
            {
                if( idx[j].size() > 0 )
                    bin( idx[j][0] ) = 1;
            }

            ftr.r( cnt ) = bin.t();
            lbl( cnt , j ) = 1.0;
            cnt++;
        }
    }

    Mati idx = MatINCRi( ftr.r() ); idx.ShuffleRows();
    ftr.SampleRows( idx ); lbl.SampleRows( idx );

    ftr.info();
    lbl.info();

    ftr.save( "datasets/features_tst" );
    lbl.save( "datasets/labels_tst" );

    return 0;
}






