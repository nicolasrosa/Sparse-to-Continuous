
#include <dirent.h>
#include <sys/types.h>
#include <cvpp/containers/matrix.h>
#include <cvpp/containers/vector.h>
#include <cvpp/containers/image.h>
#include <cvpp/properties/pose.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

using namespace cvpp;

void
filter_img( Matf& X , const int& r , const int& c )
{
    Veci idx;
    forLOOPi( X.r() )
    {   if( X(i,2) > 0.0 )
            if( X(i,0) > 0.0 && X(i,0) < c * X(i,2) &&
                X(i,1) > 0.0 && X(i,1) < r * X(i,2) ) idx.insert( i );
    }; idx.update();
    X.SampleRows( idx );
}

Matf
color_pts( const Matf& uv , const Img3c& img )
{
    Matf clr( uv.r() , 3 );
    forLOOPi( clr.r() )
    {
        clr.row(i) << img( uv(i,1) , uv(i,0) , 2 ) ,
                      img( uv(i,1) , uv(i,0) , 1 ) ,
                      img( uv(i,1) , uv(i,0) , 0 ) ;
    }
    return clr / 255.0;
}

void process_frame( const Posef& post , const Matd& velst ,
                    const Matf& TRP , const Matf& iPRT ,const int& r , const int& c ,
                    Matf& vel , Matf& xyz , Matf& all , Matf& proj , Matf& uv , Matf& dep )
{
    vel = post.o2w( velst.cl(3).toFloat() );
    proj = velst.cl(3).toFloat().appR1() * TRP;
    filter_img( proj , r , c );

    dep = proj.c(2) , uv = proj.cl(2) / dep;
    dep.AddRand( 0 , +0.25 );
    Matf xyz1 = post.o2w( ( ( uv.appR1() % ( dep        ) ).appR1() * iPRT ).cl(3) );
    Matf xyz2 = post.o2w( ( ( uv.appR1() % ( dep + 0.33 ) ).appR1() * iPRT ).cl(3) );
    Matf xyz3 = post.o2w( ( ( uv.appR1() % ( dep + 0.66 ) ).appR1() * iPRT ).cl(3) );
    Matf xyz4 = post.o2w( ( ( uv.appR1() % ( dep + 1.00 ) ).appR1() * iPRT ).cl(3) );
    Matf xyz5 = post.o2w( ( ( uv.appR1() % ( dep + 1.33 ) ).appR1() * iPRT ).cl(3) );
    Matf xyz6 = post.o2w( ( ( uv.appR1() % ( dep + 1.66 ) ).appR1() * iPRT ).cl(3) );
    Matf xyz7 = post.o2w( ( ( uv.appR1() % ( dep + 2.00 ) ).appR1() * iPRT ).cl(3) );

    xyz = xyz1;
    all = xyz1 | xyz2 | xyz3 | xyz4 | xyz5 | xyz6 | xyz7;
    all.AddRand( 0.05 );
}

void create_images( const Matf& uv , const Matf& dep ,
                    Img1f& dspf , Img1c& dspc , const double& mlt = 3.0 )
{
    dspf.setVal( 0.0 );
    forLOOPj( uv.r() )
    {
        double val = mlt * dep(j);
        if( val > 255.0 ) val = 0.0;

        float& crd = dspf( std::floor( uv(j,1) ) ,
                           std::floor( uv(j,0) ) );

        if( crd == 0 ) crd = val;
        else if( crd > val ) crd = val;
    }
    dspc = dspf.toUChar();
}

double len( const Pt3d& pt )
{
    return sqrt( std::pow( pt.x , 2.0 ) +
                 std::pow( pt.y , 2.0 ) +
                 std::pow( pt.z , 2.0 ) );
}

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

void extrapolate( const Pts3d& pts , Pts3d& npts , const double& t = 0.05 )
{
    forLOOPiii( 0 , pts.n() , 3 )
    {
        Pt3d p1 = pts[i] , p2 = pts[i+1] , p3 = pts[i+2];
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

            new_points( t1 , t2 , t , npts );
        }
    }
}

void reproject( const String& str_res , const int& r , const int& c ,
                Img1c& rest , Pts2f& respts , Vecf& resdep )
{
    rest.load( str_res );
    rest.Resize( r , c );

    forLOOPij( rest.r() , rest.c() )
    {
        int val = unsigned( rest(i,j) ) / 3;
        if( val )
        {
            respts.insert( Pt2f( j , i ) );
            resdep.insert( val );
        }
    }
    respts.update();
    resdep.update();
}

void reproject( const String& str_res , const int& r , const int& c ,
                Pts2f& respts , Vecf& resdep )
{
    Img1c rest;
    reproject( str_res , r , c , rest , respts , resdep );
}

