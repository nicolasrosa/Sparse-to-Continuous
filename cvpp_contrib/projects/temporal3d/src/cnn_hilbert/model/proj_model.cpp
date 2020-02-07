
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

void load( const String& str , Pts3d& pts1 , Pts3d& pts2 , const double& t )
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
    }
}

int main()
{
////// DATASET

    Pts3d pts1,pts2;
    load( "../data/tree_base.off" , pts1 , pts2 , 0.02 );
    pts1.mat().SpinX(-90).AddRand( 0.01 );
    pts2.mat().SpinX(-90).AddRand( 0.01 );

    double hres = 0.02;
    double rad = 0.01; int thr = 5;

    Matd Xtr = pts2.mat().clone();
    Matd Xte = MatGrid3d( Xtr.limRows( 0.5 ) , hres );

/////// GRID

    float gres = 0.1;
    Matd Xgr = MatGrid3d( Xtr.limRows() , gres );
    Matd Ygr( Xgr.r() ); Ygr.setVal(0);

    Mati siz = Xgr.gridSize();

    Matd Xmin = Xtr - Xtr.minRows();
    Mati Xidx = ( Xmin / gres ).Round().toInt();

    forLOOPi( Xtr.r() )
    {
        int ii = Xidx(i,2) + Xidx(i,1) * siz(2) + Xidx(i,0) * siz(1) * siz(2);
        Ygr(ii) = 1;
    }

////// INFO

    Xtr.info();
    Xte.info();
    Xgr.info();

////// HILBERT MAPS

    Timer t;

    HMreal hm( new HMfeatSqExp() ,
               new HMregrLogistic() );

    t.tick();
    SMatd P,M,S;
    hm.cluster( Xtr , P , M , S , rad , thr );
    t.ptick( "clustering" );

    t.tick();
    hm.add( P , M , S );
    hm.weights().setVal(1);
    t.ptick( "training" );

    t.tick();
    Matd Yte = hm.query( Xte );
    t.ptick( "querying" );

    t.tick();
    Matd srf = marching3D( Xte , Yte > 0 , Yte , 0.5 );
    t.ptick( "marching" );

    disp( srf.sum() );

/////// DRAW

    CPPlot draw("Window");
    draw[0].set3Dworld().setViewer( -4.82881 , 3.36886 , 0.3309260 ,
                                    -3.86404 , 3.27373 , 0.0856423 );

    unsigned buf_pts = draw.addBuffer3D( Xtr );
    unsigned buf_srf = draw.addBuffer3D( srf );
    unsigned buf_ctr = draw.addBuffer3D( hm.ctrs() );

    unsigned buf_grd0 = draw.addBuffer3D( Xgr.sr( Ygr == 0 ) );
    unsigned buf_grd1 = draw.addBuffer3D( Xgr.sr( Ygr == 1 ) );

    int show = 0;
    while( draw.input() )
    {
        draw[0].clear();

        switch( show )
        {
        case 0:
            draw.psc(2,WHI).pts3D( buf_pts );
            draw.psc(4,RED).pts3D( buf_ctr );
            break;
        case 1:
            draw.clr(BLU).wsurf3D( buf_srf );
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






