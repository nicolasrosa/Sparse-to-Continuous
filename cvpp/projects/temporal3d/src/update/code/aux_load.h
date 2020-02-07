#ifndef AUX_LOAD_H
#define AUX_LOAD_H

#include <cvpp/containers/matrix.h>
#include <cvpp/containers/vector.h>

using namespace cvpp;

/// Load KITTI dataset

void
loadKITTIproc( const String& file , const int& n ,
               SeqMatd& occs , SeqMatd& free , SeqMatd& grnd, Matd& circ ,
               const double& reach )
{
    Matd poses( file + "/poses.bin" );
    occs.resize( n ); free.resize( n ); grnd.resize( n );

    #pragma omp parallel for
    forLOOPi( n )
    {
        occs[i].loadBIN( file + "/pts/" + toString( i , 10 ) + ".bin" );
        grnd[i].loadBIN( file + "/grd/" + toString( i , 10 ) + ".bin" );

        Matd org = poses.r(i).cl(3);

        occs[i].SampleRows( ( ( occs[i].cl(2) - org.cl(2) ).sqsumCols().sqrt() ) < reach );
        grnd[i].SampleRows( ( ( grnd[i].cl(2) - org.cl(2) ).sqsumCols().sqrt() ) < reach );

        Pts3d pts; pts.reserve( 10 * occs[i].r() );
        forLOOPj( occs[i].r() )
        {   for( double k = 0.70 ; k <= 0.95 ; k += 0.05 )
                pts.insert( Pt3d( org(0) + k * ( occs[i](j,0) - org(0) ) ,
                                  org(1) + k * ( occs[i](j,1) - org(1) ) ,
                                  org(2) + k * ( occs[i](j,2) - org(2) ) ) );
        }
        pts.update(); pts.mat().AddRand( -0.1 , 0.0 );
        free[i] = pts.mat().clone();
    }

    circ.reset( 360 * 10 , 3 );
    forLOOPi( circ.r() )
    {
        double ang = double(i) / ( 2 * PI * 10 );
        circ.row(i) << reach * std::cos(ang) , reach * std::sin(ang) , -2.0 ;
    }
}

#endif
