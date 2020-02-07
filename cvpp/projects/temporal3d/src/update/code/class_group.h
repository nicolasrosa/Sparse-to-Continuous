#ifndef CLASS_GROUP_H
#define CLASS_GROUP_H

#include <cvpp/containers/matrix.h>
#include <cvpp/containers/vector.h>
#include "class_cluster.h"

using namespace cvpp;

class Group
{
public:

    SeqCluster clusters;
    Matd R,C,t,v;
    Pts3d pos,vel,acc;

    Group()
    {
        R = MatIDTYd( 3 );
        t = MatZEROSd( 1 , 3 );

        v = MatZEROSd( 1 , 3 );
        C = MatZEROSd( 3 , 3 );

        clusters.clear();
    }

    /// Recover cluster given index

    Cluster& operator[]( const int& i ) { return clusters[i]; }
    const Cluster& operator[]( const int& i ) const { return clusters[i]; }

    /// Recover number of clusters

    int size() const { return clusters.size(); }

    /// Push back cluster into group

    void push_back( const Cluster& cluster )
    {
        clusters.push_back( Cluster( cluster ) );
    }

    /// Set group given initial values

    void set( const SeqMatd& P , const SeqMatd& M , const SeqMatd& S , const double& w )
    {
        clusters.resize( P.size() );

        #pragma omp parallel for
        forLOOPi( P.size() ) clusters[i] = Cluster( P[i] , M[i] , S[i] , w );
    }

    /// Delete cluster given index

    void del( const int& i )
    {
        clusters.erase( clusters.begin() + i );
    }

    /// Return matrix with M values

    Matd getM() const
    {
        Matd M( clusters.size() , 3 );
        forLOOPi( M.r() ) M.row(i) = clusters[i].M.eig();
        return M;
    }

    void updM( const Matd& M )
    {
        forLOOPi( M.r() )
            clusters[i].M.eig() = M.row(i);
    }

    /// Return matrix with P values

    Matd getP() const
    {
        Matd P( 0 , 3 );
        forLOOPi( clusters.size() ) P |= clusters[i].P;
        return P;
    }

    void savePos()
    {
        Matd p = getM().meanRows();
        pos.push( Pt3d( p(0),  p(1) , p(2) ) );
        if( pos.n() > 1 ) vel.push( pos[-1] - pos[-2] );
        if( vel.n() > 1 ) acc.push( vel[-1] - vel[-2] );
    }

    Matd calcMov()
    {
        Matd v;

        int m = std::min( 5 , pos.n() );
        if( m < 5 )
        {
            v = Mat13d( 0 , 0 , 0 );
//            a = Mat13d( 0 , 0 , 0 );
        }
        else
        {
            v = vel.mat().rd(m).meanRows();
//            a = acc.mat().rd(m-1).meanRows();
        }

        return v;
    }

};

using SeqGroup = Seq<Group>;

#endif
