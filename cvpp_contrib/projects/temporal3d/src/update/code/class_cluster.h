#ifndef CLASS_CLUSTER_H
#define CLASS_CLUSTER_H

#include <cvpp/containers/matrix.h>

using namespace cvpp;

/// Cluster struct (stores pointcloud (P) , mean (M) and covariance (S)

struct Cluster
{
    Matd P,M,S;
    double wgt;

    Cluster()
    {
    }

    Cluster( const Matd& P , const Matd& M , const Matd& S , const double& wgt )
    {
        this->P = P.clone();
        this->M = M.clone();
        this->S = S.clone();
        this->wgt = wgt;
    }

    Cluster( const Cluster& cluster )
    {
        this->P = cluster.P.clone();
        this->M = cluster.M.clone();
        this->S = cluster.S.clone();
        this->wgt = cluster.wgt;
    }
};

using SeqCluster = Seq<Cluster>;

#endif
