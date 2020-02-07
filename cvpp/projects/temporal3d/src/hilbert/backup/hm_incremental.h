#ifndef HM_INCR_H
#define HM_INCR_H

#include <cvpp/algorithms/kdtree/kdtree.h>

namespace cvpp
{

class HMincr
{

public:

    int n;
    Matd ctrs,weights,qty;
    SeqMatd means,vars,covs,icovs;

public:

    HMincr()
    {
        this->n = 0;
    }

    ~HMincr()
    {

    }

    const void add( const SeqMatd& means , const SeqMatd& covs , const Matd& weights )
    {
        int nn = n + means.size();
        this->means.resize( nn ); this->vars.resize(  nn );
        this->covs.resize(  nn ); this->icovs.resize( nn );

        #pragma omp parallel for
        forLOOPii( n , nn )
        {
            int in = i - n;
            this->means[i] = means[in].clone();
            this->vars[i] = covs[in].clone();

            this->covs[i] = 5.0 * covs[in];
            this->icovs[i] = this->covs[i].inv();
        }

        this->ctrs.AppD( Matd( means ) );
        this->weights.AppD( weights );

        n = nn;
    }

    Matd query3( const Matd& pts , const int& nn = 5 )
    {
        KDtreed kdtree( ctrs );
        return query3( pts , kdtree , nn );
    }

    Matd query3( const Matd& pts , KDtreed& kdtree ,
                 const int& nn = 5 )
    {
        Matd res = MatZEROSd( pts.r() );

        SSeqi idxs; SSeqd dsts;
        kdtree.knnSearch( pts , nn , idxs , dsts );

        #pragma omp parallel for
        forLOOPi( idxs.size() )
        {
            forLOOPj( idxs[i].size() )
            {
                int k = idxs[i][j];

                Eig13d d = pts.row(i) - ctrs.row(k);
                Eig33d iS = icovs[k].eig();

                res(i) += weights(k) * std::exp( - ( d * iS * d.transpose() )(0) );
            }
        }

//        res = 1.0 / ( 1.0 + ( - res ).exp() );
        return res;
    }

    Matd query3( const Matd& pts , const Seqi& idxs )
    {
        Matd res = MatZEROSd( pts.r() );

        forLOOPij( pts.r() , idxs.size() )
        {
            int k = idxs[j];

            Eig13d d = pts.row(i) - ctrs.row(k);
            Eig33d iS = icovs[k].eig();

            res(i) += weights(k) * std::exp( - ( d * iS * d.transpose() )(0) );
        }

        return res;
    }

};

}

#endif
