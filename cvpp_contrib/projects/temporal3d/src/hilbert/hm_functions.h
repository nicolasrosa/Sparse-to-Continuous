#ifndef HM_FUNCTIONS_H
#define HM_FUNCTIONS_H

#include <cvpp/containers/matrix.h>
#include <cvpp/algorithms/clustering/quick_means.h>

using namespace cvpp;

// HM DATA
struct HMdata
{
    SeqMatd P,M,S;
    Matd W;
};

struct Model
{
    Matd R,t;
    Matd C,v;

    Model()
    {
        R.reset( 3 , 3 ).setIdentity();
        t.reset( 1 , 3 ).setVal( 0.0 );

        C.reset( 3 , 3 ).setVal( 0.0 );
        v.reset( 1 , 3 ).setVal( 0.0 );
    }

    void print()
    {
        disp( t(0) , t(1) , std::acos( R(0,0) ) );
    }
};

using SeqModel = Seq<Model>;

/// CLUSTER DATA
void cluster_data( const SeqMatd& xyz ,
                   const int& i , const double& rad ,
                   HMdata& HM , const double& w )
{
    quick_means( xyz[i] , HM.P , HM.M , HM.S , rad , 0.5 * rad , 3 );
    HM.W = w * MatONESd( HM.P.size() );
}

#endif
