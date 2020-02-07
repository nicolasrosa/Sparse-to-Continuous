#ifndef AUX_ALIGN_H
#define AUX_ALIGN_H

#include <cvpp/containers/matrix.h>

#include "class_map.h"
#include "class_group.h"

using namespace cvpp;

/// Align objects given two groups (modify second one)

void alignObjects( Map& main , Map& next ,
                   const double& thr = 2.0 )
{
    SeqMatd M1( main.size() );
    forLOOPi( M1.size() ) M1[i] = main[i].getM();

    SeqMatd M2( next.size() );
    forLOOPi( M2.size() ) M2[i] = next[i].getM();

    Matd D = PINF * MatONESd( M2.size() );
    Mati I =  -1  * MatONESi( M2.size() );

    #pragma omp parallel for
    forLOOPi( M2.size() )
    {
        KDtreed tree( M2[i] );
        forLOOPj( M1.size() )
        {
            SSeqi idxs; SSeqd dsts;
            tree.knnSearch( M1[j] , 1 , idxs , dsts );

            double min = PINF;
            forLOOPk( dsts.size() ) if( dsts[k][0] < min ) min = dsts[k][0];
            if( min < D(i) ) { D(i) = min ; I(i) = j; }
        }
    }

    SeqGroup objs( M1.size() );
    forLOOPi( M2.size() )
    {
        if( D(i) < thr ) objs[ I(i) ] = next[i];
        else objs.push_back( next[i] );
    }

    next.objs = objs;
}

#endif
