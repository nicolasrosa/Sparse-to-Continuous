#ifndef AUX_DETECT_H
#define AUX_DETECT_H

#include <cvpp/containers/matrix.h>

#include "class_map.h"
#include "class_group.h"

using namespace cvpp;

/// Detect objects from a group (recursive)

void
detectObjectsREC( const SSeqi& idxs , Mati& O , const int& i )
{
    forLOOPj( idxs[i].size() )
    {
        int k = idxs[i][j];
        if( O(k) != O(i) )
        {
            O(k) = O(i);
            detectObjectsREC( idxs , O , k );
        }
    }
}

/// Detect objects from a group (main function)

SeqGroup
detectObjects( Group& occs , const double& rad ,
               const double& mlt = 8 , const int& min = 5 )
{
    Matd M = occs.getM();
    Mati O = MatZEROSi( M.r() );

    KDtreed tree( M ); SSeqd dsts; SSeqi idxs;
    tree.radSearch( M , mlt * rad , idxs , dsts );

    int cnt = 0;
    forLOOPi( idxs.size() )
    {   if( O(i) == 0 )
        {
            O(i) = ++cnt;
            detectObjectsREC( idxs , O , i );
        }
    }

    forLOOPi( O.max() )
    {
        Mati I = O == i + 1;

        if( I.r() < min )
        {
            I.SortRows().FlipRows();
            forLOOPj( I.r() ) occs.del( I[j] );
            O.RemoveRows( I );
            forLOOPj( O.r() ) if( O(j) > i ) O(j)--;
            i--;
        }
    }

    SeqGroup objs( O.max() );
    forLOOPi( O.r() ) objs[ O(i) - 1 ].push_back( occs[i] );

    return objs;
}

#endif
