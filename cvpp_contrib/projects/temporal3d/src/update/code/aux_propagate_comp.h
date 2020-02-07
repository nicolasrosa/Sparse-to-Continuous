#ifndef AUX_PROPAGATE_H
#define AUX_PROPAGATE_H

#include <cvpp/containers/matrix.h>
#include <cvpp/algorithms/icp/libicp.h>

#include "class_map.h"
#include "class_group.h"

using namespace cvpp;

/// Propagate objects from Next to Main

void propagateObjects( Map& prev , Map& next , Map& main ,
                       const int& thrctr , const double& thrfit )
{
    SeqMatd M1p( main.size() ) , M1n( main.size() );
    auto Mobs = MatZEROSd( main.size() );

    /// Calculate movement

    #pragma omp parallel for
    forLOOPi( prev.size() )
    {
        if( prev[i].size() >= thrctr && next[i].size() >= thrctr &&
            0.8 * prev[i].size() < next[i].size() )
        {
            Matd Mp = prev[i].getM(); // Changed to P
            Matd Mn = next[i].getM(); // Changed to P

            LibICPd icp( Mn );
            icp.fit( Mp , main[i].R , main[i].t , thrfit );

            Mobs(i) = true;
        }
    }

    /// Move main groups

    forLOOPi( main.size() )
    {
        M1p[i] = main[i].getM();
        M1n[i] = icp_calc( M1p[i] , main[i].R , main[i].t );
        main[i].updM( M1n[i] );
    }

    // Remove main clusters based on next

    forLOOPi( main.size() )
    {
        if( main[i].size() >= thrctr && next[i].size() >= thrctr )
        {
            Matd Mm = main[i].getM();
            Matd Mn = next[i].getM();

            SSeqi idxs; SSeqd dsts; KDtreed tree( Mn );
            tree.knnSearch( Mm , 1 , idxs , dsts );

            Veci rems;
            forLOOPj( idxs.size() )
                if( dsts[j][0] < 0.5 * main.rad * main.rad ) rems.insert(j);
            rems.update(); rems.mat().SortRows().FlipRows();
            forLOOPj( rems.n() ) main[i].del( rems[j] );
        }
    }

    /// Copy next clusters to main

    forLOOPij( main.size() , next[i].size() )
        main[i].push_back( next[i][j] );

    forLOOPii( main.size() , next.size() )
        main.objs.push_back( next[i] );

    // Recalculate main movement

    #pragma omp parallel for
    forLOOPi( M1p.size() )
    {
        if( Mobs(i) )
        {
            Matd M2 = main[i].getM();
            if( M1p[i].r() >= thrctr && M2.r() >= thrctr )
            {
                LibICPd icp( M2 );
                icp.fit( M1p[i] , main[i].R , main[i].t , thrfit );
                Matd d = ( M2 - M1p[i] ).meanRows();

                double chg = ( main[i].v - d ).rsqsum();
                if( chg > 0.9 )
                {
                    main[i].t.setVal(0);
                    main[i].R.setIdentity();
                }
                else
                {
                    if( main[i].v(0) == 0.0 ) main[i].v = d;
                    else main[i].v = ( main[i].v + d ) / 2.0;
                    if( d.rsqsum() < 0.1 ) main[i].v.setVal(0);

                    main[i].t = main[i].v; main[i].R.setIdentity();
                }
            }
        }
        else
        {
            Matd d = ( M1n[i] - M1p[i] ).meanRows();

            if( main[i].v(0) == 0.0 ) main[i].v = d;
            else main[i].v = ( main[i].v + d ) / 2.0;
            main[i].t = main[i].v; main[i].R.setIdentity();
        }

//        main[i].t(2) = 0.0; // No Z motion
    }

    forLOOPi( main.size() ) main[i].t.print();
}

#endif

