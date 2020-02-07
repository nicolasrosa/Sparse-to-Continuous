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
    Timer t;

    SeqMatd M1p( main.size() ) , M1n( main.size() );
    auto Mobs = MatZEROSd( main.size() );

    /// Calculate movement
    t.ptick( "Calculating Movement" );

    #pragma omp parallel for
    forLOOPi( prev.size() )
    {
        if( prev[i].size() >= thrctr && next[i].size() >= thrctr &&
            0.8 * prev[i].size() < next[i].size() )
        {
            Matd Mp,Mn;
//            if( next[i].size() < 4 * thrctr )
//            {
//                Mp = prev[i].getP();
//                Mn = next[i].getP();
//            }
//            else
            {
                Mp = prev[i].getM();
                Mn = next[i].getM();
            }

            LibICPd icp( Mn );
            icp.fit( Mp , main[i].R , main[i].t , thrfit );

            Mobs(i) = true;
        }
    }

    /// Move main groups
    t.ptick( "Moving Main Groups" );

    forLOOPi( main.size() )
    {
        M1p[i] = main[i].getM();
        M1n[i] = icp_calc( M1p[i] , main[i].R , main[i].t );
        main[i].updM( M1n[i] );
    }

    // Remove main clusters based on next
    t.ptick( "Removing Main Clusters" );

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
    t.ptick( "Copying Next Clusters" );

    forLOOPij( main.size() , next[i].size() )
        main[i].push_back( next[i][j] );

    forLOOPii( main.size() , next.size() )
        main.objs.push_back( next[i] );

    // Recalculate main movement
    t.ptick( "Recalculating Main Movement" );

    #pragma omp parallel for
    forLOOPi( M1p.size() )
    {
//        main[i].v = main[i].calcMov();

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
        main[i].savePos();
    }

    forLOOPi( main.size() ) main[i].t.print();

//    // Update weights
//    t.ptick( "Updating Weights" );

//    forLOOPi( Mobs.r() )
//    {
//        forLOOPj( main[i].size() )
//            main[i][j].S -= main[i].C;

//        if( !Mobs(i) )
//        {
//            forLOOPj( main[i].size() )
//                main[i][j].wgt *= 0.95;
//        }
//        else
//        {
//            forLOOPj( main[i].size() )
//                main[i][j].wgt = 1.0;
//            main[i].C.setVal( 0.0 );
//        }

//        forLOOPj( main[i].size() )
//            main[i][j].S += main[i].C;
//    }

//    // Removing Uncertain Groups
//    t.ptick( "Removing Uncertain Groups" );

//    forLOOPi( main.size() )
//    {
//        double w = 0;
//        forLOOPj( main[i].size() ) w += main[i][j].wgt;
//        w /= main[i].size(); if( w > 0 && w < 0.25 ) main[i] = Group();
//    }

    t.ptick( "Done" );
}

#endif
