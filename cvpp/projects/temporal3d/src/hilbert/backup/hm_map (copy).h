#ifndef HM_MAP_H
#define HM_MAP_H

#include <cvpp/interfaces/cpplot.h>
#include <cvpp/algorithms/icp/libicp.h>

#include "hm_functions.h"
#include "hm_incremental.h"

namespace cvpp
{

// MAP
class Map
{

public:

    double rad;
    HMdata occs,free,stor;

    Mati O;
    SeqModel T;
    SSeqMatf mov;

public:

    Map()
    {
    }

    Map( const SeqMatf& dataoccs , SeqMatf& datafree ,
         const int& t , const double& rad )
    {
        init( dataoccs , datafree , t , rad );
    }

    void init( const SeqMatf& dataoccs , const SeqMatf& datafree ,
               const int& t , const double& rad )
    {
        this->rad = rad;
        cluster_data( dataoccs , t , rad , occs , +3.0 );
        cluster_data( datafree , t , rad , free , -1.0 );
    }

    Map copy()
    {
        Map novel;
        novel.rad = rad;

        forLOOPi( occs.M.size() )
        {
            novel.occs.P.push_back( occs.P[i].clone() );
            novel.occs.M.push_back( occs.M[i].clone() );
            novel.occs.S.push_back( occs.S[i].clone() );
        }
        novel.occs.W = occs.W.clone();

        forLOOPi( free.M.size() )
        {
            novel.free.P.push_back( free.P[i].clone() );
            novel.free.M.push_back( free.M[i].clone() );
            novel.free.S.push_back( free.S[i].clone() );
        }
        novel.free.W = free.W.clone();

        forLOOPi( stor.M.size() )
        {
            novel.stor.P.push_back( stor.P[i].clone() );
            novel.stor.M.push_back( stor.M[i].clone() );
            novel.stor.S.push_back( stor.S[i].clone() );
        }
        novel.stor.W = stor.W.clone();

        forLOOPi( T.size() )
        {
            novel.T.push_back( Model() );
            novel.T.back().R = T[i].R.clone();
            novel.T.back().t = T[i].t.clone();
        }

        novel.O = O.clone();

        return novel;
    }

////////////////// GET STUFF

    Matf getObject( const int& i )
    {
        Matf M = occs.M;
        return M.sr( O == i + 1 );
    }

////////////////// MAKE HILBERT

//    HMincr makeHM()
//    {
//        HMincr hm;
//        hm.add( occs.M , occs.S , occs.W );
//        hm.add( free.M , free.S , free.W );
//        return hm;
//    }

//    HMincr makeHMoccs()
//    {
//        HMincr hm;
//        hm.add( occs.M , occs.S , occs.W );
//        return hm;
//    }

//    HMincr makeHMfree()
//    {
//        HMincr hm;
//        hm.add( free.M , free.S , free.W );
//        return hm;
//    }

////////////////// DETECT OBJECTS

    void detectObjects( const double& rad )
    {
        O.reset( occs.M.size() );
        O.setVal( 0 );

        Matf M = occs.M;
        KDtreef tree( M ); SSeqf dsts; SSeqi idxs;
        tree.radSearch( M , 8 * rad , idxs , dsts );

        int cnt = 0;
        forLOOPi( idxs.size() )
        {   if( O(i) == 0 )
            {
                O(i) = ++cnt;
                detectObjectsREC( idxs , i , O );
            }
        }

        forLOOPi( O.max() )
        {
            Mati I = O == i + 1;

            if( I.r() < 3 )
            {
                I.SortRows().FlipRows();

                forLOOPj( I.r() )
                {
                    occs.P.erase( occs.P.begin() + I[j] );
                    occs.M.erase( occs.M.begin() + I[j] );
                    occs.S.erase( occs.S.begin() + I[j] );
                    occs.W.DelRows( I[j] );
                }

                O.RemoveRows( I );

                forLOOPj( O.r() )
                    if( O(j) > i ) O(j)--;
            }
        }

        T.resize( O.max() );
    }

    void detectObjectsREC( const SSeqi& idxs ,
                           const int& i , Mati& O )
    {
        forLOOPj( idxs[i].size() )
        {
            int k = idxs[i][j];
            if( O(k) != O(i) )
            {
                O(k) = O(i);
                detectObjectsREC( idxs , k , O );
            }
        }
    }

////////////////// DRAW STUFF

//    void drawEllipses( CPPlot& draw )
//    {
//        if( free.M.size() > 0 )
//        {
//            draw.psc(2,BLA).ellipse2D( occs.M , occs.S ).psc(5,MAG).pts2D( Matf( occs.M ) );
//            draw.psc(2,CYA).ellipse2D( free.M , free.S ).psc(5,CYA).pts2D( Matf( free.M ) );
//        }
//    }

//    void drawObjects( CPPlot& draw )
//    {
//        if( free.M.size() > 0 )
//        forLOOPi( O.max() )
//        {
//            Matf Mi = Matf( occs.M ).sr( O == i + 1 );
//            Matf Li = Mi.limRows( rad );

//            draw.lwc(3,RGBint(i)).line2D( Pt2f( Li(0,0) , Li(0,1) ) , Pt2f( Li(0,0), Li(1,1) ) ); // 0,0
//            draw.lwc(3,RGBint(i)).line2D( Pt2f( Li(1,0) , Li(0,1) ) , Pt2f( Li(1,0), Li(1,1) ) ); // 1,0
//            draw.lwc(3,RGBint(i)).line2D( Pt2f( Li(0,0) , Li(0,1) ) , Pt2f( Li(1,0), Li(0,1) ) ); // 0,1
//            draw.lwc(3,RGBint(i)).line2D( Pt2f( Li(0,0) , Li(1,1) ) , Pt2f( Li(1,0), Li(1,1) ) ); // 1,1
//        }
//    }

    void drawVelocities( CPPlot& draw )
    {
        if( free.M.size() > 0 )
        forLOOPi( O.max() )
        {
            forLOOPj( mov[i][3].r() )
                draw.lwc(3,YEL).line3D( mov[i][3].r(j) | mov[i][4].r(j) );
        }
    }

};

///////////////////////////////////////////////////////////////

/// ALIGN OBJECTS
void alignObjects( Map& main , Map& next )
{
    Matf M1 = main.occs.M;
    Matf M2 = next.occs.M;

    SSeqi idxs; SSeqf dsts;
    KDtreef tree( M1 );

    forLOOPi( next.O.max() )
    {
        Mati I2 = next.O == i + 1;
        if( I2.filled() )
        {
            Matf M2i = M2.sr( I2 );
            tree.knnSearch( M2i , 1 , idxs , dsts );

            if( dsts[0][0] < 2.5 )
            {
                int k = i + 1 , u = main.O( idxs[0][0] );
                forLOOPj( next.O.r() )
                    if( next.O(j) == k ) next.O(j) = - u;
            }
        }
    }

    int max = next.O.min() - 1;
    forLOOPi( next.O.r() )
    {   if( next.O(i) > 0 )
        {
            int k = next.O(i);
            forLOOPj( next.O.r() )
                if( next.O(j) == k ) next.O(j) = max;
            max--;
        }
    }

    main.O.Abs();
    next.O.Abs();

//    Matf L1( main.O.max() , 3 );
//    forLOOPi( L1.r() )
//            L1.r(i) = Matf( main.occs.M ).sr( main.O == i + 1 ).meanRows();

//    Matf L2( next.O.max() , 3 );
//    forLOOPi( L2.r() )
//            L2.r(i) = Matf( next.occs.M ).sr( next.O == i + 1 ).meanRows();

//    SeqModel TT( main.T.size() );

//    int cnt = -1;
//    forLOOPi( L1.r() )
//    {
//        double bst = PINF; int k = 0;
//        forLOOPj( L2.r() )
//        {
//            double dst = ( L1.r(i) - L2.r(j) ).sqsum();
//            if( dst < bst ) { bst = dst; k = j; }
//        }

//        if( bst < 2.5 )
//        {
//            TT.push_back( main.T[i] );

//            forLOOPu( main.O.r() ) if( main.O(u) == i + 1 ) main.O(u) = cnt;
//            forLOOPu( next.O.r() ) if( next.O(u) == k + 1 ) next.O(u) = cnt;
//            cnt--;
//        }

//        disp( i , k , cnt , bst );
//    }

//    main.T = TT;

//    main.O.uniqueRows().print();
//    next.O.uniqueRows().print();

//    forLOOPi( main.O.r() )
//    {
//        if( main.O(i) > 0 )
//        {
//            forLOOPjj( i + 1 , main.O.r() )
//                if( main.O(j) == main.O(i) ) main.O(j) = cnt;
//            main.O(i) = cnt; cnt--;
//        }
//    }

//    forLOOPi( next.O.r() )
//    {
//        if( next.O(i) > 0 )
//        {
//            forLOOPjj( i + 1 , next.O.r() )
//                if( next.O(j) == next.O(i) ) next.O(j) = cnt;
//            next.O(i) = cnt; cnt--;
//        }
//    }

//    main.O.uniqueRows().print();
//    next.O.uniqueRows().print();

//    main.O.Abs();
//    next.O.Abs();
}

/// PROPAGATE OBJECTS
void propagateObjects( Map& prev , Map& next , Map& main , const double& thrfit )
{
//    Mati Mobs( main.O.max() ); Mobs.setVal( 0 );
//    bool observ = next.free.M.size() > 0;

    // Prepare Matrices
    disp( "Prepare Matrices" );

    int n_main = std::max( next.O.max() , main.O.max() );

    main.mov.resize( n_main );
    forLOOPi( main.mov.size() ) main.mov[i].resize( 5 );
    while( main.T.size() < n_main ) main.T.push_back( Model() );

    // Calculate Prev-Next Movement
    disp( "Calculate Prev-Next Movement" );

    forLOOPi( prev.O.max() )
    {
        Matf M1,M2,M12;

        Mati I1 = prev.O == i + 1;
        Mati I2 = next.O == i + 1;

        if( I1.r() >= 7 && I2.r() >= 7 )
        {
            M1 = Matf( prev.occs.M ).sr( I1 );
            M2 = Matf( next.occs.M ).sr( I2 );

            LibICPf icp( M2 );
            M12 = icp.fit( M1 , main.T[i].R , main.T[i].t , thrfit );
        }

        main.mov[i][0] = M1;
        main.mov[i][1] = M2;
        main.mov[i][2] = M12;
    }

    SeqMatf mOp( main.O.max() );

    // Move Main Clusters
    disp( "Move Main Clusters" );

    forLOOPi( main.O.max() )
    {
        Mati I1 = main.O == i + 1;
        Matf M1 = Matf( main.occs.M ).sr( I1 );
        Matf M12 = icp_calc( M1 , main.T[i].R , main.T[i].t );

        forLOOPj( M12.r() )
            main.occs.M[ I1(j) ] = M12.r(j);

        main.mov[i][3] = M1;
        main.mov[i][4] = M12;

        mOp[i] = M1.clone();
    }

    // Add Main Clusters Based on Next
    disp( "Add Main Clusters Based on Next" );

    forLOOPi( next.O.max() )
    {
        Mati I2 = next.O == i + 1;
        if( I2.r() >= 0 )
        {
            Matf M2 = Matf( next.occs.M ).sr( I2 );

            forLOOPj( M2.r() )
            {
                main.occs.M.push_back( next.occs.M[ I2(j) ].clone() );
                main.occs.P.push_back( next.occs.P[ I2(j) ].clone() );
                main.occs.S.push_back( next.occs.S[ I2(j) ] + main.T[i].C );
                main.occs.W |= next.occs.W.r( I2(j) );

                Mati O( 1 ); O(0) = i + 1;
                main.O |= O;
            }
        }
    }

//    // Recalculate Movement
//    disp( "Recalculate Movement" );

//    forLOOPi( mOp.size() )
//    {
//        Matf M1 = mOp[i];

//        Mati I2 = main.O == i + 1;
//        Matf M2 = Matf( main.occs.M ).sr( I2 );

////        if( observ && Mobs(i) )
//        {
//            if( M1.r() >= 8 && M2.r() >= 8 )
//            {
//                disp( "recalculating" , i );

//                LibICPf icp( M2 );
//                icp.fit( M1 , main.T[i].R , main.T[i].t , thrfit );
//            }
//        }

////        if( !observ || !Mobs(i) )
////        {
////            disp( "filtering" , i );

////            double ang = std::acos( main.T[i].R(0,0) ) * 180 / PI;
////            if( std::fabs( ang ) < 10.0 ) main.T[i].R.setIdentity();
////        }

//        if( std::fabs( main.T[i].t(0) ) < main.rad / 2.0 ) main.T[i].t(0) = 0.0;
//        if( std::fabs( main.T[i].t(1) ) < main.rad / 2.0 ) main.T[i].t(1) = 0.0;
//    }



//    SeqMatf mOp( main.O.max() );


//    // Remove Main Clusters based on Next
//    forLOOPi( main.O.max() )
//    {
//        Mati I1 = main.O == i + 1;
//        Mati I2 = next.O == i + 1;

//        if( I1.r() >= 5 && I2.r() >= 5 )
//        {
//            Matf M1 = Matf( main.occs.M ).sr( I1 );
//            Matf M2 = Matf( next.occs.M ).sr( I2 );

//            SSeqi idxs; SSeqf dsts; KDtreef tree( M2 );
//            tree.knnSearch( M1 , 1 , idxs , dsts );

//            Veci rems;
//            forLOOPj( idxs.size() )
//                if( dsts[j][0] < 0.5 * main.rad * main.rad )
//                    rems.push( I1(j) );
//            rems.mat().SortRows().FlipRows();

//            if( rems.n() < main.O.r() )
//            {   forLOOPj( rems.n() )
//                {
//                    main.occs.P.erase( main.occs.P.begin() + rems[j] );
//                    main.occs.M.erase( main.occs.M.begin() + rems[j] );
//                    main.occs.S.erase( main.occs.S.begin() + rems[j] );
//                    main.occs.W.DelRows( rems[j] );
//                }
//                main.O.RemoveRows( rems );
//            }
//        }
//    }

//    disp( "eee" );

//    // Remove Main Clusters based on Free
//    if( observ )
//    {
//        Matf M1 = main.occs.M;
//        Matf M2 = next.free.M;

//        SSeqi idxs; SSeqf dsts; KDtreef tree( M2 );
//        tree.knnSearch( M1 , 1 , idxs , dsts );

//        Veci rems;
//        forLOOPj( idxs.size() )
//            if( dsts[j][0] < 0.5 * main.rad * main.rad )
//                rems.push( j );
//        rems.mat().SortRows().FlipRows();

//        forLOOPj( rems.n() )
//        {
//            main.occs.P.erase( main.occs.P.begin() + rems[j] );
//            main.occs.M.erase( main.occs.M.begin() + rems[j] );
//            main.occs.S.erase( main.occs.S.begin() + rems[j] );
//            main.occs.W.DelRows( rems[j] );
//        }
//        main.O.RemoveRows( rems );
//    }

//    disp( "fff" );



//    disp( "ggg" );

//    // Recalculate Movement
//    forLOOPi( mOp.size() )
//    {
//        Matf M1 = mOp[i];

//        Mati I2 = main.O == i + 1;
//        Matf M2 = Matf( main.occs.M ).sr( I2 );

//        if( observ && Mobs(i) )
//        {
//            if( M1.r() >= 8 && M2.r() >= 8 )
//            {
//                disp( "recalculating" , i );

//                LibICPf icp( M2 );
//                icp.fit( M1 , main.T[i].R , main.T[i].t , 0.2 );
//            }
//        }

//        if( !observ || !Mobs(i) )
//        {
//            disp( "filtering" , i );

//            double ang = std::acos( main.T[i].R(0,0) ) * 180 / PI;
//            if( std::fabs( ang ) < 10.0 ) main.T[i].R.setIdentity();
//        }

//        if( std::fabs( main.T[i].t(0) ) < main.rad / 2.0 ) main.T[i].t(0) = 0.0;
//        if( std::fabs( main.T[i].t(1) ) < main.rad / 2.0 ) main.T[i].t(1) = 0.0;
//    }

//    // Update Weights
//    forLOOPi( Mobs.r() )
//    {
//        Mati I = main.O == i + 1;

//        if( !observ || !Mobs(i) )
//        {
//            forLOOPj( I.r() )
//                main.occs.W( I(j) ) *= 0.9;
//        }
//        else
//        {
//            forLOOPj( I.r() )
//                main.occs.W( I(j) ) = 3.0;
//        }
//    }

//    // Calculate Velocities
//    forLOOPi( main.mov.size() )
//    {
//        main.T[i].v = ( main.mov[i][4] - main.mov[i][3] ).meanRows();

//        Mati I = main.O == i + 1;

//        forLOOPj( I.r() )
//            main.occs.S[ I(j) ] -= main.T[i].C;

//        if( !observ || !Mobs(i) )
//        {
//            Matf N( 3 , 3 );
//            N.eig() << main.T[i].v(0) / 20.0 , 0.000 , 0.000 ,
//                       0.000 , main.T[i].v(1) / 20.0 , 0.000 ,
//                       0.000 , 0.000 , main.T[i].v(2) / 20.0 ;
//            main.T[i].C += N;
//        }
//        else
//        {
//            main.T[i].C.setVal( 0.0 );
//        }

//        forLOOPj( I.r() )
//            main.occs.S[ I(j) ] += main.T[i].C;
//    }
}

///// UPDATE FREE
//void updateFree( Map& main , Map& next )
//{
//    double rad = main.rad;

//    {
//        Matf M1o = main.occs.M;
//        Matf M1f = main.free.M;

//        KDtreef treeo( M1o );
//        KDtreef treef( M1f );

//        if( next.free.M.size() > 0 )
//        {
//            Matf M2 = next.free.M;

//            SSeqi idxso,idxsf; SSeqf dstso,dstsf;
//            treeo.knnSearch( M2 , 1 , idxso , dstso );
//            treef.knnSearch( M2 , 1 , idxsf , dstsf );

//            forLOOPi( idxsf.size() )
//            {
//                if( dstsf[i][0] > 0.5 * rad * rad &&
//                    dstso[i][0] > 1.5 * rad * rad )
//                {
//                    main.free.P.push_back( next.free.P[i] );
//                    main.free.M.push_back( next.free.M[i] );
//                    main.free.S.push_back( next.free.S[i] );
//                    main.free.W |= next.free.W.r(i);
//                }
//            }
//        }
//    }

//    {
//        Matf M1o = main.occs.M;
//        Matf M1f = main.free.M;

//        KDtreef treeo( M1o );

//        SSeqi idxs; SSeqf dsts;
//        treeo.knnSearch( M1f , 1 , idxs , dsts );

//        Veci rems;
//        forLOOPi( idxs.size() )
//            if( dsts[i][0] < 1.0 * rad * rad ) rems.push( i );
//        rems.mat().SortRows().FlipRows();

//        forLOOPi( rems.n() )
//        {
//            main.stor.P.push_back( main.free.P[ rems[i] ].clone() );
//            main.stor.M.push_back( main.free.M[ rems[i] ].clone() );
//            main.stor.S.push_back( main.free.S[ rems[i] ].clone() );
//            main.stor.W |= main.free.W.r( rems[i] );

//            main.free.P.erase( main.free.P.begin() + rems[i] );
//            main.free.M.erase( main.free.M.begin() + rems[i] );
//            main.free.S.erase( main.free.S.begin() + rems[i] );
//            main.free.W.DelRows( rems[i] );
//        }
//    }
//}

///// RETURN FREE
//void returnFree( Map& main )
//{
//    double rad = main.rad;

//    Matf Mo = Matf( main.occs.M );
//    Matf Mf = Matf( main.free.M );
//    Matf Ms = main.stor.M;

//    KDtreef treeo( Mo );
//    KDtreef treef( Mf );

//    SSeqi idxso,idxsf; SSeqf dstso,dstsf;
//    treeo.knnSearch( Ms , 1 , idxso , dstso );
//    treef.knnSearch( Ms , 2 , idxsf , dstsf );

//    Veci rems;
//    forLOOPi( Ms.r() )
//    {   if( dstsf[i][1] > 1.0 * rad * rad &&
//            dstso[i][0] > 8.0 * rad * rad )
//        {
//            rems.push( i );
//        }
//    }
//    rems.mat().SortRows().FlipRows();

//    forLOOPi( rems.n() )
//    {
//        main.free.P.push_back( main.stor.P[ rems[i] ].clone() );
//        main.free.M.push_back( main.stor.M[ rems[i] ].clone() );
//        main.free.S.push_back( main.stor.S[ rems[i] ].clone() );
//        main.free.W |= main.stor.W.r( rems[i] );

//        main.stor.P.erase( main.stor.P.begin() + rems[i] );
//        main.stor.M.erase( main.stor.M.begin() + rems[i] );
//        main.stor.S.erase( main.stor.S.begin() + rems[i] );
//        main.stor.W.DelRows( rems[i] );
//    }

//}

}


#endif
