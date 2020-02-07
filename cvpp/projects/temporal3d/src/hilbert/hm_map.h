#ifndef HM_MAP_H
#define HM_MAP_H

#include <cvpp/interfaces/cpplot.h>
#include <cvpp/algorithms/icp/libicp.h>

#include <cvpp/algorithms/hilbert_maps/real/hm_real.h>
#include <cvpp/algorithms/hilbert_maps/real/feats/hm_real_feat_sqexp.h>
#include <cvpp/algorithms/hilbert_maps/real/regrs/hm_real_regr_logistic.h>
#include <cvpp/algorithms/marching_cubes/marching3D.h>

#include "hm_functions.h"

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
    SSeqMatd mov;

public:

    Map()
    {
    }

    Map( const SeqMatd& dataoccs , SeqMatd& datafree ,
         const int& t , const double& rad )
    {
        init( dataoccs , datafree , t , rad );
    }

    void init( const SeqMatd& dataoccs , const SeqMatd& datafree ,
               const int& t , const double& rad )
    {
        this->rad = rad;
        cluster_data( dataoccs , t , rad , occs , +1.0 );
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

    Matd getObject( const int& i )
    {
        Matd M = occs.M;
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

        Matd M = occs.M;
        KDtreed tree( M ); SSeqd dsts; SSeqi idxs;
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
        }

        forLOOPi( O.max() )
        {
            Mati I = O == i + 1;

            if( I.r() < 5 )
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
                i--;
            }
        }

        forLOOPi( O.max() )
        {
            Mati I = O == i + 1;
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
//            draw.psc(2,BLA).ellipse2D( occs.M , occs.S ).psc(5,MAG).pts2D( Matd( occs.M ) );
//            draw.psc(2,CYA).ellipse2D( free.M , free.S ).psc(5,CYA).pts2D( Matd( free.M ) );
//        }
//    }

    void drawObjects( CPPlot& draw )
    {
        if( free.M.size() > 0 )
        forLOOPi( O.max() )
        {
            Mati Ii = O == i + 1;
            if( Ii.r() > 5 )
            {
                Matd Mi = Matd( occs.M ).sr( Ii );
                Matd Li = Mi.limRows( 1.5 * rad );

//                Scalar clr = RGBint(i);
                Scalar clr = BLA; int thk = 1;

                draw.lwc(thk,clr).line3D( Pt3f( Li(0,0) , Li(0,1) , Li(0,2) ) , Pt3f( Li(1,0), Li(0,1) , Li(0,2) ) );
                draw.lwc(thk,clr).line3D( Pt3f( Li(0,0) , Li(1,1) , Li(0,2) ) , Pt3f( Li(1,0), Li(1,1) , Li(0,2) ) );
                draw.lwc(thk,clr).line3D( Pt3f( Li(0,0) , Li(0,1) , Li(1,2) ) , Pt3f( Li(1,0), Li(0,1) , Li(1,2) ) );
                draw.lwc(thk,clr).line3D( Pt3f( Li(0,0) , Li(1,1) , Li(1,2) ) , Pt3f( Li(1,0), Li(1,1) , Li(1,2) ) );

                draw.lwc(thk,clr).line3D( Pt3f( Li(0,0) , Li(0,1) , Li(0,2) ) , Pt3f( Li(0,0), Li(1,1) , Li(0,2) ) );
                draw.lwc(thk,clr).line3D( Pt3f( Li(1,0) , Li(0,1) , Li(0,2) ) , Pt3f( Li(1,0), Li(1,1) , Li(0,2) ) );
                draw.lwc(thk,clr).line3D( Pt3f( Li(0,0) , Li(0,1) , Li(1,2) ) , Pt3f( Li(0,0), Li(1,1) , Li(1,2) ) );
                draw.lwc(thk,clr).line3D( Pt3f( Li(1,0) , Li(0,1) , Li(1,2) ) , Pt3f( Li(1,0), Li(1,1) , Li(1,2) ) );

                draw.lwc(thk,clr).line3D( Pt3f( Li(0,0) , Li(0,1) , Li(0,2) ) , Pt3f( Li(0,0), Li(0,1) , Li(1,2) ) );
                draw.lwc(thk,clr).line3D( Pt3f( Li(1,0) , Li(0,1) , Li(0,2) ) , Pt3f( Li(1,0), Li(0,1) , Li(1,2) ) );
                draw.lwc(thk,clr).line3D( Pt3f( Li(0,0) , Li(1,1) , Li(0,2) ) , Pt3f( Li(0,0), Li(1,1) , Li(1,2) ) );
                draw.lwc(thk,clr).line3D( Pt3f( Li(1,0) , Li(1,1) , Li(0,2) ) , Pt3f( Li(1,0), Li(1,1) , Li(1,2) ) );
            }
        }
    }

    void drawVelocities( CPPlot& draw )
    {
        if( free.M.size() > 0 )
        forLOOPi( O.max() )
        {
            forLOOPj( mov[i][3].r() )
                draw.lwc(3,BLA).line3D( mov[i][3].r(j) | mov[i][4].r(j) );
        }
    }

};

///////////////////////////////////////////////////////////////

/// ALIGN OBJECTS
void alignObjects( Map& main , Map& next )
{
    Matd M1 = main.occs.M;
    Matd M2 = next.occs.M;

    SSeqi idxs; SSeqd dsts;
    KDtreed tree( M1 );

    forLOOPi( next.O.max() )
    {
        Mati I2 = next.O == i + 1;
        if( I2.filled() )
        {
            Matd M2i = M2.sr( I2 );
            tree.knnSearch( M2i , 1 , idxs , dsts );

            if( dsts[0][0] < 2.5 )
            {
                int k = i + 1 , u = main.O( idxs[0][0] );
                forLOOPj( next.O.r() )
                    if( next.O(j) == k ) next.O(j) = - u;
            }
        }
    }

    int max = std::min( - main.O.max() , next.O.min() ) - 1;
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
}

/// PROPAGATE OBJECTS
void propagateObjects( Map& prev , Map& next , Map& main , const double& thrfit )
{
    int thrctr = 8;

    Mati Mobs( main.O.max() ); Mobs.setVal( 0 );
    bool observ = next.free.M.size() > 0;

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
        Matd M1,M2,M12;

        Mati I1 = prev.O == i + 1;
        Mati I2 = next.O == i + 1;

        if( I1.r() >= thrctr && I2.r() >= thrctr &&  0.8 * I1.r() < I2.r() )
        {
            Mobs(i) = 1;

            M1 = Matd( prev.occs.M ).sr( I1 );
            M2 = Matd( next.occs.M ).sr( I2 );

            LibICPd icp( M2 );
            M12 = icp.fit( M1 , main.T[i].R , main.T[i].t , thrfit );
        }

        main.mov[i][0] = M1;
        main.mov[i][1] = M2;
        main.mov[i][2] = M12;
    }

    SeqMatd mOp( main.O.max() );

    // Move Main Clusters
    disp( "Move Main Clusters" );

    forLOOPi( main.O.max() )
    {
        Mati I1 = main.O == i + 1;
        Matd M1 = Matd( main.occs.M ).sr( I1 );
        Matd M12 = icp_calc( M1 , main.T[i].R , main.T[i].t );

        forLOOPj( M12.r() )
            main.occs.M[ I1(j) ] = M12.r(j);

        main.mov[i][3] = M1;
        main.mov[i][4] = M12;

        mOp[i] = M1.clone();
    }

    // Remove Main Clusters based on Next
    disp( "Remove Main Clusters Based on Next" );

    forLOOPi( main.O.max() )
    {
        Mati I1 = main.O == i + 1;
        Mati I2 = next.O == i + 1;

        if( I1.r() >= thrctr && I2.r() >= thrctr )
        {
            Matd M1 = Matd( main.occs.M ).sr( I1 );
            Matd M2 = Matd( next.occs.M ).sr( I2 );

            SSeqi idxs; SSeqd dsts; KDtreed tree( M2 );
            tree.knnSearch( M1 , 1 , idxs , dsts );

            Veci rems;
            forLOOPj( idxs.size() )
            {
                    if( dsts[j][0] < 0.5 * main.rad * main.rad ) rems.insert( I1(j) );
//                    if( dsts[j][0] > 5.0 * main.rad * main.rad ) rems.insert( I1(j) );
            }
            rems.update(); rems.mat().SortRows().FlipRows();

            if( rems.n() < main.O.r() )
            {   forLOOPj( rems.n() )
                {
                    main.occs.P.erase( main.occs.P.begin() + rems[j] );
                    main.occs.M.erase( main.occs.M.begin() + rems[j] );
                    main.occs.S.erase( main.occs.S.begin() + rems[j] );
                    main.occs.W.DelRows( rems[j] );
                }
                main.O.RemoveRows( rems );
            }
        }
    }

    // Add Main Clusters Based on Next
    disp( "Add Main Clusters Based on Next" );

    forLOOPi( next.O.max() )
    {
        Mati I2 = next.O == i + 1;
        if( I2.r() >= 0 )
        {
            Matd M2 = Matd( next.occs.M ).sr( I2 );

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

    // Recalculate Movement
    disp( "Recalculate Movement" );

    forLOOPi( mOp.size() )
    {
        Matd M1 = mOp[i];

        Mati I2 = main.O == i + 1;
        Matd M2 = Matd( main.occs.M ).sr( I2 );

        if( observ && Mobs(i) )
        {
            if( M1.r() >= thrctr && M2.r() >= thrctr )
            {
                LibICPd icp( M2 );
                icp.fit( M1 , main.T[i].R , main.T[i].t , thrfit );

                Matd d = ( M2 - M1 ).meanRows();

                double chg = ( main.T[i].v - d ).rsqsum();
                if( chg > 0.9 )
                {
                    main.T[i].t.setVal(0); main.T[i].R.setIdentity();
                }
                else
                {
                    if( main.T[i].v(0) == 0.0 ) main.T[i].v = d;
                    else main.T[i].v = ( main.T[i].v + d ) / 2.0;
                    if( d.rsqsum() < 0.1 ) main.T[i].v.setVal(0);

                    main.T[i].t = main.T[i].v; main.T[i].R.setIdentity();
                }
            }
        }

        if( !observ || !Mobs(i) )
        {
            Matd M1 = main.mov[i][3];
            Matd M2 = main.mov[i][4];

            Matd d = ( M2 - M1 ).meanRows();

            if( main.T[i].v(0) == 0.0 ) main.T[i].v = d;
            else main.T[i].v = ( main.T[i].v + d ) / 2.0;
            main.T[i].t = main.T[i].v; main.T[i].R.setIdentity();
        }
    }

    // Update Weights

    forLOOPi( Mobs.r() )
    {
        Mati I = main.O == i + 1;

        forLOOPj( I.r() )
            main.occs.S[ I(j) ] -= main.T[i].C;

        if( !observ || !Mobs(i) )
        {
            forLOOPj( I.r() )
                main.occs.W( I(j) ) *= 0.9;

//            Matd add = main.T[i].v.abs() / 200.0 , N( 3 , 3 );
//            N.eig() << add(0) , 0.0000 , 0.0000 ,
//                       0.0000 , add(1) , 0.0000 ,
//                       0.0000 , 0.0000 , add(2) ;
//            main.T[i].C += N;
        }
        else
        {
            forLOOPj( I.r() )
                main.occs.W( I(j) ) = 1.0;

            main.T[i].C.setVal( 0.0 );
        }

        forLOOPj( I.r() )
            main.occs.S[ I(j) ] += main.T[i].C;
    }

    // Remove Low Weight Clusters

    {
        Veci rems;
        forLOOPi( main.occs.W.r() )
            if( main.occs.W(i) < 0.20 ) rems.insert( i );
        rems.update(); rems.mat().SortRows().FlipRows();

        forLOOPi( rems.n() )
        {
            main.occs.P.erase( main.occs.P.begin() + rems[i] );
            main.occs.M.erase( main.occs.M.begin() + rems[i] );
            main.occs.S.erase( main.occs.S.begin() + rems[i] );
        }
        main.occs.W.RemoveRows( rems.mat() );
        main.O.RemoveRows( rems.mat() );
    }

    // Remove Small Objects

    forLOOPi( main.O.max() )
    {
        Mati I = main.O == i + 1;
        if( I.r() < 5 )
        {
            I.SortRows().FlipRows();
            forLOOPj( I.r() )
            {
                main.occs.P.erase( main.occs.P.begin() + I(j) );
                main.occs.M.erase( main.occs.M.begin() + I(j) );
                main.occs.S.erase( main.occs.S.begin() + I(j) );
            }
            main.occs.W.RemoveRows( I );
            main.O.RemoveRows( I );
        }
    }
}

// CREATE HILBERT MAP

void create_hm( Map& main , Matd& grd , CPPlot& draw ,
                int& buf_srf , int& buf_clr )
{
    HMreal hm( new HMfeatSqExp() , new HMregrLogistic() );
    hm.add( main.occs.P , main.occs.M , main.occs.S ); hm.weights() = main.occs.W;
    Matd prb = hm.query( grd ); Matd srf = marching3D( grd , prb > 0 , prb , 0.7 );

    KDtreed kd( Matd( main.occs.M ) ); Matd clr( srf.r() );
    SSeqi idx; SSeqd dst; kd.knnSearch( srf , 1 , idx , dst );
    for( unsigned i = 0 ; i < idx.size() ; i++ )
        if( idx[i].size() > 0 ) clr(i) = main.occs.W( idx[i][0] );

    buf_srf = draw.addBuffer3D( srf );
//    buf_clr = draw.addBufferRGBjet( srf.c(2).clone() );
    buf_clr = draw.addBufferRGBjet( clr , 0.0 , 1.0 );
}

}


#endif
