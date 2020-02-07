#ifndef CLASS_MAP_H
#define CLASS_MAP_H

#include <cvpp/interfaces/cpplot.h>
#include <cvpp/algorithms/clustering/quick_means.h>

#include "class_group.h"
#include "aux_detect.h"

using namespace cvpp;

class Map
{
public:

    double rad;
    SeqGroup objs;
    Group occs,free;

    /// Constructor

    Map()
    {
    }

    Map( const Matd& dataoccs , const Matd& datafree , const double& rad )
    {
        init( dataoccs , datafree , rad );
    }

    /// Return number of objects

    int size() const { return objs.size(); }

    /// Return object

    Group& operator[]( const int& i ) { return objs[i]; }
    const Group& operator[]( const int& i ) const { return objs[i]; }

    /// Delete cluster given index

    void del( const int& i )
    {
        objs.erase( objs.begin() + i );
    }

    /// Initializer for constructor

    void init( const Matd& dataoccs , const Matd& datafree , const double& rad )
    {
        this->rad = rad;

        SeqMatd P,M,S;

        randomise(4);
        quick_means( dataoccs , P , M , S , rad , 0.5 * rad , 3 );
        occs.set( P , M , S , +1.0 ); objs = detectObjects( occs , rad , 8.0 , 5 );
    }

    /// Return M for all groups

    Matd getM() const
    {
        int n = 0;
        forLOOPi( objs.size() ) n += objs[i].size();

        Matd M( n , 3 ); int cnt = 0;
        forLOOPij( objs.size() , objs[i].size() )
            M.row( cnt++ ) = objs[i][j].M.eig();

        return M;
    }

    void getData( SeqMatd& P , SeqMatd& M , SeqMatd& S , Matd& W )
    {
        int n = 0 , cnt = 0;
        forLOOPi( objs.size() ) n += objs[i].size();

        P.resize( n ); M.resize( n ); S.resize( n ); W.reset( n );
        forLOOPij( objs.size() , objs[i].size() )
        {
            P[ cnt ] = objs[i][j].P.clone();
            M[ cnt ] = objs[i][j].M.clone();
            S[ cnt ] = objs[i][j].S.clone();
            W( cnt++ ) = objs[i][j].wgt;
        }
    }

    /// Return O for all groups

    Mati getO() const
    {
        int n = 0;
        forLOOPi( objs.size() ) n += objs[i].size();

        Mati O( n ); int cnt = 0;
        forLOOPij( objs.size() , objs[i].size() )
            O( cnt++ ) = i + 1;

        return O;
    }

    /// Draw objects in screen

    void drawObjects( CPPlot& draw , const int& siz = 7 )
    {
        forLOOPi( objs.size() )
        {
            Matd M = objs[i].getM();
            if( M.filled() )
            {
                Scalar clr = BLA;//RGBint(i);
                draw.psc(siz,clr).pts3D( M );
            }
        }
    }

    void drawBoxes( CPPlot& draw , const int& thk = 2 )
    {
        forLOOPi( objs.size() )
        {
            Matd M = objs[i].getM();
            if( M.filled() )
            {
                Scalar clr = BLA;//RGBint(i);
                Matd Li = M.limRows( 1.5 * rad );

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

    void drawHistory( CPPlot& draw )
    {
        forLOOPi( objs.size() )
            draw.psc(3,BLA).line3D( objs[i].pos );
    }

};

#endif
