
#include <cvpp/interfaces/cpplot.h>
#include <cvpp/modules/algl.h>

using namespace cvpp;

int main()
{

    Matd o( 1 , 3 ) , v( 1 , 3 ) , w( 1 , 3 );
    o.eig() << 0.0 , 0.0 , 0.0;
    v.eig() << 4.0 , 2.0 , 5.0;
    w.eig() << 3.0 , 0.0 , 0.0;

    Matd R = algl::rotationAlign( v , w );

    Matd C( 3 , 3 );
    C.eig() << 1.0 , 0.0 , 0.0 ,
               0.0 , 0.2 , 0.0 ,
               0.0 , 0.0 , 0.2 ;

    w = w * R;
    C = R.t() * C * R;

    CPPlot draw( "Window" );
    draw[0].set3Dworld();

    while( draw.input() )
    {
        draw[0].clear();
        draw.lwc(3,YEL).line3D( o | v );
        draw.lwc(6,RED).line3D( o | w );
        draw.ellipse3D( o , C );
        draw.updateWindow(30);
    }

    return 0;
}
