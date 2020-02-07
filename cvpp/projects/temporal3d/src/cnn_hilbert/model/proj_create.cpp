
#include <cvpp/interfaces/cpplot.h>

using namespace cvpp;

int main()
{
    randomise( 2 );

    int n = 3;

    Img3u img;
    STD<Img1u> gray( n );

    img.load( "../data/black_up.png" );
    gray[0] = img.grayscale();
    img.load( "../data/black_down.png" );
    gray[1] = img.grayscale();
    img.load( "../data/black_both.png" );
    gray[2] = img.grayscale();

    Matf mat( n , gray[0].s() );
    forLOOPi( gray.size() )
    {
        Matf tmp = gray[i].toFloat() / 255 - 0.5 ;
        mat.r(i) = tmp.toRow();
    }

    int m = 1000;
    Matf images( m * 3 , mat.c() );
    Matf labels( m * 3 , n ); labels.setVal(0);

    forLOOPi( m )
    {
        images.r( n * i , n ) = mat.addRand( 0.8 );
        forLOOPj( n ) labels( n * i + j , j ) = 1;
    }

    images.save("images_tst");
    labels.save("labels_tst");

//    forLOOPi( 100 )
//    {
//        Matf rec = ( images.r(i).reshape( 50 , 50 ) + 0.5 ) * 255.0;
//        rec.Floor( 0.0 ).Ceil( 255.0 );

//        Img1u out = rec.toUChar();
//        out.save("output" + toString(i) + ".png");
//    }

    return 0;
}






