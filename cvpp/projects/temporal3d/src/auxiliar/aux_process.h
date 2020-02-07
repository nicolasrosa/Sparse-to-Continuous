
#include <cvpp/containers/matrix.h>

using namespace cvpp;

// Load KITTI Pose
Matd
loadKITTIpos( const String& file , const int& n )
{
    String pos_file = file + "/oxts/data/";

    Matd pos( n , 30 );
    forLOOPi( n )
    {
        String suf = toString( i , 10 ) + ".txt";
        std::ifstream infile( pos_file + suf );

        String line; float vals[30];
        while( std::getline( infile , line ) )
            tokenFloat( ( ' ' + line ).c_str() , vals , ' ' );

        forLOOPj( pos.c() )
            pos(i,j) = vals[j];
    }

    float lat0 = pos(0,0);
    float r = 6378137 , s = std::cos( lat0 * PI / 180.0 );
    float sr = s * r;

    Matd xyz( n , 6 );
    forLOOPi( xyz.r() )
    {
        float lat = pos(i,0) , lon = pos(i,1);
        float z = pos(i,2) , r = pos(i,3) , p = pos(i,4) , w = pos(i,5);
        float x = sr * PI * lon / 180.0;
        float y = sr * std::log( std::tan( PI * ( 90.0 + lat ) / 360.0 ) );

        xyz.row(i) << x , y , z , r , p , w;
    }

    Matd off = xyz.cl(3).r(0).clone();
    xyz.cl(3) -= off;

    return xyz;
}

// Load KITTI velodyne
SeqMatd
loadKITTIvel( const String& file , const int& n )
{
    String vel_file = file + "/velodyne_points/data/";

    int base = 1000000;
    float *data = (float*)malloc( base * sizeof(float) );

    SeqMatd vel( n );
    forLOOPi( n )
    {
        String suf = toString( i , 10 ) + ".bin";

        float *px = data + 0 , *py = data + 1 , *pz = data + 2 , *pr = data + 3;

        FILE *stream;
        stream = fopen( ( vel_file + suf ).c_str() , "rb" );
        int num = fread( data , sizeof(float) , base , stream ) / 4;
        vel[i].reset( num , 4 );

        forLOOPj( num )
        {
            vel[i].row(j) << *px , *py , *pz , *pr ;
            px += 4 ; py += 4 ; pz += 4 ; pr += 4 ;
        }

        fclose( stream );
    }

    return vel;
}

Matd
loadKITTIvel2cam( const String& file )
{
    String line;
    std::ifstream infile( file + "/calib/calib_velo_to_cam.txt" );

    float R[9] , t[3];
    while( std::getline( infile , line ) )
    {
        if( line[0] == 'R' && line[1] == ':' )
            tokenFloat( line.c_str() , R , ' ' );

        if( line[0] == 'T' && line[1] == ':' )
            tokenFloat( line.c_str() , t , ' ' );
    }

    Matd T( 4 , 4);
    T.eig() << R[0] , R[1] , R[2] , t[0] ,
               R[3] , R[4] , R[5] , t[1] ,
               R[6] , R[7] , R[8] , t[2] ,
                0.0 ,  0.0 ,  0.0 ,  1.0 ;

    return T;
}

