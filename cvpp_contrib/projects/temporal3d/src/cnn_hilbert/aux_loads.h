
#include <dirent.h>
#include <sys/types.h>
#include <cvpp/containers/matrix.h>
#include <cvpp/containers/vector.h>
#include <cvpp/containers/image.h>
#include <cvpp/properties/pose.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

using namespace cvpp;

Seq<String>
get_files( const String& dir , const int& n = 0 )
{
    DIR *dp;
    struct dirent *dirp;
    Seq<String> files;

    if( ( dp = opendir( dir.c_str() ) ) == NULL )
        disp( "Error Opening" , dir );

    while( ( dirp = readdir( dp ) ) != NULL)
    {
        String file( dirp->d_name );
        if( file[ file.size() - 4 ] == '.' )
            files.push_back( dir + file );
    }

    closedir( dp );

    std::sort( files.begin() , files.end() );
    if( n > 0 ) files.resize( n );

    return files;
}

Matf
load_vel2cam( const String& file )
{
    String line;
    std::ifstream infile( file + "/calib_velo_to_cam.txt" );

    float R[9] , t[3];
    while( std::getline( infile , line ) )
    {
        if( line[0] == 'R' && line[1] == ':' )
            tokenFloat( line.c_str() , R , ' ' );

        if( line[0] == 'T' && line[1] == ':' )
            tokenFloat( line.c_str() , t , ' ' );
    }

    Matf T( 4 , 4 );
    T.eig() << R[0] , R[1] , R[2] , t[0] ,
               R[3] , R[4] , R[5] , t[1] ,
               R[6] , R[7] , R[8] , t[2] ,
                0.0 ,  0.0 ,  0.0 ,  1.0 ;

    return T.t();
}

Matf
load_imu2vel( const String& file )
{
    String line;
    std::ifstream infile( file + "/calib_imu_to_velo.txt" );

    float R[9] , t[3];
    while( std::getline( infile , line ) )
    {
        if( line[0] == 'R' && line[1] == ':' )
            tokenFloat( line.c_str() , R , ' ' );

        if( line[0] == 'T' && line[1] == ':' )
            tokenFloat( line.c_str() , t , ' ' );
    }

    Matf T( 4 , 4 );
    T.eig() << R[0] , R[1] , R[2] , t[0] ,
               R[3] , R[4] , R[5] , t[1] ,
               R[6] , R[7] , R[8] , t[2] ,
                0.0 ,  0.0 ,  0.0 ,  1.0 ;

    return T.t();
}

void
load_cam2cam( const String& dir , Matf& K , Matf& D , Matf& R , Matf& P )
{
    String file = dir + "/calib_cam_to_cam.txt";

    String line;
    std::ifstream infile( file );

    float k[9] , d[5] , r[9] , p[12];
    while( std::getline( infile , line ) )
    {
        if( line.substr(0,4).compare( "K_02" ) == 0 )
            tokenFloat( line.substr(5).c_str() , k , ' ' );

        if( line.substr(0,4).compare( "D_02" ) == 0 )
            tokenFloat( line.substr(5).c_str() , d , ' ' );

        if( line.substr(0,9).compare( "R_rect_00" ) == 0 )
            tokenFloat( line.substr(10).c_str() , r , ' ' );

        if( line.substr(0,9).compare( "P_rect_02" ) == 0 )
            tokenFloat( line.substr(10).c_str() , p , ' ' );
    }

    K.reset( 3 , 3 );
    forLOOPij( K.r() , K.c() )
        K(i,j) = k[ i * K.c() + j ];

    D.reset( 5 );
    forLOOPi( D.r() )
        D(i) = d[ i ];

    R.reset( 4 , 4 ).setIdentity();
    forLOOPij( 3 , 3 )
        R(i,j) = r[ i * 3 + j ];
    R.blu(3) = R.blu(3).t();

    P.reset( 3 , 4 );
    forLOOPij( P.r() , P.c() )
        P(i,j) = p[ i * P.c() + j ];
    P = P.t();
}

SeqMatd
load_vel( const Seq<String>& files )
{
    int n = files.size();
    SeqMatd vels( n );

    int base = 1000000;
    float *data = (float*)malloc( base * sizeof(float) );

    forLOOPi( n )
    {
        float *px = data + 0 , *py = data + 1;
        float *pz = data + 2 , *pr = data + 3;

        FILE *stream;
        stream = fopen( files[i].c_str() , "rb" );
        int num = fread( data , sizeof(float) , base , stream ) / 4;
        vels[i].reset( num , 4 );

        forLOOPj( num )
        {
            vels[i].row(j) << *px , *py , *pz , *pr ;
            px += 4 ; py += 4 ; pz += 4 ; pr += 4 ;
        }

        fclose( stream );
    }

    return vels;
}

Matd load_vel( const String& file )
{
    Seq<String> files;
    files.push_back( file );
    SeqMatd vels = load_vel( files );
    return vels[0];
}

SeqImg3c
load_img( const Seq<String>& files )
{
    int n = files.size();
    SeqImg3c imgs( n );

    #pragma omp parallel for
    forLOOPi( n )
    {
        imgs[i].load( files[i] );
    }

    return imgs;
}

Img3c load_img( const String& file )
{
    Seq<String> files;
    files.push_back( file );
    SeqImg3c imgs = load_img( files );
    return imgs[0];
}

SeqPosef
load_pos( const Seq<String>& files )
{
    int n = files.size();

    Matf data( n , 30 );

    forLOOPi( n )
    {
        float vals[30];
        std::ifstream infile( files[i] );

        String line;
        while( std::getline( infile , line ) )
            tokenFloat( ( ' ' + line ).c_str() , vals , ' ' );

        forLOOPj( data.c() )
            data(i,j) = vals[j];
    }

    float lat0 = data(0,0);
    float r = 6378137 , s = std::cos( lat0 * PI / 180.0 );
    float sr = s * r;

    Matf xyz( n , 6 );
    forLOOPi( xyz.r() )
    {
        float lat = data(i,0) , lon = data(i,1);
        float z = data(i,2) , r = data(i,3) , p = data(i,4) , w = data(i,5);
        float x = sr * PI * lon / 180.0;
        float y = sr * std::log( std::tan( PI * ( 90.0 + lat ) / 360.0 ) );

        xyz.row(i) << x , y , z , r , p , w;


    }

    Matf off = xyz.cl(3).r(0).clone();
    xyz.cl(3) -= off;

    SeqPosef poses( n );
    forLOOPi( poses.size() )
        poses[i].setPose( xyz.r(i) );

    return poses;
}

Posef load_pos( const String& file )
{
    Seq<String> files;
    files.push_back( file );
    SeqPosef poses = load_pos( files );
    return poses[0];
}

void prep_dirs( const String& dir )
{
    String path;
    struct stat st = {0};

    path = dir ;           if( stat( path.c_str() , &st ) == -1 ) mkdir( path.c_str() , 0700 );
    path = dir + "/imgs";  if( stat( path.c_str() , &st ) == -1 ) mkdir( path.c_str() , 0700 );
    path = dir + "/disp1"; if( stat( path.c_str() , &st ) == -1 ) mkdir( path.c_str() , 0700 );
    path = dir + "/disp2"; if( stat( path.c_str() , &st ) == -1 ) mkdir( path.c_str() , 0700 );
    path = dir + "/xyz";   if( stat( path.c_str() , &st ) == -1 ) mkdir( path.c_str() , 0700 );
    path = dir + "/all";   if( stat( path.c_str() , &st ) == -1 ) mkdir( path.c_str() , 0700 );
    path = dir + "/uvxyz"; if( stat( path.c_str() , &st ) == -1 ) mkdir( path.c_str() , 0700 );
}

void save_data( const String& suf , const String& dir , const String& path ,
                const Img3c& img , const Img1c& dsp1 , const Img1c& dsp2 ,
                const Matf& xyz , const Matf& uvxyz , const Matf& all )
{
    int nn = path.length() , n = 0;
    while( path[ nn - n ] != '/' ) n++; n--;
    String name = path.substr( nn - n , n -4 );

    String sufname = suf + name;
    disp( sufname );

    img.saveIMG(   dir + "/imgs/"  + sufname + ".png" );
    dsp1.saveIMG(  dir + "/disp1/" + sufname + ".png" );
    dsp2.saveIMG(  dir + "/disp2/" + sufname + ".png" );
    xyz.saveBIN(   dir + "/xyz/"   + sufname + ".bin" );
    all.saveBIN(   dir + "/all/"   + sufname + ".bin" );
    uvxyz.saveBIN( dir + "/uvxyz/" + sufname + ".bin" );
}

