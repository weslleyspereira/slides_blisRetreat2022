#include <cblas.h>
#include <math.h>
#include <float.h>
#include <complex.h>

int main() {

    {
        const float x[] = { 0, NAN, 2 };
        printf( "cblas_isamax([ 0, NaN, 2 ]) = %d\n", (int) cblas_isamax( 3, x, 1 ) );
        const float y[] = { NAN, 0, 2 };
        printf( "cblas_isamax([ NaN, 0, 2 ]) = %d\n", (int) cblas_isamax( 3, y, 1 ) );
    }
    {
        const float OV = FLT_MAX;
        const float _Complex x[] = { OV+OV*I, INFINITY };
        printf( "cblas_icamax([ OV+i*OV, Inf+i*0 ]) = %d\n", (int) cblas_icamax( 2, x, 1 ) );
        const float _Complex y[] = { .6*OV+.6*OV*I, .7*OV+.7*OV*I };
        printf( "cblas_icamax([ .6*OV+i*.6*OV, .7*OV+i*.7*OV ]) = %d\n", (int) cblas_icamax( 2, y, 1 ) );
    }
    {
        float A[] = { NAN, 1, 0, 2 };
        float b[] = { 0, 1 };
        
        int ipiv = cblas_isamax( 2, A, 1 );
        if( ipiv != 0 )
            cblas_sswap( 2, &A[0], 2, &A[1], 2 );
        cblas_sscal( 1, 1/A[0], &A[1], 1 );
        printf( "PA = [ %f, %f; %f, %f ]\n", A[0], A[1], A[2], A[3] );

        cblas_sger( CblasColMajor, 1, 1, -1, &A[1], 1, &A[2], 2, &A[3], 2 );
        printf( "L = [ 1, 0; %f, 1 ]\n", A[1] );
        printf( "U = [ %f, %f; 0, %f ]\n", A[0], A[2], A[3] );
        printf( "ipiv = %d\n", ipiv );

        if( ipiv != 0 )
            cblas_sswap( 1, &b[0], 1, &b[1], 1 );

        cblas_strsm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, 2, 2, 1, A, 2, b, 2);
        printf( "y = [ %f; %f ]\n", b[0], b[1] );

        cblas_strsm( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, 2, 2, 1, A, 2, b, 2);
        printf( "x = [ %f; %f ]\n", b[0], b[1] );
    }

    return 0;
}