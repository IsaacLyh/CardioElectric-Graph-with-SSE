//
// Timer definition
// Make successive calls and take a difference to get the elapsed time.
//

#include <cstdlib>
#include <iostream>
using namespace std;

#include <sys/time.h>
const double kMicro = 1.0e-6;
double getTime()
{
    struct timeval TV;

    const int RC = gettimeofday(&TV, NULL);
    if(RC == -1)
    {
        cout << "ERROR: Bad call to gettimeofday\n";
        return(-1);
    }

    return( ((double)TV.tv_sec) + kMicro * ((double)TV.tv_usec) );

}  // end getTime()
