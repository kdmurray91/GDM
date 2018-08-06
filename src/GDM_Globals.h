//
//
// GDM_Globals.h
//
//
#ifndef __GDM_GLOBALS_H__
#define __GDM_GLOBALS_H__


///////////////////////////////////////////////////////////////////////////////////
//
//           NOTE: Uncomment the following line if compiling for Windows
//
//#define WINDOWS
//
//              otherwise uncomment the following line if compiling for UNIX
//
//#undef WINDOWS
//
//////////////////////////////////////////////////////////////////////////////////

#ifdef WINDOWS
	#include <stdio.h>
	#include <stdlib.h>
	#include <io.h>
	#include <fcntl.h>
	#include <sys/stat.h>
	#include <math.h>	
	#include <float.h>
#else
	#include <stdio.h>
	#include <unistd.h>
	#include <stdlib.h>
	#include <math.h>
	#include <time.h>
	//#include <sys/io.h>
	#include <string.h>
	#include <fcntl.h>
	#include <sys/stat.h>
	#include <float.h>
	#include <limits.h>
#endif


#endif // __GDM_GLOBALS_H__
