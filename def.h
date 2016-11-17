//-----------------------------------------------------------------------------------
//
//	Basic definitions
//
//	Anton M. Krivtsov
//
//	13.06.2001
//
//-----------------------------------------------------------------------------------

#ifndef __DEF_H__
#define __DEF_H__

//-----------------------------------------------------------------------------------

//#define FLT double
#define FLT float
#define DBL double
#define FLA const FLT&

//-----------------------------------------------------------------------------------
//	Definitions for compartibility with MFC
//-----------------------------------------------------------------------------------

#define LONG  int
#define SHORT int
#define BYTE  int
#define BOOL  int
#define afxDump cout
#define TRUE  1
#define FALSE 0
#define NULL  0

//#define min(a, b)  (((a) < (b)) ? (a) : (b))
//#define max(a, b)  (((a) > (b)) ? (a) : (b))

//-----------------------------------------------------------------------------------
//	Definitions for compartibility between DOS and UNIX versions of C++
//-----------------------------------------------------------------------------------

//	#define UNIX	// Use this definition for compilation under UNIX

//-----------------------------------------------------------------------------------

#endif //__DEF_H__

