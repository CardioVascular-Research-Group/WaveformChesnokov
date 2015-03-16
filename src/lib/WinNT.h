/* contains a subset of the lines found in c:\Program Files (x86)\Microsoft SDKs\Windows\v7.0A\Include\WinNT.h
 * which are needed to compile the original MS Visual C++ version of Chesnokov.
 * Ideally, all the lines here can be removed as Microsoft specific code is replaced.
 * WinNT.h
 *
 *  Created on: Mar 9, 2015
 *      Author: root
 */

#ifndef WINNT_H_
#define WINNT_H_

//
// Void
//

//typedef void *PVOID;
//typedef void * POINTER_64;
//typedef void * PVOID64;

//
// Handle to an Object
//

//#ifdef STRICT
	typedef void *HANDLE;
	#if 0 && (_MSC_VER > 1000)
		#define DECLARE_HANDLE(name) struct name##__; typedef struct name##__ *name
	#else
		#define DECLARE_HANDLE(name) struct name##__{int unused;}; typedef struct name##__ *name
	#endif
//#else
//	typedef PVOID HANDLE;
//	#define DECLARE_HANDLE(name) typedef HANDLE name
//#endif
typedef HANDLE *PHANDLE;



#endif /* WINNT_H_ */
