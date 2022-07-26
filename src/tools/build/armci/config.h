/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.ac by autoheader.  */

/* TODO */
/* #undef ACC_SMP */

/* TODO */
/* #undef ALLOC_MUNMAP */

/* TODO */
#define ALLOW_PIN 1

/* upper bound for ARMCI_DEFAULT_SHMMAX environment variable */
#define ARMCI_DEFAULT_SHMMAX_UBOUND 131072

/* Define to 1 if GPC calls are enabled */
#define ARMCI_ENABLE_GPC_CALLS 0

/* TODO */
/* #undef ARMCI_GROUP */

/* Define if profiling is enabled */
/* #undef ARMCI_PROFILE */

/* Define to 1 on BlueGene/L systems */
/* #undef BGL */

/* Define to 1 on BlueGene/P systems */
/* #undef BGP */

/* Define to name of caxpy routine to call from C */
#define BLAS_CAXPY caxpy

/* Define to name of ccopy routine to call from C */
#define BLAS_CCOPY ccopy

/* Define to name of daxpy routine to call from C */
#define BLAS_DAXPY daxpy

/* Define to name of dcopy routine to call from C */
#define BLAS_DCOPY dcopy

/* Define to name of saxpy routine to call from C */
#define BLAS_SAXPY saxpy

/* Define to name of scopy routine to call from C */
#define BLAS_SCOPY scopy

/* Define to sizeof(INTEGER) used to compile BLAS */
#define BLAS_SIZE 8

/* Define to name of zaxpy routine to call from C */
#define BLAS_ZAXPY zaxpy

/* Define to name of zcopy routine to call from C */
#define BLAS_ZCOPY zcopy

/* Define to 1 on Cray XT systems using Catamount */
/* #undef CATAMOUNT */

/* TODO */
/* #undef CLUSTER */

/* TODO */
/* #undef CRAY_SHMEM */

/* Define to 1 on Cray SV1 systems */
/* #undef CRAY_SV1 */

/* Define to 1 on Cray SV2 systems */
/* #undef CRAY_SV2 */

/* Define to 1 on Cray T3E systems */
/* #undef CRAY_T3E */

/* for Gemini */
/* #undef CRAY_UGNI */

/* TODO */
/* #undef CRAY_XT */

/* Define to 1 on Cray YMP systems */
/* #undef CRAY_YMP */

/* Define to 1 on Cygnus systems */
/* #undef CYGNUS */

/* Define to 1 on Cygwin systems */
/* #undef CYGWIN */

/* TODO */
#define DATA_SERVER 1

/* Define to 1 on DEC OSF */
/* #undef DECOSF */

/* TODO */
/* #undef DOELAN4 */

/* Define if checkpointing is enabled */
/* #undef ENABLE_CHECKPOINT */

/* Define if tracing is enabled */
/* #undef ENABLE_TRACE */

/* Define to 1 on fujitsu systems */
/* #undef FUJITSU_VPP */

/* Define to 1 on fujitsu systems */
/* #undef FUJITSU_VPP64 */

/* CPP symbol for function name, if available */
#define FUNCTION_NAME __func__

/* for Gemini */
/* #undef GEMINI */

/* TODO */
/* #undef GM */

/* set to 1 if ARMCI has armci_group_comm function */
#define HAVE_ARMCI_GROUP_COMM 1

/* set to 1 if ARMCI has ARMCI_Group.comm member */
/* #undef HAVE_ARMCI_GROUP_COMM_MEMBER */

/* set to 1 if ARMCI has ARMCI_Initialized function */
#define HAVE_ARMCI_INITIALIZED 1

/* set to 1 if ARMCI has armci_msg_finalize function */
#define HAVE_ARMCI_MSG_FINALIZE 1

/* set to 1 if ARMCI has armci_msg_init function */
#define HAVE_ARMCI_MSG_INIT 1

/* set to 1 if ARMCI has armci_notify function */
#define HAVE_ARMCI_NOTIFY 1

/* set to 1 if ARMCI has armci_stride_info_init function */
#define HAVE_ARMCI_STRIDE_INFO_INIT 1

/* Define to 1 if you have the <assert.h> header file. */
#define HAVE_ASSERT_H 1

/* Define to 1 if using external BLAS library */
#define HAVE_BLAS 1

/* Define to 1 if you have the `bzero' function. */
#define HAVE_BZERO 1

/* Define to 1 if you have the <c_asm.h> header file. */
#define HAVE_C_ASM_H 0

/* Define to 1 if you have the <dirent.h> header file, and it defines `DIR'.
   */
#define HAVE_DIRENT_H 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if you don't have `vprintf' but do have `_doprnt.' */
/* #undef HAVE_DOPRNT */

/* Define to 1 if you have the <errno.h> header file. */
#define HAVE_ERRNO_H 1

/* Define to 1 if you have the `fastbcopy' function. */
#define HAVE_FASTBCOPY 0

/* Define to 1 if you have the <fcntl.h> header file. */
#define HAVE_FCNTL_H 1

/* Define to 1 if you have the <float.h> header file. */
#define HAVE_FLOAT_H 1

/* Define to 1 if you have the 'floor' function. */
#define HAVE_FLOOR 1

/* Define to 1 if you have the `fork' function. */
#define HAVE_FORK 1

/* Define to 1 if you have the `ftruncate' function. */
#define HAVE_FTRUNCATE 1

/* Define to 1 if you have the `gethostbyname' function. */
#define HAVE_GETHOSTBYNAME 1

/* Define to 1 if you have the `getpagesize' function. */
#define HAVE_GETPAGESIZE 1

/* Define to 1 if you have the `gettimeofday' function. */
#define HAVE_GETTIMEOFDAY 1

/* Define to 1 if you have the <ia64/sys/inline.h> header file. */
#define HAVE_IA64_SYS_INLINE_H 0

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the <malloc.h> header file. */
#define HAVE_MALLOC_H 1

/* Define to 1 if you have the <math.h> header file. */
#define HAVE_MATH_H 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the `memset' function. */
#define HAVE_MEMSET 1

/* Define to 1 if you have a working `mmap' system call. */
#define HAVE_MMAP 1

/* Define to 1 if you have the <mpp/shmem.h> header file. */
#define HAVE_MPP_SHMEM_H 0

/* Define to 1 if you have the `munmap' function. */
#define HAVE_MUNMAP 1

/* Define to 1 if you have the <ndir.h> header file, and it defines `DIR'. */
/* #undef HAVE_NDIR_H */

/* Define to 1 if you have the <netdb.h> header file. */
#define HAVE_NETDB_H 1

/* Define to 1 if you have the <netinet/in.h> header file. */
#define HAVE_NETINET_IN_H 1

/* Define to 1 if you have the <netinet/tcp.h> header file. */
#define HAVE_NETINET_TCP_H 1

/* set to 1 if libonesided has onesided_fadd (added in v1.6) */
/* #undef HAVE_ONESIDED_FADD */

/* set to 1 if libonesided has onesided_mem_htflush (added in v1.5) */
/* #undef HAVE_ONESIDED_MEM_HTFLUSH */

/* Define to 1 if you have the 'pow' function. */
#define HAVE_POW 1

/* Define to 1 if you have the <process.h> header file. */
#define HAVE_PROCESS_H 0

/* Define to 1 if you have the `pthread_getconcurrency' function. */
#define HAVE_PTHREAD_GETCONCURRENCY 0

/* Define to 1 if you have the `pthread_setconcurrency' function. */
#define HAVE_PTHREAD_SETCONCURRENCY 0

/* Define to 1 if the system has the type `ptrdiff_t'. */
#define HAVE_PTRDIFF_T 1

/* Define to 1 if you have the <rpc/rpc.h> header file. */
#define HAVE_RPC_RPC_H 1

/* Define to 1 if you have the <rpc/types.h> header file. */
#define HAVE_RPC_TYPES_H 1

/* Define to 1 if you have the <rpc/xdr.h> header file. */
#define HAVE_RPC_XDR_H 1

/* Define to 1 if you have the `select' function. */
#define HAVE_SELECT 1

/* Define to 1 if you have the <setjmp.h> header file. */
#define HAVE_SETJMP_H 1

/* Define to 1 if you have the <signal.h> header file. */
#define HAVE_SIGNAL_H 1

/* Define to 1 if you have the `socket' function. */
#define HAVE_SOCKET 1

/* Define to 1 if you have the 'sqrt' function. */
#define HAVE_SQRT 1

/* Define to 1 if you have the <stdarg.h> header file. */
#define HAVE_STDARG_H 1

/* Define to 1 if stdbool.h conforms to C99. */
#define HAVE_STDBOOL_H 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdio.h> header file. */
#define HAVE_STDIO_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the `strchr' function. */
#define HAVE_STRCHR 1

/* Define to 1 if you have the `strdup' function. */
#define HAVE_STRDUP 1

/* Define to 1 if you have the `strerror' function. */
#define HAVE_STRERROR 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the `strstr' function. */
#define HAVE_STRSTR 1

/* Define to 1 if you have the `strtol' function. */
#define HAVE_STRTOL 1

/* Define to 1 if you have the <sys/atomic_op.h> header file. */
#define HAVE_SYS_ATOMIC_OP_H 0

/* Define to 1 if you have the <sys/dir.h> header file, and it defines `DIR'.
   */
/* #undef HAVE_SYS_DIR_H */

/* Define to 1 if you have the <sys/errno.h> header file. */
#define HAVE_SYS_ERRNO_H 1

/* Define to 1 if you have the <sys/file.h> header file. */
#define HAVE_SYS_FILE_H 1

/* Define to 1 if you have the <sys/ipc.h> header file. */
#define HAVE_SYS_IPC_H 1

/* Define to 1 if you have the <sys/mman.h> header file. */
#define HAVE_SYS_MMAN_H 1

/* Define to 1 if you have the <sys/ndir.h> header file, and it defines `DIR'.
   */
/* #undef HAVE_SYS_NDIR_H */

/* Define to 1 if you have the <sys/param.h> header file. */
#define HAVE_SYS_PARAM_H 1

/* Define to 1 if you have the <sys/select.h> header file. */
#define HAVE_SYS_SELECT_H 1

/* Define to 1 if you have the <sys/sem.h> header file. */
#define HAVE_SYS_SEM_H 1

/* Define to 1 if you have the <sys/shm.h> header file. */
#define HAVE_SYS_SHM_H 1

/* Define to 1 if you have the <sys/socket.h> header file. */
#define HAVE_SYS_SOCKET_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/syscall.h> header file. */
#define HAVE_SYS_SYSCALL_H 1

/* Define to 1 if you have the <sys/systemcfg.h> header file. */
#define HAVE_SYS_SYSTEMCFG_H 0

/* Define to 1 if you have the <sys/time.h> header file. */
#define HAVE_SYS_TIME_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <sys/uio.h> header file. */
#define HAVE_SYS_UIO_H 1

/* Define to 1 if you have the <sys/wait.h> header file. */
#define HAVE_SYS_WAIT_H 1

/* Define this if weak aliases may be created with #pragma weak */
#define HAVE_SYS_WEAK_ALIAS_PRAGMA 1

/* Define to 1 if you have the <time.h> header file. */
#define HAVE_TIME_H 1

/* define if sys/sem.h has union semun */
/* #undef HAVE_UNION_SEMUN */

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if you have the `vfork' function. */
#define HAVE_VFORK 1

/* Define to 1 if you have the <vfork.h> header file. */
/* #undef HAVE_VFORK_H */

/* Define to 1 if you have the `vprintf' function. */
#define HAVE_VPRINTF 1

/* Define to 1 if you have the <windows.h> header file. */
#define HAVE_WINDOWS_H 0

/* Define to 1 if you have the <winsock.h> header file. */
#define HAVE_WINSOCK_H 0

/* Define to 1 if `fork' works. */
#define HAVE_WORKING_FORK 1

/* Define to 1 if `vfork' works. */
#define HAVE_WORKING_VFORK 1

/* Define to 1 if you have the `xdr_char' function. */
#define HAVE_XDR_CHAR 1

/* Define to 1 if you have the `_acquire_lock' function. */
#define HAVE__ACQUIRE_LOCK 0

/* Define to 1 if the system has the type `_Bool'. */
#define HAVE__BOOL 1

/* Define to 1 if you have the `_lock_try' function. */
#define HAVE__LOCK_TRY 0

/* Define to 1 on hitachi systems */
/* #undef HITACHI */

/* Define to 1 on HP-UX systems */
/* #undef HPUX */

/* Define to 1 on 64bit HP-UX systems */
/* #undef HPUX64 */

/* Define to 1 on IBM SP systems */
/* #undef IBM */

/* Define to 1 on 64bit IBM SP systems */
/* #undef IBM64 */

/* Define to 1 on ??? systems */
/* #undef INTERIX */

/* TODO */
/* #undef LAPI */

/* Define to 1 on 64bit IBM systems with LAPI */
/* #undef LAPI64 */

/* TODO */
/* #undef LAPI_RDMA */

/* for Gemini */
/* #undef LIBONESIDED */

/* Define to 1 on generic Linux systems */
#define LINUX 1

/* Define to 1 on generic 64bit Linux systems */
#define LINUX64 1

/* Define to the sub-directory where libtool stores uninstalled libraries. */
#define LT_OBJDIR ".libs/"

/* Define to 1 on OSX systems */
/* #undef MACX */

/* Define to 1 on 64bit OSX systems */
/* #undef MACX64 */

/* TODO */
/* #undef MELLANOX */

/* TODO */
/* #undef MPI3 */

/* TODO */
/* #undef MPI_MT */

/* TODO */
/* #undef MPI_PR */

/* TODO */
/* #undef MPI_PT */

/* TODO */
/* #undef MPI_SPAWN */

/* TODO */
/* #undef MPI_TS */

/* Use MPI for messaging */
#define MSG_COMMS_MPI 1

/* Use TCGMSG for messaging */
/* #undef MSG_COMMS_TCGMSG */

/* Use TCGMSG (ipcv4.0) for messaging */
/* #undef MSG_COMMS_TCGMSG4 */

/* Use TCGMSG (ipcv5.0) for messaing */
/* #undef MSG_COMMS_TCGMSG5 */

/* TODO */
/* #undef MULTI_CTX */

/* TODO */
/* #undef NB_NONCONT */

/* Define to 1 if assertions should be disabled. */
/* #undef NDEBUG */

/* Define to 1 on NEC systems */
/* #undef NEC */

/* Define to 1 on 64bit NEC systems */
/* #undef NEC64 */

/* whether to wait until the last moment to call ARMCI_Init() in TCGMSG-MPI */
#define NEED_DELAY_TCGMSG_MPI_STARTUP 1

/* Creates memfenc macro */
#define NEED_MEM_SYNC 1

/* for portals, enable new malloc */
/* #undef NEW_MALLOC */

/* this was defined unconditionally when using LAPI for tcgmsg 5 */
/* #undef NOTIFY_SENDER */

/* define when inline asm is not supported */
/* #undef NO_I386ASM */

/* TODO */
#define OPENIB 1

/* Name of package */
#define PACKAGE "armci"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "hpctools@pnl.gov"

/* Define to the full name of this package. */
#define PACKAGE_NAME "Aggregate Remote Memory Copy Interface (ARMCI)"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "Aggregate Remote Memory Copy Interface (ARMCI) 1.5"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "armci"

/* Define to the home page for this package. */
#define PACKAGE_URL "http://www.emsl.pnl.gov/docs/parsoft/armci/"

/* Define to the version of this package. */
#define PACKAGE_VERSION "1.5"

/* TODO */
#define PEND_BUFS 1

/* TODO */
/* #undef PORTALS */

/* TODO */
/* #undef POSIX_THREADS */

/* TODO */
#define PTHREADS 1

/* TODO */
/* #undef QUADRICS */

/* TODO */
/* #undef REPORT_SHMMAX */

/* Define as the return type of signal handlers (`int' or `void'). */
#define RETSIGTYPE void

/* Define to the type of arg 1 for `select'. */
#define SELECT_TYPE_ARG1 int

/* Define to the type of args 2, 3 and 4 for `select'. */
#define SELECT_TYPE_ARG234 (fd_set *)

/* Define to the type of arg 5 for `select'. */
#define SELECT_TYPE_ARG5 (struct timeval *)

/* TODO */
#define SERVER_THREAD 1

/* Define to 1 on ??? systems */
/* #undef SGI */

/* Define to 1 on ??? systems */
/* #undef SGITFP */

/* Define to 1 on ??? systems */
/* #undef SGI_N32 */

/* The size of `char', as computed by sizeof. */
#define SIZEOF_CHAR 1

/* The size of `double', as computed by sizeof. */
#define SIZEOF_DOUBLE 8

/* The size of `float', as computed by sizeof. */
#define SIZEOF_FLOAT 4

/* The size of `int', as computed by sizeof. */
#define SIZEOF_INT 4

/* The size of `long', as computed by sizeof. */
#define SIZEOF_LONG 8

/* The size of `long long', as computed by sizeof. */
#define SIZEOF_LONG_LONG 8

/* The size of `short', as computed by sizeof. */
#define SIZEOF_SHORT 2

/* The size of `void*', as computed by sizeof. */
#define SIZEOF_VOIDP 8

/* TODO */
/* #undef SOCKETS */

/* Define to 1 on Solaris systems */
/* #undef SOLARIS */

/* Define to 1 on 64bit Solaris systems */
/* #undef SOLARIS64 */

/* TODO */
/* #undef SPAWN_CRAY_XT */

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Defined if we want this system to use SYSV shared memory */
#define SYSV 1

/* deprecated, use MSG_COMMS_TCGMSG */
/* #undef TCGMSG */

/* remote shell for TCGMSG */
#define TCGMSG_RSH "/usr/bin/rsh"

/* Gather timing information for TCGMSG */
#define TCGMSG_TIMINGS 1

/* turn on thread safety */
/* #undef THREAD_SAFE */

/* Enable extensions on AIX 3, Interix.  */
#ifndef _ALL_SOURCE
# define _ALL_SOURCE 1
#endif
/* Enable GNU extensions on systems that have them.  */
#ifndef _GNU_SOURCE
# define _GNU_SOURCE 1
#endif
/* Enable threading extensions on Solaris.  */
#ifndef _POSIX_PTHREAD_SEMANTICS
# define _POSIX_PTHREAD_SEMANTICS 1
#endif
/* Enable extensions on HP NonStop.  */
#ifndef _TANDEM_SOURCE
# define _TANDEM_SOURCE 1
#endif
/* Enable general extensions on Solaris.  */
#ifndef __EXTENSIONS__
# define __EXTENSIONS__ 1
#endif


/* TODO */
#define VAPI 1

/* Version number of package */
#define VERSION "1.5"

/* TODO */
/* #undef VIA */

/* Define to 1 if on MINIX. */
/* #undef _MINIX */

/* Define to 2 if the system does not provide POSIX.1 features except with
   this defined. */
/* #undef _POSIX_1_SOURCE */

/* Define to 1 if you need to in order for `stat' and other things to work. */
/* #undef _POSIX_SOURCE */

/* TODO */
#define _REENTRANT 1

/* Define for Solaris 2.5.1 so the uint32_t typedef from <sys/synch.h>,
   <pthread.h>, or <semaphore.h> is not used. If the typedef were allowed, the
   #define below would cause a syntax error. */
/* #undef _UINT32_T */

/* Define for Solaris 2.5.1 so the uint64_t typedef from <sys/synch.h>,
   <pthread.h>, or <semaphore.h> is not used. If the typedef were allowed, the
   #define below would cause a syntax error. */
/* #undef _UINT64_T */

/* Define for Solaris 2.5.1 so the uint8_t typedef from <sys/synch.h>,
   <pthread.h>, or <semaphore.h> is not used. If the typedef were allowed, the
   #define below would cause a syntax error. */
/* #undef _UINT8_T */

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif

/* Define to the type of a signed integer type of width exactly 16 bits if
   such a type exists and the standard includes do not define it. */
/* #undef int16_t */

/* Define to the type of a signed integer type of width exactly 32 bits if
   such a type exists and the standard includes do not define it. */
/* #undef int32_t */

/* Define to the type of a signed integer type of width exactly 64 bits if
   such a type exists and the standard includes do not define it. */
/* #undef int64_t */

/* Define to `long int' if <sys/types.h> does not define. */
/* #undef off_t */

/* Define to `int' if <sys/types.h> does not define. */
/* #undef pid_t */

/* Define to the equivalent of the C99 'restrict' keyword, or to
   nothing if this is not supported.  Do not define if restrict is
   supported directly.  */
#define restrict __restrict
/* Work around a bug in Sun C++: it does not support _Restrict or
   __restrict__, even though the corresponding Sun C compiler ends up with
   "#define restrict _Restrict" or "#define restrict __restrict__" in the
   previous line.  Perhaps some future version of Sun C++ will work with
   restrict; if so, hopefully it defines __RESTRICT like Sun C does.  */
#if defined __SUNPRO_CC && !defined __RESTRICT
# define _Restrict
# define __restrict__
#endif

/* Define to `unsigned int' if <sys/types.h> does not define. */
/* #undef size_t */

/* Define to `int' if <sys/types.h> does not define. */
/* #undef ssize_t */

/* Define to the type of an unsigned integer type of width exactly 16 bits if
   such a type exists and the standard includes do not define it. */
/* #undef uint16_t */

/* Define to the type of an unsigned integer type of width exactly 32 bits if
   such a type exists and the standard includes do not define it. */
/* #undef uint32_t */

/* Define to the type of an unsigned integer type of width exactly 64 bits if
   such a type exists and the standard includes do not define it. */
/* #undef uint64_t */

/* Define to the type of an unsigned integer type of width exactly 8 bits if
   such a type exists and the standard includes do not define it. */
/* #undef uint8_t */

/* Define as `fork' if `vfork' does not work. */
/* #undef vfork */

/* Define to empty if the keyword `volatile' does not work. Warning: valid
   code using `volatile' can become incorrect without. Disable with care. */
/* #undef volatile */
