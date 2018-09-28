/*
 * util.h
 */

#ifndef UTIL_H
#define UTIL_H

#include <GASPI.h>
#include <stddef.h> /* size_t */
#include <stdio.h> /* printf */
#include <stdlib.h> /* EXIT_FAILURE */
#include <float.h> /* DBL_EPSILON */

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

#if defined(DEBUG)
#  define DBG_MSG printf
#  define DBG_PRINT 1
#else
#  define DBG_MSG(...)
#  define DBG_PRINT 0
#endif

#if defined SCHEDULE_DYNAMIC
#  if defined SET_CHUNKSIZE
#    define SCHEDULE schedule(dynamic, SET_CHUNKSIZE)
#  else
#    define SCHEDULE schedule(dynamic)
#  endif
#elif defined SCHEDULE_GUIDED
#  if defined SET_CHUNKSIZE
#    define SCHEDULE schedule(guided, SET_CHUNKSIZE)
#  else
#    define SCHEDULE schedule(guided)
#  endif
#elif defined SCHEDULE_AUTO
# define SCHEDULE schedule(auto)
#elif defined SCHEDULE_RUNTIME
# define SCHEDULE schedule(runtime)
#else
#  if defined SET_CHUNKSIZE
#    define SCHEDULE schedule(static, SET_CHUNKSIZE)
#  else
#    define SCHEDULE schedule(static)
#  endif
#endif

#define STR(token) #token
#define EXPSTR(macro) STR(macro)
#define SCHEDSTRING EXPSTR(SCHEDULE)

#define CHECK(expr) \
  if(!(expr)) \
  { \
    printf("Error: '%s' [%s:%i]\n", #expr, __FILE__, __LINE__); \
    exit(EXIT_FAILURE); \
  }

#define GASPICHECK(f...) \
  do \
  { \
    const gaspi_return_t r = f; \
    if (r != GASPI_SUCCESS) \
    { \
      printf("Error: '%s' [%s:%i]: %i\n", #f, __FILE__, __LINE__, r); \
      exit(EXIT_FAILURE); \
    }\
  } while(0)

/*******************************************************************************
*
*******************************************************************************/
void  check_free(void *ptr);
void *check_malloc(size_t bytes);
void *check_calloc(size_t number, size_t bytes);
void *check_realloc(void *old, size_t bytes);

#endif /* UTIL_H */
