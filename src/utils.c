#include "fid.h"

void fatal(const char * format, ...)
{
  va_list argptr;
  va_start(argptr, format);
  vfprintf(stderr, format, argptr);
  va_end(argptr);
  fprintf(stderr, "\n");
  exit(1);
}

void * xmalloc(size_t size, const size_t alignment)
{
  void * t = NULL;
  posix_memalign(&t, alignment, size);

  if (t==NULL)
    fatal("Unable to allocate enough memory.");

  return t;
}
