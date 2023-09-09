/*
 * crypto_core/try-notest.c version 20180223
 * D. J. Bernstein
 * Public domain.
 * Auto-generated by trygen-notest.py; do not edit.
 */

#include "crypto_core.h"
#include "try.h"

const char *primitiveimplementation = crypto_core_implementation;

#ifdef SMALL
#define LOOPS 512
#else
#define LOOPS 4096
#endif

static unsigned char *h;
static unsigned char *n;
static unsigned char *k;
static unsigned char *c;
static unsigned char *h2;
static unsigned char *n2;
static unsigned char *k2;
static unsigned char *c2;
#define hlen crypto_core_OUTPUTBYTES
#define nlen crypto_core_INPUTBYTES
#define klen crypto_core_KEYBYTES
#define clen crypto_core_CONSTBYTES

void preallocate(void)
{
}

void allocate(void)
{
  unsigned long long alloclen = 0;
  if (alloclen < crypto_core_OUTPUTBYTES) alloclen = crypto_core_OUTPUTBYTES;
  if (alloclen < crypto_core_INPUTBYTES) alloclen = crypto_core_INPUTBYTES;
  if (alloclen < crypto_core_KEYBYTES) alloclen = crypto_core_KEYBYTES;
  if (alloclen < crypto_core_CONSTBYTES) alloclen = crypto_core_CONSTBYTES;
  h = alignedcalloc(alloclen);
  n = alignedcalloc(alloclen);
  k = alignedcalloc(alloclen);
  c = alignedcalloc(alloclen);
  h2 = alignedcalloc(alloclen);
  n2 = alignedcalloc(alloclen);
  k2 = alignedcalloc(alloclen);
  c2 = alignedcalloc(alloclen);
}

void predoit(void)
{
}

void doit(void)
{
  crypto_core(h,n,k,c);
}