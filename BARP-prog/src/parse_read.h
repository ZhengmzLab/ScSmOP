#ifndef _PARSE_READ_H_
#define _PARSE_READ_H_

#include "zlib.h"
#include "kstring.h"
#include "kseq.h"
#include "unistd.h"

KSEQ_INIT(gzFile, gzread)

typedef struct fastqhandler
{
    char * rd1, *rd2, *rd3, *rd4;
    int nth;
    size_t flag;
    gzFile r1, r2, r3, r4;
    kseq_t * ks1, *ks2, *ks3, *ks4;
    arrstr_k * Read1List, * Read2List, *Read3List, *Read4List;
}fqhd_k;


typedef struct fastqseq
{
    kstring_t rdn, rds, rdc, rdq;
} fqseq_k;

typedef struct fastqread
{
    fqseq_k rd1, rd2, rd3, rd4;
}fqrd_k;


typedef struct fastqpool
{
    size_t count;
    fqrd_k * fqrd;
}fqpool_k;

typedef struct fseqstr
{
    char rdn[600]; 
    char rds[600];
    char rdc[600];
    char rdq[600];
    size_t nl;
    size_t sl;
    size_t cl;
    size_t ql;
}fqseqstr_k;

typedef struct parsedfastqread
{
    fqseqstr_k rd1, rd2, rd3, rd4;
    int fb;
    char type;
}pfqrd_k;

#endif
