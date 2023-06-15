#ifndef _PARSE_SAMFIELD_H_
#define _PARSE_SAMFIELD_H_

#include "kstring.h"
#include "kseq.h"
#include "unistd.h"

typedef struct outputstruct
{
    uint64_t start;
    uint64_t end;
    int32_t tid;
}opts_k;

typedef struct optsar
{
    size_t l;
    opts_k * opts;
}optar_k;


int addOptToOpts(optar_k * opts, opts_k opt);

#endif