#ifndef _PARSE_CONFIG_H_
#define _PARSE_CONFIG_H_

#include "config_hash.h"
#include "kutils.h"

typedef struct ReadPos
{
    int rd, pos;
} rdp_k;

// ord assign a number for the barcode used in identifier. if DPM = 1; RPM = 2; then DPM-RPM-DPM == 1-2-1; for GENOME ord assigned with -1. rec: name of the barcode array
typedef struct Barcode
{
    size_t sp, la, minlen, maxlen, mist, ord;
    char * na;
    kh_m8_t * hash; 
    arrstr_k * rec;
    int bighash, random;
} bar_k;
// most number of barcode;

typedef struct BarcodeHole
{
    size_t l;
    bar_k * bar;
} barh_k;

typedef struct BarcodeChain
{
    rdp_k rdp;
    size_t l, maxb;
    int rw;
    barh_k * barh;
    char gtype;
} barch_k;

typedef struct identifier
{
    size_t l, q;
    int anchor;
    arwl_k * v;
    char * na;
}id_k;

// l is the number of barcode chain; idl is the number of idf, 
typedef struct barcodechainlist
{
    size_t sonl, l, idl;
    barch_k * barch;
    arwl_k * jump;
    id_k * idf;
    arrstr_k * sonch;
}barcl_k;

barcl_k * parseConfig(char * fname);
void barclFree(barcl_k * s);
#endif
