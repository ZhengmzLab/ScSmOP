#ifndef _CONFIG_HASH_H_
#define _CONFIG_HASH_H_

#include "khash.h"
#include "kutils.h"


// q and number pair q: 1, 2, 3;
typedef struct qnpa_k
{
    arwl_k qn;
    size_t q;
} qnpa_k;



typedef struct cawl_K
{
    size_t l;
    qnpa_k * mtar;
}cawl_k;

typedef struct harv_k
{
    size_t l;
    qnpa_k * qnpa_k;
}harv_k;



// count : current identifier have how many reads under the iden, idfcount : id of the identifier, mark : if the identifier have more than one kind of cluster
typedef struct idfv_k
{
    size_t count, idfcount;
    int mark;
    // arwl_k * svar;
} hidv_k;

typedef struct fhav_k
{
    arwl_k * rec, * q, * mist;
    size_t count;
}fhav_k;

KHASH_MAP_INIT_INT64(m64, size_t)
KHASH_MAP_INIT_INT64(m32, harv_k *)
KHASH_MAP_INIT_INT64(m16, hidv_k *)
KHASH_MAP_INIT_INT64(m8, fhav_k *)

//hash table array
typedef struct harh_k
{
    size_t l;
    kh_m64_t * ecar;
    kh_m32_t ** mtar;
}harh_k;

typedef struct BarAndRec
{
    size_t rec;
    uint64_t bar;
} bandrec_k;

void combinationFromAToB(size_t Cur, size_t n, size_t k, sstack_k * Stack, arwl_k * Result);
arwl_k * cNMinsArray(size_t s, size_t mist);
pkarwl_k * encodeSeqForHash(char * tag, size_t mist);
// uint64_t encodeSeqForReadEct(char * tag);
pkarwl_k encodeSeqForReadMis(char * tag, size_t mist, uint64_t OriValue);
void initLocker();
void destroyLocker();
harh_k * initParseHashTable(const size_t mist);
void hashDestroy(harh_k * har);
size_t querySeqBest(char * seq, harh_k * har, size_t mist);
int addSeqToHashTableXC(char * tag, harh_k * har, size_t n, size_t mist);
size_t queryIdentifier(idvp_k v, kh_m16_t * har, size_t * count, size_t idflen);
harh_k * hashCopy(harh_k * har);
void traverseHash(harh_k * har);
void initMistComArray();
kh_m8_t * initParseHashTableB();
void destroyMistComArray();
void hashDestroyB(kh_m8_t * har);
bandrec_k querySeqBestB(char * seq, char * qual, kh_m8_t * har, size_t mist, int mim_qual, int dropN);
bandrec_k querySeqBestD(char * seq, char * qual, kh_m8_t * har, size_t mist, int mim_qual, int dropN, arrstr_k * NameArr);
kh_m8_t * hashCopyB(kh_m8_t * har);
arwl_k * addSeqToHashTableXCB(char * tag, kh_m8_t * har, size_t n, size_t mist);
bandrec_k querySeqBestC(char * seq, char * qual, kh_m8_t * har, size_t mist, int min_qual, int dropN, int BarMode);
void initPrekQArray();
void destroyPrekQArray();
int writeHash(kh_m16_t * har, char * fout);
int pushIdentifierFromFile(kh_m16_t * har, char * fi);
int writeBarcodeHash(kh_m8_t * har, char * fout);
int pushBarcodeFromFile(kh_m8_t * har, char * fi);
int addCountForBarcode(kh_m8_t * har, uint64_t bar);
#endif
