#include <stdlib.h>
#include <stdio.h>
#include "kstring.h"
#include "khash.h"
#include "kutils.h"
#include "config_hash.h"
#include <pthread.h>
#include <math.h>

#define BASE_A 0b001
#define BASE_T 0b010
#define BASE_G 0b011
#define BASE_C 0b100
#define BASE_N 0b101
#define BASE_KIND 5
#define mUnit 3

// #define BASE_A 0b01
// #define BASE_T 0b10
// #define BASE_G 0b11
// #define BASE_C 0b00
// #define mUnit 2
#define MAXREADLENGTH 22
#define MAXMISMATCH 5

pthread_mutex_t locker;

arwl_k * MistCombinationArray[MAXREADLENGTH][MAXMISMATCH];
arwl_k * MistFlagArray[MAXREADLENGTH];
arwl_k * MistQarArray[MAXREADLENGTH][MAXMISMATCH];
prek_k * PrekQArray[MAXREADLENGTH][MAXMISMATCH];


void qnpaFree(qnpa_k * s)
{
    if(s)
    {
        if(s->qn.ar)
        {
            free(s->qn.ar);
            s->qn.ar = NULL;
        }
    }
}

void harvFree(harv_k * s)
{
    if(s)
    {
        if(s->qnpa_k)
        {
            qnpa_k * curqnpa = s->qnpa_k;
            for(int i = 0; i < s->l; i++)
            {
                qnpaFree(curqnpa);
                // free(curqnpa);
                curqnpa++;
            }
            free(s->qnpa_k);
        }
    }
}

// void hidvFree(hidv_k * s)
// {
//     if(s)
//     {
//         kaFree(s->svar);
//     }
// }

uint64_t encodeBase(char c)
{
    switch(c)
    {
        case 'A':
        case 'a':
            return BASE_A;
            break;
        case 'T':
        case 't':
            return BASE_T;
            break;
        case 'G':
        case 'g':
            return BASE_G;
            break;
        case 'C':
        case 'c':
            return BASE_C;
            break;
        case 'N':
        case 'n':
            return BASE_N;
            break;
        default:
            error("Non-DNA code \"%c\" found in barcode!", c);
            exit(1);
    }
}

uint64_t encodeBaseDropN(char c)
{
    switch(c)
    {
        case 'A':
        case 'a':
            return BASE_A;
            break;
        case 'T':
        case 't':
            return BASE_T;
            break;
        case 'G':
        case 'g':
            return BASE_G;
            break;
        case 'C':
        case 'c':
            return BASE_C;
            break;
        case 'N':
        case 'n':
            return -1;
            break;
        default:
            warnings("Non-DNA code \"%c\" found in barcode!", c);
            return -1;
    }
    return -1;
}

char * Decode(uint64_t st, size_t l)
{
    char * c = (char *) calloc(l+1, sizeof(char));
    uint8_t x;
    for (size_t i = l; i > 0; i--)
    {
        x = st & ((1 << mUnit) - 1);
        switch (x)
        {
            case BASE_A:
                c[i-1] = 'A';
                goto MoveRight;
            case BASE_C:
                c[i-1] = 'C';
                goto MoveRight;
            case BASE_G:
                c[i-1] = 'G';
                goto MoveRight;
            case BASE_T:
                c[i-1] = 'T';
                goto MoveRight;
            case BASE_N:
                c[i-1] = 'N';
                goto MoveRight;
            default :
                c[i-1] = 'B';
                printf("No DNA Base found\n");
                goto MoveRight;
        } 
        MoveRight:
            st = st >> mUnit;
    }
    c[l] = '\0';
    return c;
}
// return a array of [1, n, n(n-1), n(n-1)(n-mist)]
arwl_k * nMinsArray(size_t s, size_t mist)
{
    size_t ArrayLen = mist + 1, t = 1;
    arwl_k * ret = (arwl_k *) calloc(1, sizeof(arwl_k));
    ret->l = ArrayLen;
    ret->ar = (size_t *) calloc(ArrayLen, sizeof(size_t));
    ret->ar[0] = 1;
    
    for(int i = 1; i < ArrayLen; i++)
    {
        t = t * s;
        ret->ar[i] = t;
        s--;
    }
    return ret;
}

arwl_k * cNMinsArray(size_t s, size_t mist)
{
    size_t ArrayLen = mist + 1;
    arwl_k * ret = (arwl_k *) calloc(1, sizeof(arwl_k));
    ret->l = ArrayLen;
    ret->ar = (size_t *) calloc(ArrayLen, sizeof(size_t));
    ret->ar[0] = 1;

    for(int i = 1; i < ArrayLen; i++)
    {
        ret->ar[i] = cmnCount(s, i);
    }
    return ret;
}

// Get Cnm combination
void combinationFromAToB(size_t Cur, size_t n, size_t k, sstack_k * Stack, arwl_k * Result)
{
    size_t TopElement;
    if(Stack->length == k)
    {
        for(int i = 0; i < Stack->length; i++)
        {
            // printf("Stack.element[%d]: %ld\n", i, Stack->element[i]);
            addEleToArwlRep(Result, Stack->element[i]);
        }
        return;
    }
    for(; Cur <= n;)
    {
        // printf("i : %ld\n", i);
        pushNumStack(Stack, Cur);
        combinationFromAToB(++Cur, n, k, Stack, Result);
        popNumStack(Stack, &TopElement);
    }
    return;
}

void initMistComArray()
{
    sstack_k * Stack = (sstack_k *) calloc(1, sizeof(sstack_k));
    for(int i = 0; i < MAXREADLENGTH; i++)
    {
        MistFlagArray[i] = (arwl_k *) calloc(1, sizeof(arwl_k));
        MistFlagArray[i]->ar = (size_t *) calloc(i, sizeof(size_t));
        MistFlagArray[i]->l = i;
        for(int k = 0; k < i; k++) MistFlagArray[i]->ar[k] = 1 << k;

        for(int j = 0; j < MAXMISMATCH; j++)
        {
            initNumStack(Stack);
            MistQarArray[i][j] = cNMinsArray(i, j);
            MistCombinationArray[i][j] = (arwl_k *) calloc(1, sizeof(arwl_k));
            combinationFromAToB(1, i, j, Stack, MistCombinationArray[i][j]);
        }
    }
    free(Stack);
}

// q array 111111 and pkarraylen stored at TmpAr[0].v oriq stored at TmpAr[1].v
void initPrekQArray()
{
    for(int b = 0; b < MAXREADLENGTH; b++)
    {
        size_t QLen = b;
        for(size_t mist = 0; mist < MAXMISMATCH; mist++)
        {
            size_t PKArrayLen = 0;
            arwl_k * Qar = MistQarArray[QLen][mist];
            size_t Power = 1, TmpCount = 0;
            for(int i = 0; i < Qar->l; i++)
            {
                TmpCount = TmpCount + Qar->ar[i];
                PKArrayLen = PKArrayLen + Qar->ar[i] * Power;
                Power = Power * BASE_KIND;
            }
            Power = (size_t) Power/BASE_KIND;
            PrekQArray[QLen][mist] = (prek_k *) calloc(TmpCount, sizeof(prek_k));
            prek_k * TmpAr = PrekQArray[QLen][mist];
            prek_k * curpk = TmpAr;
            arwl_k * FlagArray = MistFlagArray[QLen];
            size_t Oriq = (1 << QLen) - 1;
            for(int i = 1; i <= mist; i++)
            {
                arwl_k * Result = MistCombinationArray[QLen][i];
                for(int j = 0; j < (Result->l/i); j++)
                {
                    size_t Tmpq = Oriq;
                    for(int k = 0; k < i; k++)
                    {
                        Tmpq = Tmpq ^ FlagArray->ar[(Result->ar[(j * i) + k] - 1)];
                    }
                    (++curpk)->q = Tmpq;
                }
            }
            TmpAr[0].v = PKArrayLen;
            TmpAr[1].v = Oriq;
        }
    }
}

void destroyPrekQArray()
{
    for(int i = 0; i < MAXREADLENGTH; i++)
    {
        for(size_t mist = 0; mist < MAXMISMATCH; mist++)
        {
            free(PrekQArray[i][mist]);
        }
    }
}

void initLocker()
{
    if(pthread_mutex_init(&locker, NULL) != 0)
    {
        printf("Locker initialize failed!\n");
    }
}

void destroyLocker()
{
    pthread_mutex_destroy(&locker);
}

pkarwl_k * encodeSeqForHash(char * tag, size_t mist)
{
    if(mist < 0) mist = 0;
    size_t QLen = strlen(tag), PKArrayLen = 0;
    arwl_k * Qar = MistQarArray[QLen][mist];
    // pk array length for the seq
    for( int i = 0; i < Qar->l; i++) PKArrayLen = PKArrayLen + Qar->ar[i];
    // pk array construct
    pkarwl_k * ret = (pkarwl_k *) calloc(1, sizeof(pkarwl_k));
    ret->ar = (prek_k *) calloc(PKArrayLen, sizeof(prek_k));
    ret->l = PKArrayLen;
    ret->qar = Qar;
    ret->qlar = (size_t *) calloc(PKArrayLen, sizeof(size_t));

    uint64_t OriValue = 0;
    for(int i = 0; i < QLen; i++) OriValue = (OriValue << mUnit | encodeBase(tag[i]));
    ret->ar[0].v = OriValue;
    size_t Oriq = (1 << QLen) - 1;
    ret->ar[0].q = Oriq;
    ret->qlar[0] = QLen;
    // mismatch array
    size_t * m = ret->qlar;
    prek_k * curpk = ret->ar;
    arwl_k * FlagArray = MistFlagArray[QLen];
    for(int i = 1; i <= mist; i++)
    {
        arwl_k * Result = MistCombinationArray[QLen][i];
        for(int j = 0; j < (Result->l/i); j++)
        {
            size_t Tmpq = Oriq;
            for(int k = 0; k < i; k++)
            {
                Tmpq = Tmpq ^ FlagArray->ar[(Result->ar[(j * i) + k] - 1)];
            }
            *(++m) = QLen - i;
            (++curpk)->q = Tmpq;
        }
    }
    
    curpk = ret->ar;
    for(int i = 1; i < ret->l; i++)
    {
        size_t MisFlag = 1;
        uint64_t ValueMisFlag = 0b111;
        uint64_t TmpValue = 0;
        for(int j = 0; j < QLen; j++)
        {
            TmpValue = ret->ar[i].q & MisFlag ? (TmpValue | (OriValue & ValueMisFlag)) : TmpValue;
            MisFlag = MisFlag << 1;
            ValueMisFlag = ValueMisFlag << mUnit;
        }
        (++curpk)->v = TmpValue;
    }
    return ret;
}

pkarwl_k * encodeSeqForHashB(char * tag, size_t mist)
{
    if(mist < 0) mist = 0;
    size_t QLen = strlen(tag), PKArrayLen = 0;
    arwl_k * Qar = MistQarArray[QLen][mist];
    // pk array length for the seq
    size_t Power = 1, TmpCount = 0;
    for( int i = 0; i < Qar->l; i++)
    {
        TmpCount = TmpCount + Qar->ar[i];
        PKArrayLen = PKArrayLen + Qar->ar[i] * Power;
        Power = Power * BASE_KIND;
    }
    Power = (size_t) Power/BASE_KIND;
    // pk array construct
    pkarwl_k * ret = (pkarwl_k *) calloc(1, sizeof(pkarwl_k));
    ret->ar = (prek_k *) calloc(PKArrayLen, sizeof(prek_k));
    ret->l = PKArrayLen;
    ret->qar = Qar;
    ret->qlar = (size_t *) calloc(PKArrayLen, sizeof(size_t));
    uint64_t OriValue = 0;
    for(int i = 0; i < QLen; i++) OriValue = (OriValue << mUnit | encodeBase(tag[i]));
    ret->ar[0].v = OriValue;
    size_t Oriq = (1 << QLen) - 1;
    ret->ar[0].q = Oriq;
    ret->qlar[0] = QLen;
    // mismatch array
    size_t * TmpQlar = (size_t *) calloc(TmpCount, sizeof(size_t));
    size_t * m = TmpQlar;
    prek_k * TmpAr = (prek_k *) calloc(TmpCount, sizeof(prek_k));
    prek_k * curpk = TmpAr;
    arwl_k * FlagArray = MistFlagArray[QLen];
    for(int i = 1; i <= mist; i++)
    {
        arwl_k * Result = MistCombinationArray[QLen][i];
        for(int j = 0; j < (Result->l/i); j++)
        {
            size_t Tmpq = Oriq;
            for(int k = 0; k < i; k++)
            {
                Tmpq = Tmpq ^ FlagArray->ar[(Result->ar[(j * i) + k] - 1)];
            }
            *(++m) = QLen - i;
            (++curpk)->q = Tmpq;
        }
    }
    
    curpk = ret->ar;
    int CI = 1;
    size_t TmpPower = 1;
    for(int i = 1; i <= mist; i++)
    {
        TmpPower = TmpPower * BASE_KIND;
        for(int m = 0; m < Qar->ar[i]; m++)
        {
            size_t MisFlag = 1;
            uint64_t ValueMisFlag = 0b111;
            uint64_t BaseArray[BASE_KIND] = {0b001, 0b010, 0b011, 0b100, 0b000};
            uint64_t * ValueArray = (uint64_t *) calloc(TmpPower, sizeof(uint64_t));
            size_t InnerPower = BASE_KIND;

            for(int j = 0; j < QLen; j++)
            {
                if(TmpAr[CI].q & MisFlag)
                {
                    for(int k = 0; k < TmpPower; k++)
                    {
                        ValueArray[k] = ValueArray[k] | (OriValue & ValueMisFlag);
                    }
                }
                else
                {
                    size_t RepTime = (size_t) TmpPower/InnerPower;
                    size_t k = 0, ICU = 0;
                    while(k < TmpPower)
                    {
                        for(int TmpCounter = 0; TmpCounter < RepTime; TmpCounter++)
                        {
                            ValueArray[k] = ValueArray[k] | BaseArray[ICU];
                            k++;
                        }
                        ICU++;
                        if(ICU%BASE_KIND == 0) ICU = 0;
                    }
                    InnerPower = InnerPower * BASE_KIND;
                }
                MisFlag = MisFlag << 1;
                ValueMisFlag = ValueMisFlag << mUnit;
                for(int j = 0; j < 4; j++) BaseArray[j] = BaseArray[j] << mUnit;
            }
            for(int j = 0; j < TmpPower; j++)
            {
                ++curpk;
                curpk->v = ValueArray[j];
                curpk->q = TmpAr[CI].q;
            }
            free(ValueArray);
            CI++;
        }
    }
    free(TmpAr);
    free(TmpQlar);
    return ret;
}

static inline uint64_t encodeSeqForReadEct(char * tag, int dropN)
{
    size_t QLen = strlen(tag);
    uint64_t OriValue = 0;
    if(dropN != 1)
    {
        for(int i = 0; i < QLen; i++) 
        {
            OriValue = (OriValue << mUnit | encodeBase(tag[i]));
        }
    }
    else
    {
        for(int i = 0; i < QLen; i++)
        {
            if(encodeBaseDropN(tag[i]) != -1)
            {
                OriValue = (OriValue << mUnit | encodeBase(tag[i]));
            }
            else
            {
                OriValue = 0;
                return OriValue;
            }
        }
    }
    return OriValue;
}

size_t queryMisSeq(char * tag, size_t mist, uint64_t OriValue, harh_k * har)
{
    size_t TmpRet = SIZE_MAX;
    size_t QLen = strlen(tag), PKArrayLen = 0;
    arwl_k * Qar = MistQarArray[QLen][mist];
    PKArrayLen = Qar->ar[mist];
    arwl_k * FlagArray = MistFlagArray[QLen];
    arwl_k * Result = MistCombinationArray[QLen][mist];
    for(int i = 0; i < PKArrayLen; i++)
    {
        size_t Tmpq = (1 << QLen) - 1;
        for(int k = 0; k < mist; k++) Tmpq = Tmpq ^ FlagArray->ar[(Result->ar[(i * mist) + k] - 1)];
        size_t MisFlag = 1;
        uint64_t ValueMisFlag = 0b111;
        uint64_t TmpValue = 0;
        for(int j = 0; j < QLen; j++)
        {
            TmpValue = Tmpq & MisFlag ? (TmpValue | (OriValue & ValueMisFlag)) : TmpValue;
            MisFlag = MisFlag << 1;
            ValueMisFlag = ValueMisFlag << mUnit;
        }

        khint_t k = kh_get(m32, har->mtar[mist - 1], TmpValue);
        if(k != kh_end(har->mtar[mist - 1]))
        {
            harv_k * curhv;
            curhv = kh_value(har->mtar[mist - 1], k);
            for(int o = 0; o < curhv->l; o++)
            {
                if(Tmpq == curhv->qnpa_k[o].q)
                {
                    // if there are conflict in barcode
                    // if(curhv->qnpa_k->qn.l > 1)
                    // {
                    //     TmpRet = curhv->qnpa_k->qn.ar[0];
                    //     return TmpRet;
                    // }
                    // else
                    // {
                    //     TmpRet = curhv->qnpa_k->qn.ar[0];
                    //     return TmpRet;
                    // }
                    TmpRet = curhv->qnpa_k->qn.ar[0];
                    return TmpRet;
                }
            }
        }
    }

    return TmpRet;
}

harh_k * initParseHashTable(const size_t mist)
{
    harh_k * har = (harh_k *) calloc(1, sizeof(harh_k));
    har->l = mist + 1;
    har->ecar = kh_init(m64);
    har->mtar = (kh_m32_t **) calloc(mist, sizeof(kh_m32_t *));

    for(int i = 0; i < mist; i++)
    {
        har->mtar[i] = kh_init(m32);
    }
    return har;
}

kh_m8_t * initParseHashTableB()
{
    kh_m8_t * har = kh_init(m8);
    return har;
}

void hashDestroy(harh_k * har)
{
    if(har)
    {
        kh_destroy_m64(har->ecar);
        // int count = 0;
        for(int i = 0; i < har->l - 1; i++)
        {
            // printf("I: %d\n", i);
            if(har->mtar[i])
            {
                khint_t k;
                for( k = 0; k < kh_end(har->mtar[i]); ++k)
                {
                    if(kh_exist(har->mtar[i], k))
                    {
                        harvFree(kh_value(har->mtar[i], k));
                        free(kh_value(har->mtar[i], k));
                    }
                }
            }
            kh_destroy_m32(har->mtar[i]);
        }
        free(har->mtar);
        har->mtar = NULL;
    }
}

void hashDestroyB(kh_m8_t * har)
{
    if(har)
    {
        khint_t k;
        for(k = 0; k < kh_end(har); ++k)
        {
            if(kh_exist(har, k))
            {
                fhav_k * hav = kh_value(har, k);
                kaFree(hav->mist);
                free(hav->mist);
                kaFree(hav->q);
                free(hav->q);
                kaFree(hav->rec);
                free(hav->rec);
            }
        }
    }
    kh_destroy_m8(har);
}

void traverseHash(harh_k * har)
{
    if(har)
    {
        printf("traversing hash table....\n");
        printf("Hash table size: %d\n", kh_size(har->ecar));
        khint_t k;
        for( k = 0; k < kh_end(har->ecar); ++k)
        {
            if(kh_exist(har->ecar, k))
            {
                printf("value: %ld\n", kh_value(har->ecar, k));
            }
        }

        for(int i = 0; i < har->l - 1; i++)
        {
            for(k = 0; k < kh_end(har->mtar[i]); ++k)
            {
                if(kh_exist(har->mtar[i], k))
                {
                    harv_k * OriHarv = kh_value(har->mtar[i], k);
                    printf("harv length : %ld\n", OriHarv->l);
                    printf("harv qnpa.qn.l %ld\n", OriHarv->qnpa_k->qn.l);
                    printf("harv q: %ld\n", OriHarv->qnpa_k->q);
                }
            }
        }
    }
}

int writeBarcodeHash(kh_m8_t * har, char * fout)
{
    FILE * fo = fopen(fout, "w");
    fhav_k * ret;
    khint_t k;
    uint64_t key;
    for( k = 0; k < kh_end(har); ++k)
    {
        if(kh_exist(har, k))
        {
            key = kh_key(har, k);
            ret = kh_val(har, k);
            // key barcode q count
            fprintf(fo, "%lu\t%ld\t%ld\t%ld\t%s\n", key, ret->rec->ar[0], ret->q->ar[0], ret->count, Decode(key, 16));  
        }
    }
    fclose(fo);
    return 0;
}

harh_k * hashCopy(harh_k * har)
{
    harh_k * target;
    if(har)
    {
        target = initParseHashTable(har->l - 1);
        khint_t k, q;
        int absent;
        for( k = 0; k < kh_end(har->ecar); ++k)
        {
            if(kh_exist(har->ecar, k))
            {
                // printf("%d\n", count++);
                q = kh_put(m64, target->ecar, kh_key(har->ecar, k), &absent);
                kh_value(target->ecar, q) = kh_value(har->ecar, k);
            }
        }
        for(int i = 0; i < har->l - 1; i++)
        {
            for(k = 0; k < kh_end(har->mtar[i]); ++k)
            {
                if(kh_exist(har->mtar[i], k))
                {
                    q = kh_put(m32, target->mtar[i], kh_key(har->mtar[i], k), &absent);
                    harv_k * OriHarv = kh_value(har->mtar[i], k);
                    harv_k * NewHarv = (harv_k *) calloc(1, sizeof(harv_k));
                    NewHarv->qnpa_k = (qnpa_k *) calloc(1, sizeof(qnpa_k));
                    NewHarv->l = OriHarv->l;
                    NewHarv->qnpa_k->qn.l = OriHarv->qnpa_k->qn.l;
                    NewHarv->qnpa_k->qn.ar = (size_t *) calloc(NewHarv->qnpa_k->qn.l, sizeof(size_t));
                    memmove(NewHarv->qnpa_k->qn.ar, OriHarv->qnpa_k->qn.ar, NewHarv->qnpa_k->qn.l * sizeof(size_t));
                    NewHarv->qnpa_k->q = OriHarv->qnpa_k->q;
                    kh_value(target->mtar[i], q) = NewHarv;
                }
            }
        }
    }
    else
    {
        target = NULL;
    }
    return target;
}

kh_m8_t * hashCopyB(kh_m8_t * har)
{
    kh_m8_t * target;
    if(har)
    {
        target = kh_init(m8);
        khint_t k, q;
        int absent;
        for(k = 0; k < kh_end(har); ++k)
        {
            if(kh_exist(har, k))
            {
                q = kh_put(m8, target, kh_key(har, k), &absent);
                fhav_k * Orihav = kh_value(har, k);
                fhav_k * Newhav = (fhav_k *) calloc(1, sizeof(Newhav));
                Newhav->q = (arwl_k *) calloc(1, sizeof(arwl_k));
                Newhav->mist = (arwl_k *) calloc(1, sizeof(arwl_k));
                Newhav->rec = (arwl_k *) calloc(1, sizeof(arwl_k));
                for(int i = 0; i < Orihav->q->l; i++) 
                {
                    addEleToArwlRep(Newhav->q, Orihav->q->ar[i]);
                    addEleToArwlRep(Newhav->mist, Orihav->mist->ar[i]);
                    addEleToArwlRep(Newhav->rec, Orihav->rec->ar[i]);
                }
                kh_value(target, q) = Newhav;
            }
        }
    }
    else
    {
        target = NULL;
    }
    return target;
}

int addQnpaPairToHarv(harv_k * s, qnpa_k * m)
{
    qnpa_k * Tmp = (qnpa_k *) calloc(s->l + 1, sizeof(qnpa_k));
    memmove(Tmp, s->qnpa_k, s->l * sizeof(qnpa_k));
    memmove(Tmp + s->l, m, sizeof(qnpa_k));
    free(s->qnpa_k);
    free(m);
    s->qnpa_k = Tmp;
    s->l = s->l + 1;
    return 0;
}

int addSeqToHashTableXC(char * tag, harh_k * har, size_t n, size_t mist)
{
    pkarwl_k * k_h;
    k_h = encodeSeqForHash(tag, mist);
    int absent;
    prek_k * curpk = k_h->ar;
    khint_t k;
    k = kh_get(m64, har->ecar, curpk->v);

    if(k == kh_end(har->ecar))
    {
        k = kh_put(m64, har->ecar, curpk->v, &absent);
        kh_value(har->ecar, k) = n;
    }
    else
    {
        error("%s have existed! You may want to check the barcode list!\n", Decode(curpk->v, k_h->qlar[0]));
        kpFree(k_h);
        free(k_h);
        return 0;
    }
    for(int i = 0; i < mist; i++)
    {
        for(int j = 0; j < k_h->qar->ar[i+1]; j++)
        {
            curpk++;
            qnpa_k * curqk = (qnpa_k *) calloc(1, sizeof(qnpa_k));
            curqk->q = curpk->q;
            curqk->qn.l = 1;

            k = kh_get(m32, har->mtar[i], curpk->v);
            if(k == kh_end(har->mtar[i]))
            {
                harv_k * newhv = (harv_k *) calloc(1, sizeof(harv_k));
                newhv->qnpa_k = (qnpa_k *) calloc(1, sizeof(qnpa_k));
                newhv->l = 1;
                newhv->qnpa_k->qn.ar = (size_t *) calloc(1, sizeof(size_t));
                newhv->qnpa_k->qn.ar[0]  = n;
                newhv->qnpa_k->q = curpk->q;
                newhv->qnpa_k->qn.l = 1;
                k = kh_put(m32, har->mtar[i], curpk->v, &absent);
                kh_value(har->mtar[i], k) = newhv;
            }
            else
            {
                harv_k * curhv;
                curhv = kh_value(har->mtar[i], k);
                int o = 0;
                for(int o = 0; o < curhv->l; o++)
                {
                    if(curpk->q == curhv->qnpa_k[o].q)
                    {
                        addEleToArwlXC(&curhv->qnpa_k[o].qn, n);
                        break;
                    }
                }
                if(o == curhv->l)
                {
                    curqk->qn.ar = (size_t *) calloc(1, sizeof(size_t));
                    curqk->qn.ar[0] = n;
                    addQnpaPairToHarv(curhv, curqk);
                }
                kh_value(har->mtar[i], k) = curhv;
            }
            free(curqk);
        }
    }
    kpFree(k_h);
    free(k_h);
    return 0;
}

// this version is to add all mismatched sequence to hash table
arwl_k * addSeqToHashTableXCB(char * tag, kh_m8_t * har, size_t n, size_t mist)
{
    pkarwl_k * k_h;
    arwl_k * ret = (arwl_k *) calloc(1, sizeof(arwl_k));
    k_h = encodeSeqForHashB(tag, mist);
    int absent;
    prek_k * curpk = k_h->ar;
    // printf("Value : %lu\n", curpk->v);
    khint_t k;
    size_t tn = strlen(tag);
    size_t cmist;
    for(int i = 0; i < k_h->l; i++)
    {
        cmist = tn - k_h->qlar[i];
        k = kh_get(m8, har, curpk->v);
        if(k == kh_end(har))
        {
            k = kh_put(m8, har, curpk->v, &absent);
            fhav_k * newhv = (fhav_k *) calloc(1, sizeof(fhav_k));
            newhv->mist = (arwl_k *) calloc(1, sizeof(arwl_k));
            newhv->rec = (arwl_k *) calloc(1, sizeof(arwl_k));
            newhv->q = (arwl_k *) calloc(1, sizeof(arwl_k));
            addEleToArwlXC(newhv->mist, cmist);
            addEleToArwlXC(newhv->q, curpk->q);
            addEleToArwlXC(newhv->rec, n);
            kh_value(har, k) = newhv;
            addEleToArwlXC(ret, n);
        }
        else
        {
            fhav_k * curhv;
            curhv = kh_value(har, k);
            for(int j = 0; j < curhv->rec->l; j++)
            {
                addEleToArwlXC(ret, curhv->rec->ar[j]);
            }
            int Top = addEleToArwlXC(ret, n);
            if(Top != -1)
            {
                int Pos = addEleToArwlXCSorted(curhv->mist, mist);
                insertEleToArwl(curhv->q, curpk->q, Pos);
                insertEleToArwl(curhv->rec, n, Pos);
            }
            kh_value(har, k) = curhv;
        }
        curpk++;
    }
    kpFree(k_h);
    free(k_h);
    return ret;
}
// random sequence as barcode, no mismatch allowed.
arwl_k * addSeqToHashTableXCD(char * tag, kh_m8_t * har, size_t n, size_t mist, arrstr_k * NameArr)
{
    pthread_mutex_lock(&locker);
    addStrToArrstrRepXC(NameArr, tag);
    pkarwl_k * k_h;
    arwl_k * ret = (arwl_k *) calloc(1, sizeof(arwl_k));
    k_h = encodeSeqForHashB(tag, mist);
    int absent;
    prek_k * curpk = k_h->ar;
    khint_t k;
    size_t tn = strlen(tag);
    size_t cmist;
    for(int i = 0; i < k_h->l; i++)
    {
        cmist = tn - k_h->qlar[i];
        k = kh_get(m8, har, curpk->v);
        if(k == kh_end(har))
        {
            k = kh_put(m8, har, curpk->v, &absent);
            fhav_k * newhv = (fhav_k *) calloc(1, sizeof(fhav_k));
            newhv->mist = (arwl_k *) calloc(1, sizeof(arwl_k));
            newhv->rec = (arwl_k *) calloc(1, sizeof(arwl_k));
            newhv->q = (arwl_k *) calloc(1, sizeof(arwl_k));
            addEleToArwlXC(newhv->mist, cmist);
            addEleToArwlXC(newhv->q, curpk->q);
            addEleToArwlXC(newhv->rec, n);
            kh_value(har, k) = newhv;
            addEleToArwlXC(ret, n);
        }
        curpk++;
    }
    kpFree(k_h);
    free(k_h);
    pthread_mutex_unlock(&locker);
    return ret;
}
size_t querySeqBest(char * seq, harh_k * har, size_t mist)
{

    uint64_t EctValue = encodeSeqForReadEct(seq, 0);
    size_t ret = SIZE_MAX;

    khint_t k = kh_get(m64, har->ecar, EctValue);
    if(k != kh_end(har->ecar))
    {
        ret = kh_value(har->ecar, k);
        return ret;
    }
    
    for(int i = 0; i < mist; i++)
    {
        size_t TmpRet = queryMisSeq(seq, i+1, EctValue, har);
        if(TmpRet != SIZE_MAX) return TmpRet;
    }

    return ret;
}

size_t queryMisSeqB(char * tag, size_t mist, uint64_t OriValue, kh_m8_t * har)
{
    size_t TmpRet = SIZE_MAX;
    size_t QLen = strlen(tag), PKArrayLen = 0;
    arwl_k * Qar = MistQarArray[QLen][mist];
    PKArrayLen = Qar->ar[mist];
    arwl_k * FlagArray = MistFlagArray[QLen];
    arwl_k * Result = MistCombinationArray[QLen][mist];

    for(int i = 0; i < PKArrayLen; i++)
    {
        size_t Tmpq = (1 << QLen) - 1;
        for(int k = 0; k < mist; k++) Tmpq = Tmpq ^ FlagArray->ar[(Result->ar[(i * mist) + k] - 1)];
        size_t MisFlag = 1;
        uint64_t ValueMisFlag = 0b111;
        uint64_t TmpValue = 0;
        for(int j = 0; j < QLen; j++)
        {
            TmpValue = Tmpq & MisFlag ? (TmpValue | (OriValue & ValueMisFlag)) : TmpValue;
            MisFlag = MisFlag << 1;
            ValueMisFlag = ValueMisFlag << mUnit;
        }
        khint_t k = kh_get(m8, har, TmpValue);
        if(k != kh_end(har))
        {
            fhav_k * hav = kh_value(har, k);
            TmpRet = hav->rec->ar[0];
            hav->count++;
            // if(hav->count > 0) printf("%ld\n", hav->count);
            return TmpRet;
        }
    }
    return TmpRet;
}

bandrec_k querySeqBestB(char * seq, char * qual, kh_m8_t * har, size_t mist, int min_qual, int dropN)
{
    bandrec_k lxd;
    lxd.rec = SIZE_MAX;
    if(min_qual > 0)
    {
        size_t SumQual = 0;
        for(int i = 0; i < strlen(qual); i++)
        {
            SumQual = SumQual + qual[i] - 33;
        }
        // printf("SumQual : %ld - MinQual : %ld\n", SumQual, (min_qual * strlen(qual)));
        if(SumQual < (min_qual * strlen(qual))) return lxd;
    }
    uint64_t EctValue = encodeSeqForReadEct(seq, dropN);
    if(EctValue > 0) 
    {
        fhav_k * hav;
        khint_t k = kh_get(m8, har, EctValue);
        if(k != kh_end(har))
        {
            hav = kh_value(har, k);
            lxd.rec = hav->rec->ar[0];
            // lxd.bar = EctValue;
            return lxd;
        }
    }
    return lxd;
}

// random sequence as barcode, no mismatch allowed.
bandrec_k querySeqBestD(char * seq, char * qual, kh_m8_t * har, size_t mist, int min_qual, int dropN, arrstr_k * NameArr)
{
    bandrec_k lxd;
    lxd.rec = SIZE_MAX;
    if(min_qual > 0)
    {
        size_t SumQual = 0;
        for(int i = 0; i < strlen(qual); i++)
        {
            SumQual = SumQual + qual[i] - 33;
        }
        // printf("SumQual : %ld - MinQual : %ld\n", SumQual, (min_qual * strlen(qual)));
        if(SumQual < (min_qual * strlen(qual))) return lxd;
    }
    uint64_t EctValue = encodeSeqForReadEct(seq, dropN);
    if(EctValue > 0) 
    {
        fhav_k * hav;
        khint_t k = kh_get(m8, har, EctValue);
        if(k != kh_end(har))
        {
            hav = kh_value(har, k);
            lxd.rec = hav->rec->ar[0];
            // lxd.bar = EctValue;
            return lxd;
        }
        else
        {
            size_t count = kh_size(har);
            arwl_k * Conflict = addSeqToHashTableXCD(seq, har, count, mist, NameArr);
            lxd.rec = count;
            if(count >= SIZE_MAX) error("Barcode number exceeded limits: %zu", SIZE_MAX);
            kaFree(Conflict);
            free(Conflict);
            return lxd;
        }
    }
    return lxd;
}

int mycmp(const void * a, const void * b)
{
    return (*(size_t *) a - *(size_t *) b);
}

// decide 10x rna-seq
int DeciderTenXR(arwl_k * s)
{
    
    qsort(s->ar, s->l, sizeof(size_t), mycmp);
    float * ProVec = (float *) calloc(s->l, sizeof(float));
    float Sum = 0;
    for(int i = 0; i < s->l; i++)
    {
        ProVec[i] = pow(10, -(((float)s->ar[i])/10));
        Sum = Sum + ProVec[i];
    }
    for(int i = 0; i < s->l; i++)
    {
        ProVec[i] = ProVec[i]/Sum;
        if(ProVec[i] >= 0.975) 
        {
            free(ProVec);
            return 1;
        }
    }
    free(ProVec);
    return 0;
}

// if mismatch at pos poss less than 90 return 0 else return 1
int DeciderTenXATAC(arwl_k * s, arwl_k * c)
{
    // if(AllSumQ == 1) return 1;
    // if(s->ar[pos] < 12) return 1;
    float * ProVec = (float *) calloc(s->l, sizeof(float));
    float Sum = 0;
    for(int i = 0; i < s->l; i++)
    {
        ProVec[i] = pow(10, -(((float)s->ar[i])/10)) * (c->ar[i] > 0 ? c->ar[i] : 1);
        Sum = Sum + ProVec[i];
    }
    for(int i = 0; i < s->l; i++)
    {
        ProVec[i] = ProVec[i]/Sum;
        if(ProVec[i] >= 0.975) 
        {
            free(ProVec);
            return i;
        }
    }
    free(ProVec);
    return -1;
}

// This function is aimed to 10x with quality selection.
bandrec_k querySeqBestC(char * seq, char * qual, kh_m8_t * har, size_t mist, int min_qual, int dropN, int BarMode)
{
    bandrec_k lxd = {0};
    lxd.rec = SIZE_MAX;
    if(min_qual > 0)
    {
        size_t SumQual = 0;
        for(int i = 0; i < strlen(qual); i++)
        {
            SumQual = SumQual + qual[i] - 33;
        }
        if(SumQual < (size_t) (min_qual * strlen(qual))) return lxd;
    }
    uint64_t EctValue = encodeSeqForReadEct(seq, dropN);
    if(EctValue > 0)
    {
        fhav_k * hav;
        khint_t k = kh_get(m8, har, EctValue);
        if(k != kh_end(har))
        {
            hav = kh_value(har, k);
            lxd.rec = hav->rec->ar[0];
            lxd.bar = EctValue;
            return lxd;
        }
        else
        {
            if(mist == 0) return lxd;
            // barcode number in whitelistll
            size_t AllSumQ = 0;
            arwl_k * QualArr = (arwl_k *) calloc(1, sizeof(arwl_k));
            arwl_k * CountArr;
            arwl_k * RecArr;
            if(BarMode == 2)
            {
                CountArr = (arwl_k *) calloc(1, sizeof(arwl_k));
                RecArr = (arwl_k *) calloc(1, sizeof(arwl_k));
            }
            // arwl_k * Result = (arwl_k *) calloc(1, sizeof(arwl_k));
            size_t QLen = strlen(seq);
            uint64_t OriValue = EctValue;
            prek_k * TmpAr = PrekQArray[QLen][mist];
            size_t Oriq = TmpAr[1].v;
            arwl_k * Qar = MistQarArray[QLen][mist];
            int CI = 1;
            size_t TmpPower = 1;
            size_t MinsQual = SIZE_MAX;
            for(int i = 1; i <= mist; i++)
            {
                TmpPower = TmpPower * BASE_KIND;
                for(int m = 0; m < Qar->ar[i]; m++)
                {
                    size_t MisFlag = 1;
                    uint64_t ValueMisFlag = 0b111;
                    uint64_t BaseArray[BASE_KIND] = {BASE_A, BASE_C, BASE_G, BASE_T, BASE_N};
                    size_t InnerPower = BASE_KIND;
                    // loop each q
                    uint64_t ValueArray[128] = {0};
                    for(int j = 0; j < QLen; j++)
                    {
                        if(TmpAr[CI].q & MisFlag)
                        {
                            for(int k = 0; k < TmpPower; k++)
                            {
                                ValueArray[k] = ValueArray[k] | (OriValue & ValueMisFlag);
                            }
                        }
                        else
                        {
                            size_t RepTime = (size_t) TmpPower/InnerPower;
                            size_t k = 0, ICU = 0;
                            while(k < TmpPower)
                            {
                                for(int TmpCounter = 0; TmpCounter < RepTime; TmpCounter++)
                                {
                                    ValueArray[k] = ValueArray[k] | BaseArray[ICU%BASE_KIND];
                                    k++;
                                }
                                ICU++;
                            }
                            InnerPower = InnerPower * BASE_KIND;
                        }
                        MisFlag = MisFlag << 1;
                        ValueMisFlag = ValueMisFlag << mUnit;
                        for(int p = 0; p < BASE_KIND; p++) BaseArray[p] = BaseArray[p] << mUnit;
                    }
                    for(int p = 0; p < TmpPower; p++)
                    {
                        // summary of the the exact barcode missing quality. 1 mist match is one
                        size_t SumQual = 0;
                        size_t TmpValue = SIZE_MAX;
                        size_t MistFlag = 1;
                        k = kh_get(m8, har, ValueArray[p]);
                        if(k != kh_end(har))
                        {
                            hav = kh_value(har, k);
                            TmpValue = hav->rec->ar[0];
                            AllSumQ++;
                            size_t MistPos = Oriq ^ TmpAr[CI].q;
                            for(int j = QLen - 1; j >= 0; j--)
                            {
                                if(MistPos & MistFlag)
                                {
                                    SumQual = SumQual + qual[j] - 33;
                                }
                                MistFlag = MistFlag << 1;
                            }
                            if(SumQual < MinsQual)
                            {
                                MinsQual = SumQual;
                                lxd.rec = TmpValue;
                                lxd.bar = ValueArray[p];
                            }
                            // printf("Tmp : %lu\n", ValueArray[p]);
                            // printf("Value : %ld\n", hav->rec->ar[0]);
                            // printf("WhiteListSeq : %s\t", Decode(ValueArray[p], 16));
                            // printf("Qual: %ld\t", SumQual);
                            // printf("Ect Count : %ld\t", hav->count);
                            // ret always represent for the min quality barcode.
                            // mode 0 : default : find and stop; 1 : 10x rna seq; 2 10x atac seq; 3 minus quality
                            switch(BarMode)
                            {
                                case 0: 
                                {
                                    kaFree(QualArr);
                                    free(QualArr);
                                    return lxd;
                                }
                                case 1:
                                {
                                    addEleToArwlRep(QualArr, SumQual);
                                    break;
                                }
                                case 2:
                                {
                                    addEleToArwlRep(CountArr, hav->count);
                                    addEleToArwlRep(QualArr, SumQual);
                                    addEleToArwlRep(RecArr, TmpValue);
                                    break;
                                }
                                case 3:
                                {
                                    break;
                                }
                            }

                        }
                    }
                    CI++;
                }
                switch(BarMode)
                {
                    case 2 :
                    {
                        int Pos = DeciderTenXATAC(QualArr, CountArr);
                        if( Pos >= 0) lxd.rec = RecArr->ar[Pos]; else lxd.rec = SIZE_MAX;
                        break;
                    }
                    case 1 :
                    {
                        if(!DeciderTenXR(QualArr)) lxd.rec = SIZE_MAX;
                        break;
                    }
                }
                if(lxd.rec != SIZE_MAX) 
                {
                    kaFree(QualArr);
                    free(QualArr);
                    return lxd;
                }
            }
            if(BarMode == 2)
            {
                kaFree(CountArr);
                kaFree(RecArr);
                free(CountArr);
                free(RecArr);
            }
            kaFree(QualArr);
            free(QualArr);
            return lxd;
        }
    }
    return lxd;
}

int pushIdentifierFromFile(kh_m16_t * har, char * fi)
{
    FILE * f = fopen(fi, "r");
    if(f == NULL) error("idb file do not exist\n");
    kstring_t line = {0, 0, NULL};
    khint_t k;
    uint64_t key = 0;
    char * over;
    int absent;
    size_t idfcount = 0;
    size_t maxidfcount = 0;
    while(kgetline(&line, (kgets_func *) fgets, (FILE *) f) != EOF)
    {
        arrstr_k * StrArr = splitStrByChar(line.s, '\t');
        idfcount = strtol(StrArr->mem + StrArr->offset_ar[1], &over, 10);
        if(maxidfcount < idfcount) maxidfcount = idfcount;
        key = strtoul(StrArr->mem + StrArr->offset_ar[0], &over, 10);
        k = kh_get(m16, har, key);
        if(k == kh_end(har))
        {
            k = kh_put(m16, har, key, &absent);
            hidv_k * hidv = (hidv_k *) calloc(1, sizeof(hidv_k));
            hidv->count = 0;
            hidv->idfcount = idfcount;
            hidv->mark = 0;
            kh_value(har, k) = hidv;
        }
        arrstrFree(StrArr);
        free(StrArr);
        ks_free(&line);
    }
    return maxidfcount;
}

int pushBarcodeFromFile(kh_m8_t * har, char * fi)
{
    FILE * f = fopen(fi, "r");
    if(f == NULL) error("Barcode file error\n");
    kstring_t line = {0, 0, NULL};
    khint_t k;
    uint64_t key = 0;
    char * over;
    int absent;
    size_t bccount = 0;
    size_t q = 0;
    size_t count = 0;
    while(kgetline(&line, (kgets_func *) fgets, (FILE *) f) != EOF)
    {
        arrstr_k * StrArr = splitStrByChar(line.s, '\t');
        bccount = atoi(StrArr->mem + StrArr->offset_ar[1]);
        key = strtoul(StrArr->mem + StrArr->offset_ar[0], &over, 10);
        q = atoi(StrArr->mem + StrArr->offset_ar[2]);
        count = atoi(StrArr->mem + StrArr->offset_ar[3]);
        k = kh_get(m8, har, key);
        if(k != kh_end(har))
        {
            fhav_k * fhav = kh_value(har, k);
            fhav->count = count;
        } 
        else
        {
            k = kh_put(m8, har, key, &absent);
            fhav_k * fhav = (fhav_k *) calloc(1, sizeof(fhav_k));
            fhav->rec = (arwl_k *) calloc(1, sizeof(arwl_k));
            addEleToArwlXC(fhav->rec, bccount);
            fhav->mist = (arwl_k *) calloc(1, sizeof(arwl_k));
            addEleToArwlXC(fhav->mist, 0);
            fhav->q = (arwl_k *) calloc(1, sizeof(arwl_k));
            addEleToArwlXC(fhav->q, q);
            fhav->count = count;
            kh_val(har, k) = fhav;
        }
        arrstrFree(StrArr);
        free(StrArr);
        ks_free(&line);
    }
    return 1;
}

int addCountForBarcode(kh_m8_t * har, uint64_t bar)
{
    khint_t k;
    k = kh_get(m8, har, bar);
    if(k != kh_end(har))
    {
        // printf("bar : %lu\n", bar);
        fhav_k * fhav = kh_value(har, k);
        fhav->count++;
        return 0;
    }
    return 1;
}

size_t queryIdentifier(idvp_k v, kh_m16_t * har, size_t * count, size_t idflen)
{
    hidv_k * ret;
    size_t result = SIZE_MAX;
    pthread_mutex_lock(&locker);
    khint_t k = kh_get(m16, har, v.cv);
    int absent;
    if(k != kh_end(har))
    {
        ret = kh_value(har, k);
        // if(ret->mark != 1)
        // {
        //     if(addEleToArwlXC(ret->svar, v.cv) == 1)
        //     {
        //         if(ret->svar->l == idflen) ret->mark = 1;
        //     }
        // }
        result = ret->idfcount;
        ret->count++;
    }
    else
    {
        // printf("New Iden\n");
        k = kh_put(m16, har, v.cv, &absent);
        hidv_k * hidv = (hidv_k *) calloc(1, sizeof(hidv_k));
        hidv->count = 1;
        // hidv->svar = (arwl_k *) calloc(1, sizeof(arwl_k));
        // hidv->svar->ar = (size_t *) calloc(1, sizeof(size_t));
        // addEleToArwlRep(hidv->svar, v.sv);
        hidv->idfcount = ++(*count);
        hidv->mark = 0;
        kh_value(har, k) = hidv;
        result = hidv->idfcount;
    }
    pthread_mutex_unlock(&locker);
    return result;
}

void idfHashDestroy(kh_m16_t * har)
{
    if(har)
    {
        // khint_t k;
        // for(k = 0; k < kh_end(har); ++k)
        // {
        //     if(kh_exist(har, k))
        //     {
        //         hidv_k * curhidv = kh_value(har, k);
        //         hidvFree(curhidv);
        //         free(curhidv);
        //     }
        // }
        kh_destroy_m16(har);
    }
}

void destroyMistComArray()
{
    for(int i = 0; i < MAXREADLENGTH; i++)
    {
        kaFree(MistFlagArray[i]);
        free(MistFlagArray[i]);
        for(int j = 0; j < MAXMISMATCH; j++)
        {
            kaFree(MistCombinationArray[i][j]);
            free(MistCombinationArray[i][j]);
            kaFree(MistQarArray[i][j]);
            free(MistQarArray[i][j]);
        }
    }
}
// void countIdenNumber(kh_m16_t * har)
// {
//     if(har)
//     {
//         khint_t k;
//         for(k = 0; k < kh_end(har); ++k)
//         {
//             if(kh_exist(har, k))
//             {
//                 hidv_k * curhidv = kh_value(har, k);
//                 printf(curhidv->odar->count);
//             }
//         }
//     }
// }

// int main(int argc, char const *argv[])
// {
//     char * seq = "GTGTCAAGCACCGCT";
//     initMistComArray();
//     pkarwl_k * k_h = encodeSeqForHash(seq, 2);
//     printf("k_h : %ld\n", k_h->l);
//     for(int i = 0; i < k_h->l; i++)
//     {
//         printf("%lu\n", k_h->ar[i].v);
//     }
// }
//     // sstack_k * Stack = (sstack_k *) calloc(1, sizeof(sstack_k));
//     // initNumStack(Stack);
//     // arwl_k * Result = (arwl_k *) calloc(1, sizeof(arwl_k));
//     // combinationFromAToB(1, 5, 1, Stack, Result);
//     // for(int i = 0; i < Result->l; i++)
//     // {
//     //     printf("Number: %ld\n", Result->ar[i]);
//     // }
//     // seqMistPKXC(13410, 5);
//     // initMistComArray();
//     // char * TAG = "TCATGTCTTCCGATCT";

//     char * TAG1 = "ABCDEFG";
//     char * TAG2 = "MNPPK";
//     for(int i = 0; i < 1000000000; i++)
//     {
//         printf("i: %d\n", i);
//         char * tmp = strConcate(TAG1, TAG2, "-");
//         printf("%s\n", tmp);
//         free(tmp);
//     }
//     // while (1)
//     // {
//     //     pkarwl_k * k_h = encodeSeqForReadMis(TAG, 2);
//     //     // printf("k_h Length: %ld\n", k_h->l);
//     //     // for(int i = 0; i < k_h->l; i++){
//     //     // //     // char * seq = Decode(k_h->ar[i].v, k_h->qlar[i]);
//     //     //     // printf("qlar[%d]: %ld\n", i, k_h->qlar[i]);
//     //     // //     // for(int j = 0; j < 3; j++) 
//     //     // //     // {
//     //     // //         // for( int m = 0; m < k_h->qar->ar[j]; m++)
//     //     // //         // {
//     //     //             printf("Q: %ld \t ", k_h->ar[i].q);
//     //     //             printf("V under Q: %lu\n", k_h->ar[i].v);
//     //     // //         // }
//     //     // }

//     //     //     // free(seq);
//     //     // }
//     //     kpFree(k_h);
//     //     free(k_h);
//     // }
//     // destroyMistComArray();
// }
    // printf("TAG:%s\n", TAG);
    // while (1)
    // {
    //     /* code */
    //     // char * flag = calloc(1, sizeof(char));
    //     // printf("%p\n", flag);

    //     arrstr_k * sar = splitStrByStr(TAG, "CN");
    //                     // free(flag);
    //     // for(int i = 0; i < sar->l; i++)
    //     // {
    //     //     printf("%s\n", sar->ar[i]);
    //     // }
    //     arrstrFree(sar);
    //     free(sar);

    // }
    
//     return 0;
// }

// int main()
// {
//     char 
//     size_t s = 5, mist = 2;
//     pkarwl_k * tmp = encodeSeq(s, 2);
//     for(int i = 0; i < tmp->l; i++)
//     {
//         printf("%ld\n", tmp->ar[i]);
//     }
// }

