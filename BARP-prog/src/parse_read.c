#include "parse_config.h"
#include "parse_read.h"
#include <ctype.h>
#include "kutils.h"
#include <pthread.h>
#include <fcntl.h>
#include <getopt.h>
#include <sys/stat.h>

// barcode kind unit like DPM RPM ODD
#define BUNIT 4
// barcode name unit like DPM-1 DPM-2
#define NUM_THREAD 5
#define PACKAGE_VERSION "0.0.1"
#define OUTPUTTPYENUM 16
#define BUFFER_SIZE 1048576


extern char * optarg;
extern int optind, opterr, optopt;

typedef struct FragTypeCounter
{
    size_t DNA, RNA, NDNR, FullDNA, FullRNA, FullNDNR;
} FragTypeStat;

typedef struct DoubleArrstr
{
    arrstr_k * Read1;
    arrstr_k * Read2;
    arrstr_k * Read3;
    arrstr_k * Read4;
    arwl_k * fb;
    arwl_k * w1, * w2, * w3, * w4;
}darr_k;

typedef struct ParseSeqStruct
{

    fqpool_k fqpool;
    FragTypeStat * ftc;
    // all barcode identify
    int fb, torder, pread;
    int dropN, min_qual;
    int barmode, round, readwrite;
    FILE * fp[OUTPUTTPYENUM];
    barcl_k * barcl;
    arwl_k * idfcount;
    kh_m16_t ** har;
    arwl_k ** barstat;
    darr_k * DoubleRead;
    barc_k * barcounter;
} parseq_k;

typedef struct DoubleArwl
{
    size_t l;
    size_t * a, * b;
} darwl_k;

typedef struct FatherSonStat
{
    darwl_k son;
    size_t father;
}fass_k;

pthread_mutex_t lock;

int FqseqClear(fqseq_k fqseq)
{
    if(fqseq.rds.s) free(fqseq.rds.s);
    if(fqseq.rdq.s) free(fqseq.rdq.s);
    // if(fqseq.rdn.s) free(fqseq.rdn.s);
    // if(fqseq.rdc.s) free(fqseq.rdc.s);
    return 0;
}

void fqseqCopy(fqseq_k * s, kseq_t * m)
{
    memmove(&s->rdc, &m->comment, sizeof(kstring_t));
    memmove(&s->rdn, &m->name, sizeof(kstring_t));
    memmove(&s->rdq, &m->qual, sizeof(kstring_t));
    memmove(&s->rds, &m->seq, sizeof(kstring_t));
    // printf("%s\n", m->name.s);
    s->rdc.s = (char *) calloc(m->comment.l + 1, sizeof(char));
    memmove(s->rdc.s, m->comment.s, m->comment.l);
    s->rdn.s = (char *) calloc(m->name.l + 1, sizeof(char));
    memmove(s->rdn.s, m->name.s, m->name.l);
    s->rds.s = (char *) calloc(m->seq.l + 1, sizeof(char));
    memmove(s->rds.s, m->seq.s, m->seq.l);
    s->rdq.s = (char *) calloc(m->qual.l + 1, sizeof(char));
    memmove(s->rdq.s, m->qual.s, m->qual.l);
}

static int Usage()
{
    fprintf(stdout, "\n");
    fprintf(stdout, "Usage: BARP idb -1 [Read 1] -2 [Read 2] -c [Config file] -@ [Thread] -o [Prefix] --mode [default]\n"
    "Options:\n"
    "  -1           FILE    Fastq read 1\n"
    "  -2           FILE    Fastq read 2\n"
    "  -3           FILE    Fastq read 3\n"
    "  -4           FILE    Fastq read 4\n"
    "  -c           FILE    Config file include barcode chain and white list\n"
    "  -@           INT     Thread number [1]\n"
    "  -o           STR     Prefix of output file\n"
    "  -s           INT     Chunk size [10000]\n"
    "  -a                   Output all reads instead of fully barcoded read\n"
    "  --dropN      STR     Treat Seq barcode with N as Not_Found\n"
    "  --mode       STR     Mode when collision happens: \"default\", \"rnat\", \"atact\", \"minq\" [\"default\"]\n"
    "  --min-mistq  STR     Minimal mismatched base quality to consider the base as sequencing error\n"
    );
    fprintf(stdout, "\n");
    return -1;
}

int openFastq(fqhd_k * fqh)
{
    // printf("Openning FASTQ file.\n");
    fqh->rd1 = fqh->Read1List->mem + fqh->Read1List->offset_ar[fqh->nth];
    fqh->rd2 = fqh->Read2List->mem + fqh->Read2List->offset_ar[fqh->nth];
    fqh->rd3 = fqh->Read3List->mem + fqh->Read3List->offset_ar[fqh->nth];
    fqh->rd4 = fqh->Read4List->mem + fqh->Read4List->offset_ar[fqh->nth];

    for(int i = 0; i < 4; i++)
    {
        size_t Flag = 1 << i;
        if((fqh->flag & Flag) == Flag)
        {
            switch(i)
            {
                case 0:
                {
                    if(fqh->nth > 0) gzclose(fqh->r1);
                    fqh->r1 = gzopen(fqh->rd1, "r");
                    if(fqh->r1) 
                    {
                        if(fqh->nth > 0) kseq_destroy(fqh->ks1);
                        fqh->ks1 = kseq_init(fqh->r1); 
                    }
                    else 
                    {
                        fprintf(stderr,ANSI_COLOR_RED "[error]  " ANSI_COLOR_RESET ANSI_COLOR_MAGENTA  "No such file or file error: %s."ANSI_COLOR_RESET "\n", fqh->rd1);
                        exit(1);
                    }
                    break;
                }
                case 1:
                {
                    if(fqh->nth > 0) gzclose(fqh->r2);
                    fqh->r2 = gzopen(fqh->rd2, "r");
                    if(fqh->r2) 
                    {
                        if(fqh->nth > 0) kseq_destroy(fqh->ks2);
                        fqh->ks2 = kseq_init(fqh->r2); 
                    }
                    else 
                    {
                        fprintf(stderr,ANSI_COLOR_RED "[error]  " ANSI_COLOR_RESET ANSI_COLOR_MAGENTA  "No such file or file error: %s."ANSI_COLOR_RESET "\n", fqh->rd2);
                        exit(1);
                    }
                    break;
                }
                case 2:
                {
                    if(fqh->nth > 0) gzclose(fqh->r3);
                    fqh->r3 = gzopen(fqh->rd3, "r");
                    if(fqh->r3) 
                    {
                        if(fqh->nth > 0) kseq_destroy(fqh->ks3);
                        fqh->ks3 = kseq_init(fqh->r3); 
                    }
                    else 
                    {
                        fprintf(stderr,ANSI_COLOR_RED "[error]  " ANSI_COLOR_RESET ANSI_COLOR_MAGENTA  "No such file or file error: %s."ANSI_COLOR_RESET "\n", fqh->rd3);
                        exit(1);
                    }
                    break;
                }
                case 3:
                {
                    if(fqh->nth > 0) gzclose(fqh->r4);
                    fqh->r4 = gzopen(fqh->rd4, "r");
                    if(fqh->r4) 
                    {
                        if(fqh->nth > 0) kseq_destroy(fqh->ks4);
                        fqh->ks4 = kseq_init(fqh->r4); 
                    }
                    else 
                    {
                        fprintf(stderr,ANSI_COLOR_RED "[error]  " ANSI_COLOR_RESET ANSI_COLOR_MAGENTA  "No such file or file error: %s."ANSI_COLOR_RESET "\n", fqh->rd4);
                        exit(1);
                    }
                    break;
                }
                default:
                {
                    break;
                }
            }
        }
    }
    return 0;
}

fqpool_k loadFastq(fqhd_k * fqh, size_t PoolSize)
{
    fqpool_k ret = {0};
    ret.fqrd = (fqrd_k *) calloc(PoolSize, sizeof(fqrd_k));;
    int r1l, r2l, r3l, r4l;
    size_t count = 0;
    fqrd_k * Curfqrd = ret.fqrd;

    LoadLoop:
    {
        while(count < PoolSize)
        {
            // printf("CUR : %p\n", Curfqrd);
            if(fqh->ks1 != NULL) 
            {
                if((r1l = kseq_read(fqh->ks1)) >= 0) fqseqCopy(&Curfqrd->rd1, fqh->ks1); else break;
            }
            if(fqh->ks2 != NULL) 
            {
                if((r2l = kseq_read(fqh->ks2)) >= 0) fqseqCopy(&Curfqrd->rd2, fqh->ks2); else break;
            }
            if(fqh->ks3 != NULL) 
            {
                if((r3l = kseq_read(fqh->ks3)) >= 0) fqseqCopy(&Curfqrd->rd3, fqh->ks3); else break;
            }
            if(fqh->ks4 != NULL) 
            {
                if((r4l = kseq_read(fqh->ks4)) >= 0) fqseqCopy(&Curfqrd->rd4, fqh->ks4); else break;
            }
            Curfqrd++;
            count++;
        }
        if(count < PoolSize)
        {
            if(fqh->nth == fqh->Read1List->l - 1)
            {
                goto LoadOver;
            } 
            else
            {
                fqh->nth++;
                if(openFastq(fqh)) 
                {
                    error("Open FASTQ file failed.");
                }
                else
                {
                    goto LoadLoop;
                }
            }
        }
        else
        {
            goto LoadOver;
        }
    }

    LoadOver:
    // printf("Load Over\n");
    // for(int i = 0; i < count; i++)
    // {
    //     printf("%s\n", ret.fqrd->rd1.rdn.s);
    // }
    ret.count = count;
    // printf("RET : %p\n", ret.fqrd);
    return ret;
}

void fqseqPoolFree(fqpool_k seqPool)
{
    // printf("Free Memory: %ld\n", seqPool.count);
    for(int i = 0; i < seqPool.count; i++)
    {
        if(seqPool.fqrd[i].rd1.rdc.s) free(seqPool.fqrd[i].rd1.rdc.s);
        if(seqPool.fqrd[i].rd1.rdn.s) free(seqPool.fqrd[i].rd1.rdn.s);
        if(seqPool.fqrd[i].rd1.rds.s) free(seqPool.fqrd[i].rd1.rds.s);
        if(seqPool.fqrd[i].rd1.rdq.s) free(seqPool.fqrd[i].rd1.rdq.s);
        if(seqPool.fqrd[i].rd2.rdc.s) free(seqPool.fqrd[i].rd2.rdc.s);
        if(seqPool.fqrd[i].rd2.rdn.s) free(seqPool.fqrd[i].rd2.rdn.s);
        if(seqPool.fqrd[i].rd2.rdq.s) free(seqPool.fqrd[i].rd2.rdq.s);
        if(seqPool.fqrd[i].rd2.rds.s) free(seqPool.fqrd[i].rd2.rds.s);
        if(seqPool.fqrd[i].rd3.rdc.s) free(seqPool.fqrd[i].rd3.rdc.s);
        if(seqPool.fqrd[i].rd3.rdn.s) free(seqPool.fqrd[i].rd3.rdn.s);
        if(seqPool.fqrd[i].rd3.rdq.s) free(seqPool.fqrd[i].rd3.rdq.s);
        if(seqPool.fqrd[i].rd3.rds.s) free(seqPool.fqrd[i].rd3.rds.s);
        if(seqPool.fqrd[i].rd4.rdc.s) free(seqPool.fqrd[i].rd4.rdc.s);
        if(seqPool.fqrd[i].rd4.rdn.s) free(seqPool.fqrd[i].rd4.rdn.s);
        if(seqPool.fqrd[i].rd4.rdq.s) free(seqPool.fqrd[i].rd4.rdq.s);
        if(seqPool.fqrd[i].rd4.rds.s) free(seqPool.fqrd[i].rd4.rds.s);
    }
    seqPool.count = 0;
    free(seqPool.fqrd);
}

static inline void mergeFile(char * input, int fo)
{
    size_t ReadLen;
    unsigned char buffer[BUFFER_SIZE];
    int fi = open(input, O_RDONLY);
    while (1)
    {
        ReadLen = read(fi, buffer, BUFFER_SIZE);
        if(ReadLen > 0)
        {
            write(fo, buffer, ReadLen);
        }
        else
        {
            break;
        }
    }
    close(fi);
}

// because slash in name will be cutted from star
static inline void removeSlashInName(char * s)
{
    char * ptr = s;
    while(1)
    {
        if(* ptr == '/')
        {
            *ptr = '\0';
            break;
        } 
        if(* ptr == '\0') break;
    }
}

static inline int writeToResult(arrstr_k * s, pfqrd_k ret, fqrd_k * fqrd, kstring_t * name, int read, int fb)
{
    switch(read)
    {
        case 1:
        {
            char * nptr;
            nptr = ret.rd1.rdn;
            * nptr++ = '@';
            strcpy(nptr, fqrd->rd1.rdn.s);
            nptr = nptr + fqrd->rd1.rdn.l;
            if(name->l > 0)
            {
                removeSlashInName(nptr);
                strcpy(nptr, "|||");
                nptr = nptr + 3;
                strcpy(nptr, name->s);
                nptr = nptr + name->l;   
            }
            *nptr = '\0';
            addStrToArrstrRepXC(s, ret.rd1.rdn);
            // write seq
            if(fb == 1) 
            {
                addStrToArrstrRepXC(s, ret.rd1.rds);
            }
            else
            {
                addStrToArrstrRepXC(s, fqrd->rd1.rds.s);
            }
            // write comment
            nptr = ret.rd1.rdc;
            *nptr++ = '+';
            strcpy(nptr, fqrd->rd1.rdc.s);
            *nptr = '\0';
            addStrToArrstrRepXC(s, ret.rd1.rdc);
            if(fb == 1) 
            {
                addStrToArrstrRepXC(s, ret.rd1.rdq);
            }
            else
            {
                // write quality
                addStrToArrstrRepXC(s, fqrd->rd1.rdq.s);
            }
            break;
        }
        case 2:
        {
            char * nptr;
            nptr = ret.rd2.rdn;
            * nptr++ = '@';
            strcpy(nptr, fqrd->rd2.rdn.s);
            nptr = nptr + fqrd->rd2.rdn.l;
            if(name->l > 0)
            {
                removeSlashInName(nptr);
                strcpy(nptr, "|||");
                nptr = nptr + 3;
                strcpy(nptr, name->s);
                nptr = nptr + name->l;   
            }
            *nptr = '\0';
            addStrToArrstrRepXC(s, ret.rd2.rdn);
            // write seq
            if(fb == 1) 
            {
                addStrToArrstrRepXC(s, ret.rd2.rds);
            }
            else
            {
                addStrToArrstrRepXC(s, fqrd->rd2.rds.s);
            }
            // write comment
            nptr = ret.rd2.rdc;
            *nptr++ = '+';
            strcpy(nptr, fqrd->rd2.rdc.s);
            *nptr = '\0';
            addStrToArrstrRepXC(s, ret.rd2.rdc);
            if(fb == 1) 
            {
                addStrToArrstrRepXC(s, ret.rd2.rdq);
            }
            else
            {
                // write quality
                addStrToArrstrRepXC(s, fqrd->rd2.rdq.s);
            }
            // printf("Name : %s\n", ret.rd2.rdn);
            break;
        }
        case 3:
        {
            char * nptr;
            nptr = ret.rd3.rdn;
            * nptr++ = '@';
            strcpy(nptr, fqrd->rd3.rdn.s);
            nptr = nptr + fqrd->rd3.rdn.l;
            if(name->l > 0)
            {
                removeSlashInName(nptr);
                strcpy(nptr, "|||");
                nptr = nptr + 3;
                strcpy(nptr, name->s);
                nptr = nptr + name->l;   
            }
            *nptr = '\0';
            addStrToArrstrRepXC(s, ret.rd3.rdn);
            // write seq
            if(fb == 1) 
            {
                addStrToArrstrRepXC(s, ret.rd3.rds);
            }
            else
            {
                addStrToArrstrRepXC(s, fqrd->rd3.rds.s);
            }
            // write comment
            nptr = ret.rd3.rdc;
            *nptr++ = '+';
            strcpy(nptr, fqrd->rd3.rdc.s);
            *nptr = '\0';
            addStrToArrstrRepXC(s, ret.rd3.rdc);
            if(fb == 1) 
            {
                addStrToArrstrRepXC(s, ret.rd3.rdq);
            }
            else
            {
                // write quality
                addStrToArrstrRepXC(s, fqrd->rd3.rdq.s);
            }
            break;
        }
        case 4:
        {
            char * nptr;
            nptr = ret.rd4.rdn;
            * nptr++ = '@';
            strcpy(nptr, fqrd->rd4.rdn.s);
            nptr = nptr + fqrd->rd4.rdn.l;
            if(name->l > 0)
            {
                removeSlashInName(nptr);
                strcpy(nptr, "|||");
                nptr = nptr + 3;
                strcpy(nptr, name->s);
                nptr = nptr + name->l;   
            }
            *nptr = '\0';
            addStrToArrstrRepXC(s, ret.rd4.rdn);
            // write seq
            if(fb == 1) 
            {
                addStrToArrstrRepXC(s, ret.rd4.rds);
            }
            else
            {
                addStrToArrstrRepXC(s, fqrd->rd4.rds.s);
            }
            // write comment
            nptr = ret.rd4.rdc;
            *nptr++ = '+';
            strcpy(nptr, fqrd->rd4.rdc.s);
            *nptr = '\0';
            addStrToArrstrRepXC(s, ret.rd4.rdc);
            if(fb == 1) 
            {
                addStrToArrstrRepXC(s, ret.rd4.rdq);
            }
            else
            {
                // write quality
                addStrToArrstrRepXC(s, fqrd->rd4.rdq.s);
            }
            break;
        }
        default:
            break;
    }
    return 0;
}

void * parseSeq(void * arg)
{
    parseq_k SeqConfig = *((parseq_k *) arg);
    
    int a = SeqConfig.fb;
    barcl_k * barcl = SeqConfig.barcl;
    kh_m16_t ** hart = SeqConfig.har;
    arwl_k * idfcount = SeqConfig.idfcount;
    arwl_k ** BarStat = SeqConfig.barstat;
    int min_qual = SeqConfig.min_qual;
    int dropN = SeqConfig.dropN;
    char tmpbar[600] = {0};
    char tmpqual[600] = {0};
    darr_k * DoubleRead = SeqConfig.DoubleRead;
    FragTypeStat * FTC = SeqConfig.ftc;
    arrstr_k * ResultPoolR1 = DoubleRead->Read1;
    arrstr_k * ResultPoolR2 = DoubleRead->Read2;
    arrstr_k * ResultPoolR3 = DoubleRead->Read3;
    arrstr_k * ResultPoolR4 = DoubleRead->Read4;
    int TMPReadWrite = SeqConfig.readwrite;
    int BarMode = SeqConfig.barmode;
    int Round = SeqConfig.round;
    // a barcode order; b barcode counter
    barc_k * BarcodeCounter = SeqConfig.barcounter;
    for(int i = 0; i < barcl->l; i++) BarStat[i] = (arwl_k *) calloc(1, sizeof(arwl_k));
    // printf("Here\n");
    // printf("Parsing Seq : %ld\n", SeqConfig.fqpool.count);
    for(int p = 0; p < SeqConfig.fqpool.count; p++)
    {
        pfqrd_k ret = {0};
        fqrd_k * fqrd = &SeqConfig.fqpool.fqrd[p];
        fqseq_k curseq;
        fqseqstr_k * tmpseq;
        arrstr_k * Result = (arrstr_k *) calloc(1, sizeof(arrstr_k));
        ret.fb = 1;
        // var used to identify if reads belong to one identifier; bvar to check the read is barcoded by barcode of identifier.
        size_t q = 0;
        barch_k * curbarch;
        arwl_k * var, * bvar;
        if(barcl->idf)
        {
            var = (arwl_k *) calloc(1, sizeof(arwl_k));
            bvar = (arwl_k *) calloc(1, sizeof(arwl_k));
        }
        // write read to Read 1 , 2 ,3...
        int ReadWrite = 0;
        for(int i = 0; i < barcl->l; i++)
        {
            arwl_k * CurBarStat = BarStat[i];
            size_t CurBarchainStat = 0;
            curbarch = &barcl->barch[i];
            ReadWrite = ReadWrite | curbarch->rw;
            arrstr_k * CurbarChainResult = (arrstr_k *) calloc(1, sizeof(arrstr_k));

            // current barcode chain situation: var: barcode type; bvar: barcode.
            arwl_k * curbarvar, * curbarbvar;
            if(barcl->idf)
            {
                curbarvar = (arwl_k *) calloc(1, sizeof(arwl_k));
                curbarbvar = (arwl_k *) calloc(1, sizeof(arwl_k));
            }

            switch(curbarch->rdp.rd)
            {
                case 1:
                {
                    // printf("Read 1: %s\n", curseq.rds.s);
                    curseq = fqrd->rd1;
                    tmpseq = &ret.rd1;
                    break;
                }
                case 2:
                {  
                    curseq = fqrd->rd2;
                    tmpseq = &ret.rd2;
                    break;
                }
                case 3:
                {
                    curseq = fqrd->rd3;
                    tmpseq = &ret.rd3;
                    break;
                }
                case 4:
                {
                    curseq = fqrd->rd4;
                    tmpseq = &ret.rd4;
                    break;
                }
                default:
                {
                    fprintf(stderr, ANSI_COLOR_RED "[error]  " ANSI_COLOR_RESET ANSI_COLOR_MAGENTA "Currently SCSM only supports for pair end and single end sequence!" ANSI_COLOR_RESET "\n");
                    exit(1);
                    break;
                }
            }
            // printf("%d\n", SeqConfig.torder);
            char * sptr = curseq.rds.s, * qptr = curseq.rdq.s;
            sptr = movePtrAlongStr(sptr, curseq.rds.s, curbarch->rdp.pos - 1);
            qptr = movePtrAlongStr(qptr, curseq.rdq.s, curbarch->rdp.pos - 1);
            if(sptr && qptr)
            {
                // printf("%s - %s\n", SeqConfig.fqpool.fqrd[0].rd1.rdn.s, sptr);
                for(int j = 0; j < curbarch->l; j++)
                {
                    char * presptr = sptr, * preqptr = qptr;
                    // barcode length.
                    size_t zcp = 0;
                    char * CurbarResult = (char *) calloc(1, 10 * sizeof(char));
                    char NOT_FOUND[] = "NOT_FOUND";
                    memmove(CurbarResult, NOT_FOUND, 10);
                    size_t RecResult;
                    bar_k curbar;
                    size_t EndGenome = 0;
                    for(int m = 0; m < curbarch->barh[j].l; m++)
                    {
                        curbar = curbarch->barh[j].bar[m];
                        sptr = movePtrAlongStr(presptr, curseq.rds.s, curbar.sp);
                        qptr = movePtrAlongStr(preqptr, curseq.rdq.s, curbar.sp);
                        RecResult = SIZE_MAX;
                        if(sptr && qptr)
                        {
                            size_t z = curbar.minlen;
                            // current barcode is not genome.
                            if(curbar.ord != SIZE_MAX)
                            {
                                int k = 0;
                                for(k = 0; k <= curbar.la; k++)
                                {
                                    if(sptr && (strlen(sptr) >= z))
                                    {
                                        for(z = curbar.minlen; z <= curbar.maxlen; z++)
                                        {
                                            zcp = z;
                                            // partial sequence of the read 1
                                            strncpy(tmpbar, sptr, z * sizeof(char));
                                            tmpbar[z] = '\0';
                                            if(curbar.bighash || min_qual > 0 || BarMode > 0) 
                                            {
                                                strncpy(tmpqual, qptr, z * sizeof(char));
                                                tmpqual[z] = '\0';
                                            }

                                            if(curbar.hash)
                                            {
                                                bandrec_k QueryRet;
                                                // dense barcode white list 
                                                if(curbar.bighash || BarMode > 0) 
                                                {
                                                    // first round 
                                                    if(Round == 0 && BarMode == 2)
                                                    {
                                                        QueryRet = querySeqBestC(tmpbar,tmpqual,curbar.hash, 0, min_qual, dropN, BarMode);
                                                    }
                                                    else
                                                    {
                                                        // second round 
                                                        QueryRet = querySeqBestC(tmpbar,tmpqual,curbar.hash, curbar.mist, min_qual, dropN, BarMode);
                                                    }
                                                    RecResult = QueryRet.rec;
                                                }
                                                // sparse barcode white list. 
                                                else 
                                                {
                                                    if(curbar.random == 1)
                                                    {
                                                        // query random barcode; only support sparse barcode. 
                                                        QueryRet = querySeqBestD(tmpbar, tmpqual, curbar.hash, curbar.mist, min_qual, dropN, curbar.rec);
                                                        RecResult = QueryRet.rec;
                                                    }
                                                    else
                                                    {
                                                        QueryRet = querySeqBestB(tmpbar, tmpqual, curbar.hash, curbar.mist, min_qual, dropN);
                                                        RecResult = QueryRet.rec;
                                                    }
                                                }
                                                
                                                if(RecResult != SIZE_MAX)
                                                {
                                                    if(Round != 1) addBarToBackRep(BarcodeCounter, QueryRet.bar);
                                                    if(curbar.random == 1) 
                                                    {
                                                        CurbarResult = realloc(CurbarResult, (z + 1 + strlen(curbar.na) + 1) * sizeof(char));
                                                        strcpy(CurbarResult, curbar.na);
                                                        strcpy(CurbarResult + strlen(curbar.na), ":");
                                                        strcpy(CurbarResult + strlen(curbar.na) + 1, tmpbar);
                                                        CurbarResult[z + strlen(curbar.na) + 1] = '\0';
                                                    }
                                                    else
                                                    {
                                                        CurbarResult = (char *) realloc(CurbarResult, (strlen(curbar.rec->mem + curbar.rec->offset_ar[RecResult]) + 1) * sizeof(char));
                                                        strcpy(CurbarResult, curbar.rec->mem + curbar.rec->offset_ar[RecResult]);
                                                        CurbarResult[strlen(curbar.rec->mem + curbar.rec->offset_ar[RecResult])] = '\0';
                                                    }
                                                    z++;
                                                    zcp = z;
                                                    break;
                                                }
                                            }
                                            else
                                            {
                                                // UMI ...
                                                free(CurbarResult);
                                                CurbarResult = strConcate(curbar.na, tmpbar, ":");
                                                z++;
                                                zcp = z;
                                                RecResult = SIZE_MAX - 1;
                                                break;
                                            }
                                        }
                                    }
                                    else
                                    {
                                        // fprintf(stderr, ANSI_COLOR_YELLOW "[warnings]  Barcode: %s position exceeded read length at %s! Considered as NOT_FOUND!" ANSI_COLOR_RESET "\n", curbar.na, curseq.rdn.s);
                                        break;
                                    }

                                    if(RecResult != SIZE_MAX) break;
                                    if(k != curbar.la)
                                    {
                                        sptr = movePtrAlongStr(sptr, curseq.rds.s, 1);
                                        qptr = movePtrAlongStr(qptr, curseq.rdq.s, 1);
                                    }
                                }
                                
                                if(RecResult != SIZE_MAX) goto HAVEHASH;
                                if(m == curbarch->barh[j].l - 1) goto HAVEHASH;
                            }
                            else
                            {
                                if(curbar.minlen != SIZE_MAX)
                                {
                                    tmpseq->sl = curbar.minlen;
                                    tmpseq->ql = curbar.minlen;
                                    strncpy(tmpseq->rds, sptr, curbar.minlen);
                                    tmpseq->rds[curbar.minlen] = '\0';
                                    strncpy(tmpseq->rdq, qptr, curbar.minlen);
                                    tmpseq->rdq[curbar.minlen] = '\0';
                                    sptr = movePtrAlongStr(sptr, curseq.rds.s, curbar.minlen);
                                    qptr = movePtrAlongStr(qptr, curseq.rdq.s, curbar.minlen);
                                }
                                else
                                {
                                    tmpseq->sl = strlen(sptr);
                                    tmpseq->ql = strlen(sptr);
                                    strncpy(tmpseq->rds, sptr, strlen(sptr));
                                    tmpseq->rds[strlen(sptr)] = '\0';
                                    strncpy(tmpseq->rdq, qptr, strlen(qptr));
                                    tmpseq->rdq[strlen(qptr)] = '\0';
                                    sptr = movePtrAlongStr(sptr, curseq.rds.s, strlen(sptr));
                                    qptr = movePtrAlongStr(qptr, curseq.rdq.s, strlen(qptr));
                                }
                                goto NOTHAVEHASH;
                            }
                        }
                        else
                        {
                            if(curbar.minlen > 0) 
                            {
                                if(curbar.hash)
                                {
                                    RecResult = SIZE_MAX;
                                }
                                else
                                {
                                    RecResult = 1;
                                    EndGenome = 1;
                                }
                                // warnings("Barcode chain %d start position exceeded read length at %s!", i, curseq.rdn.s);
                                continue;
                            }
                            else
                            {
                                switch (curbarch->rdp.rd)
                                {
                                    // means the genome do not have read length, do not write.
                                    case 1:
                                        ReadWrite = ReadWrite & (~1);
                                        break;
                                    case 2:
                                        ReadWrite = ReadWrite &  (~(1 << 1));
                                        break;
                                    case 3:
                                        ReadWrite = ReadWrite &  (~(1 << 2));
                                        break;
                                    case 4:
                                        ReadWrite = ReadWrite &  (~(1 << 3));
                                        break;
                                    default:
                                        break;
                                }
                                goto NOTHAVEHASH;
                            }
                        }
                    }

                    HAVEHASH:
                    {
                        sptr = movePtrAlongStr(sptr, curseq.rds.s, zcp-1);
                        qptr = movePtrAlongStr(qptr, curseq.rdq.s, zcp-1);

                        // barcode found
                        if(RecResult != SIZE_MAX)
                        {
                            if(EndGenome != 1) addStrToArrstrRepXC(CurbarChainResult, CurbarResult);
                            // printf("Result:%s - %s - %d - %ld\n",curseq.rdn.s, CurbarResult, j, curbarch->l);
                            free(CurbarResult);
                            CurBarchainStat = (CurBarchainStat << 1) | 1;
                            if(barcl->idf)
                            {
                                addEleToArwlRep(curbarbvar, RecResult);
                                addEleToArwlRep(curbarvar, curbar.ord);
                                if(RecResult != SIZE_MAX) q = (q << 1) | 1; else q = (q << 1) | 0;
                            }

                            if(j < curbarch->l - 1)
                            {
                                // next barcode
                                continue;
                            }
                            // current barcode chain find over, go to the next chain of the same genome or common barcode chain.
                            else
                            {
                                // current barcode chain is not the last barcode chain
                                if( i != barcl->l - 1)
                                {
                                    // if next barcode chain is not the same genome and is not common barcode chain, go to common barcode chain.
                                    if((barcl->barch[i].gtype != barcl->barch[i + 1].gtype) && (barcl->barch[i + 1].gtype != 0))
                                    {
                                        int f = 0;
                                        // go to the common barcode chain
                                        for(f = 0; f < barcl->jump->l; f++)
                                        {
                                            if(barcl->barch[f].gtype == 0)
                                            {
                                                i = f;
                                                break;
                                            }
                                        }
                                        // no common barcode chain
                                        if(f == barcl->jump->l)
                                        {
                                            i = barcl->l;
                                            goto ENDUP;
                                        }
                                        else
                                        {
                                            goto ENDUP;
                                        }
                                    }
                                    else
                                    {
                                        goto ENDUP;
                                    }
                                }
                                else
                                {
                                    goto ENDUP;
                                }
                            }
                        }
                        // barcode not found
                        else
                        {
                            // current barcode chain is not the final barcode chain
                            if(i != barcl->l - 1)
                            {
                                // printf("Not found: %s - %d\n", curseq.rdn.s, i);
                                // if next genome is common genome or this genome is the final genome, as next genome is neither current genome nor common genome, then this genome need to find over.
                                int NextGenome = 0, f = i + 1;
                                for(; f < barcl->l - 1; f++)
                                {
                                    if(barcl->barch[f].gtype != curbarch->gtype || barcl->barch[f].gtype != 0)
                                    {
                                        NextGenome = 1;
                                        break;
                                    }
                                }
                                if(NextGenome == 0)
                                {
                                    ret.fb = 0;
                                    q = q << 1 | 0;
                                    CurBarchainStat = (CurBarchainStat << 1 ) | 0;
                                    addStrToArrstrRepXC(CurbarChainResult, CurbarResult);
                                    free(CurbarResult);
                                    if(barcl->idf)
                                    {
                                        addEleToArwlRep(curbarvar, curbar.ord);
                                        addEleToArwlRep(curbarbvar, RecResult);
                                    }
                                    if(j < curbarch->l - 1) 
                                    {
                                        continue;
                                    }
                                    else
                                    {
                                        goto ENDUP;
                                    }
                                }
                                // next genome is not common genome
                                else
                                {
                                    // go to next genome, the rest barcode chain of the genome is considered as not found.
                                    CurBarchainStat = CurBarchainStat << (curbarch->l - j);
                                    // printf("Seq : %s - CurbarchainStat : %ld\n", curseq.rdn.s, CurBarchainStat);
                                    q = 0;
                                    if(barcl->idf)
                                    {
                                        kaFree(curbarbvar);
                                        free(curbarbvar);
                                        kaFree(curbarvar);
                                        free(curbarvar);
                                        curbarvar = NULL;
                                        curbarbvar = NULL;
                                    }
                                    free(CurbarResult);
                                    CurbarChainResult->l = 0;
                                    i = barcl->jump->ar[i] - 1;
                                    // printf("Go to next genome\n");
                                    // printf("Curbar: %s\n", curbar.na);
                                    // break j: do not loop current barcode chain.
                                    break;
                                }
                            }
                            else
                            {
                                ret.fb = 0;
                                CurBarchainStat = (CurBarchainStat << 1 ) | 0;
                                
                                addStrToArrstrRepXC(CurbarChainResult, CurbarResult);
                                free(CurbarResult);

                                if(barcl->idf)
                                {
                                    addEleToArwlRep(curbarvar, curbar.ord);
                                    addEleToArwlRep(curbarbvar, RecResult);
                                }

                                if(j < curbarch->l - 1) 
                                {
                                    continue;
                                }
                                else
                                {
                                    goto ENDUP;
                                }
                            }
                        }
                    }

                    NOTHAVEHASH:
                    {
                        CurBarchainStat = (CurBarchainStat << 1) | 1;
                        free(CurbarResult);
                        if(j < curbarch->l - 1)
                        {
                            // next barcode
                            continue;
                        }
                        // current barcode chain find over, go to the next chain of the same genome or common barcode chain.
                        else
                        {
                            // current genome is not the last barcode chain
                            if( i != barcl->l - 1)
                            {
                                // if next barcode chain is not the same genome and not common barcode chain, go to common barcode chain.
                                if((barcl->barch[i].gtype != barcl->barch[i + 1].gtype) && (barcl->barch[i + 1].gtype != 0))
                                {
                                    int f = 0;
                                    // go to the common barcode chain
                                    for(f = 0; f < barcl->jump->l; f++)
                                    {
                                        if(barcl->barch[f].gtype == 0)
                                        {
                                            i = f;
                                            break;
                                        }
                                    }
                                    // no common barcode chain
                                    if(f == barcl->jump->l)
                                    {
                                        i = barcl->l;
                                        goto ENDUP;
                                    }
                                    else
                                    {
                                        goto ENDUP;
                                    }
                                }
                            }
                            goto ENDUP;
                        }
                    }
                }
            }
            else
            {
                if(curseq.rdn.s) warnings("Barcode chain %d start position exceeded read length at %s!", i, curseq.rdn.s); else error("Read format wrong or barcode specified at non-exist read.\n");
                // fprintf(stderr, ANSI_COLOR_RED "[error]  " ANSI_COLOR_RESET ANSI_COLOR_MAGENTA "Barcode chain %d start position exceeded read length at %s!" ANSI_COLOR_RESET "\n", i, curseq.rdn.s);
                // exit(1);
            }

            ENDUP:
            {
                // printf("Ending Up: %d - %ld\n", ko, CurBarchainStat);
                addEleToArwlRep(CurBarStat, CurBarchainStat);
                // printf("CurBarChain : %ld\n", CurbarChainResult->l);
                if(CurbarChainResult->l > 0)
                {
                    overLayArrstr(Result, CurbarChainResult);
                }
                else
                {
                    arrstrFree(CurbarChainResult);
                    free(CurbarChainResult);
                }
                
                if(barcl->idf)
                {
                    overLayArwl(var, curbarvar);
                    overLayArwl(bvar, curbarbvar);
                }
            }
        }
        
        kstring_t * name = (kstring_t *) calloc(1, sizeof(kstring_t));
        ks_initialize(name);
        
        for(int i = 0; i < Result->l; i++)
        {
            if(i == Result->l - 1)
            {
                kputs(Result->mem + Result->offset_ar[i], name);
            }
            else
            {
                kputs(Result->mem + Result->offset_ar[i], name);
                kputs("|", name);
            }
        }
        arrstrFree(Result);
        free(Result);
        if(ret.fb == 1 || a == 1)
        {
            int IDFlag = 1;
            // printf("Length : %ld\n", var->l);
            for( int i = 0; i < barcl->idl; i++)
            {
                kh_m16_t * har = hart[i];
                size_t Gflag = 0;
                char NOT_FOUND_1[] = "NOT_FOUND";
                char * CurIdenResult = (char *) calloc(10, sizeof(char));
                memmove(CurIdenResult, NOT_FOUND_1, 10 * sizeof(char));
                id_k curidf = barcl->idf[i];
                char Gtype = '0';
                if(strcmp(curidf.na, "DNA") == 0)
                {
                    Gflag = 1;
                    Gtype = 'D';
                }
                else if(strcmp(curidf.na, "RNA") == 0)
                {
                    Gflag = 1;
                    Gtype = 'R';
                }
                
                if(Gflag == 0)
                {
                    if(var->l < curidf.l)
                    {
                        if(IDFlag == 1)
                        {
                            kputs("|||", name);
                            IDFlag = 0;
                        }
                        else
                        {
                            kputs("|", name);
                        }
                        kputs(CurIdenResult, name);
                        free(CurIdenResult);
                        // if(!Gflag) fprintf(stderr, ANSI_COLOR_YELLOW "[warnings]  Identifier: %s position exceeded barcode chain length at %s! Considered as NOT_FOUND!" ANSI_COLOR_RESET "\n", curidf.na, curseq.rdn.s);
                        continue;
                    }
                }

                size_t curq = curidf.q << (var->l - curidf.l), flag = 1 << (var->l - 1), vpos = 0;
                size_t curres = SIZE_MAX;

                // printf("curq - q: %ld - %ld\n", curq, q);
                if(curidf.anchor == 1) curidf.l = var->l;


                idvp_k v = {0};
                size_t vflag = (((1 << BUNIT) - 1) << (BUNIT * (var->l - 1)));
                int MoveTime = 0;
                size_t VUNIT = barcl->barch[0].maxb;
                for(int j = 0; j < (var->l - curidf.l) + 1; j++)
                {
                    if((curq & q) == curq)
                    {
                        for(int k = 0; k < curidf.l; k++)
                        {
                            if((curq & flag) > 0)
                            {
                                // the barcode hole have more than one kind of barcode : DPM/RPM ov means the rest of the barcode chain.
                                size_t ori = curidf.v->ar[0] & vflag, set = 0;
                                for(int o = 0; o < curidf.v->l; o++)
                                {
                                    if((curidf.v->ar[o] & vflag) != ori)
                                    {
                                        set = 1;
                                    }
                                }
                                MoveTime++;
                                if(set)
                                {
                                    // sv : DPM - Y - ODD -EVEN / RPM - Y - ODD - EVEN
                                    // ov : DPM1 - Y2 - ODD3 - EVEN4 
                                    // cv : Y2 - ODD3 - EVEN4
                                    v.sv = (v.sv << BUNIT) | var->ar[vpos];
                                    // v.ov = (v.ov << VUNIT) | bvar->ar[vpos];
                                }
                                else
                                {
                                    v.cv = (v.cv << VUNIT) | bvar->ar[vpos];
                                    v.sv = (v.sv << BUNIT) | var->ar[vpos];
                                    // v.ov = (v.ov << VUNIT) | bvar->ar[vpos];
                                }
                                if((MoveTime * VUNIT) >= ULLONG_MAX) 
                                {
                                    fprintf(stderr, ANSI_COLOR_RED "[error]  " ANSI_COLOR_RESET ANSI_COLOR_MAGENTA "Identifier exceeded the max capacity!" ANSI_COLOR_RESET "\n");
                                    exit(1);
                                }
                            }
                            flag = flag >> 1;
                            vflag = vflag >> BUNIT;
                            vpos++;
                        }
                        size_t idffb = 0;
                        for(int o = 0; o < curidf.v->l; o++)
                        {
                            if(v.sv == curidf.v->ar[o])
                            {
                                idffb = 1;
                                break;
                            }
                        }
                        // printf("Found Pattern\n");
                        if(idffb > 0)
                        {
                            if(Gflag)
                            {
                                ret.type = Gtype;
                                break;
                            }
                            else
                            {
                                curres = queryIdentifier(v, har, &idfcount->ar[i], curidf.v->l);
                                // curres = 1;
                                break;
                            }
                        }
                    }
                    else
                    {
                        vpos++;
                        flag = flag >> 1;
                        curq = curq >> 1;
                        vflag = vflag >> BUNIT;
                    }
                }
                if(!Gflag)
                {
                    if(curres != SIZE_MAX)
                    {
                        char * number = convertSizetToStr(curres);
                        free(CurIdenResult);
                        CurIdenResult = strConcate(curidf.na, number, "_");
                        free(number);
                    }
                    if(IDFlag == 1)
                    {
                        kputs("|||", name);
                        IDFlag = 0;
                    }
                    else
                    {
                        kputs("|", name);
                    }
                    kputs(CurIdenResult, name);
                }
                free(CurIdenResult);
            }
        }
        // printf("ret.type: %d\n", ret.type);
        if(barcl->idf)
        {
            kaFree(var);
            free(var);
            kaFree(bvar);
            free(bvar);
        }

        int Write = 0;
        pfqrd_k * resfqrd = &ret;
        // is atac ; first round; not fully barcode.
        if(Round == 0 && BarMode == 2 && resfqrd->fb == 0)
        {
            if(DoubleRead->w1 != NULL) addEleToArwlRep(DoubleRead->w1, 12);
            if(DoubleRead->w2 != NULL) addEleToArwlRep(DoubleRead->w2, 13);
            if(DoubleRead->w3 != NULL) addEleToArwlRep(DoubleRead->w3, 14);
            if(DoubleRead->w4 != NULL) addEleToArwlRep(DoubleRead->w4, 15);
            Write = 1;
        }
        else
        {
            switch (resfqrd->type)
            {
                case 'D':
                {
                    FTC->DNA++;
                    if(resfqrd->fb)
                    {
                        FTC->FullDNA++;
                        Write = 1;
                    }
                    else
                    {
                        if(a == 1)
                        {
                            Write = 1;
                        }
                    }
                    if(Write)
                    {
                        if((ReadWrite & 1) > 0) addEleToArwlRep(DoubleRead->w1, 0);
                        if((ReadWrite & 2) > 0) addEleToArwlRep(DoubleRead->w2, 1);
                        if((ReadWrite & 4) > 0) addEleToArwlRep(DoubleRead->w3, 2);
                        if((ReadWrite & 8) > 0) addEleToArwlRep(DoubleRead->w4, 3);
                    }
                    break;
                }
                case 'R':
                {
                    FTC->RNA++;
                    if(resfqrd->fb)
                    {
                        FTC->FullRNA++;
                        Write = 1;
                    }
                    else
                    {
                        if(a == 1)
                        {
                            Write = 1;
                        }
                    }
                    if(Write) 
                    {
                        if((ReadWrite & 1) > 0) addEleToArwlRep(DoubleRead->w1, 4);
                        if((ReadWrite & 2) > 0) addEleToArwlRep(DoubleRead->w2, 5);
                        if((ReadWrite & 4) > 0) addEleToArwlRep(DoubleRead->w3, 6);
                        if((ReadWrite & 8) > 0) addEleToArwlRep(DoubleRead->w4, 7);
                    }
                    break;
                }
                default:
                {
                    FTC->NDNR++;
                    if(resfqrd->fb)
                    {
                        FTC->FullNDNR++;
                        Write = 1;
                    }
                    else
                    {
                        if(a == 1)
                        {   
                            Write = 1;
                        }
                    }
                    if(Write)
                    {
                        if((ReadWrite & 1) > 0) addEleToArwlRep(DoubleRead->w1, 8);
                        if((ReadWrite & 2) > 0) addEleToArwlRep(DoubleRead->w2, 9);
                        if((ReadWrite & 4) > 0) addEleToArwlRep(DoubleRead->w3, 10);
                        if((ReadWrite & 8) > 0) addEleToArwlRep(DoubleRead->w4, 11);
                    }
                    break;
                }
            }
        }
        if(Write)
        {
            if(BarMode == 2 && Round == 0 && resfqrd->fb == 0)
            {
                ReadWrite = TMPReadWrite;
                name->l = 0;
            }

            // printf("Read Write : %d\n", ReadWrite);
            if(ResultPoolR1 && ((ReadWrite & 1) > 0))
            {
                writeToResult(ResultPoolR1, ret, fqrd, name, 1, resfqrd->fb);
            }
            if(ResultPoolR2 && ((ReadWrite & 2) > 0))
            {
                writeToResult(ResultPoolR2, ret, fqrd, name, 2, resfqrd->fb);
            }
            if(ResultPoolR3 && ((ReadWrite & 4) > 0))
            {
                writeToResult(ResultPoolR3, ret, fqrd, name, 3, resfqrd->fb);
            }
            if(ResultPoolR4 && ((ReadWrite & 8) > 0))
            {
                writeToResult(ResultPoolR4, ret, fqrd, name, 4, resfqrd->fb);
            }
        }
        ks_free(name);
        free(name);
    }
    fqseqPoolFree(SeqConfig.fqpool);
    return NULL;
}

arrstr_k * parseIdenFile(char * idfstr)
{
    arrstr_k * Ret = (arrstr_k *) calloc(1, sizeof(arrstr_k));
    arrstr_k * Tmp = splitStrByChar(idfstr, '-');
    for(int i = 0; i < Tmp->l; i++)
    {
        arrstr_k * FileTmp = splitStrByChar(Tmp->mem + Tmp->offset_ar[i], ':');
        if(FileTmp->l != 2) error("-b option error! Should be like [Identifier]:[IDB FILE], for example: CELL:CELL.idb!\n");
        else
        {
            addStrToArrstrRepXC(Ret, FileTmp->mem + FileTmp->offset_ar[0]);
            addStrToArrstrRepXC(Ret, FileTmp->mem + FileTmp->offset_ar[1]);
        }
        arrstrFree(FileTmp);
        free(FileTmp);
    }
    arrstrFree(Tmp);
    free(Tmp);
    return Ret;
}

// int writeOutputToFile(darr_k * DoubleRead, FILE ** OutputFileArr, int BarMode, int Round)
// {
//     if(DoubleRead->w1 != NULL)
//     {
//         // printf("Double Read : %ld - %ld\n", DoubleRead->fb->l, DoubleRead->w1->l);
//         for(int i = 0; i < DoubleRead->Read1->l; i++)
//         {
//             // printf("Write To: %ld\n", DoubleRead->fb->ar[(int) i/4]);
//             if(BarMode == 2 && Round == 0 && DoubleRead->fb->ar[(int) i/4] != 1)
//             {
//                 // printf("Write To: %ld\n", DoubleRead->fb->ar[(int) i/4]);
//                 fprintf(OutputFileArr[6], "%s\n", DoubleRead->Read1->mem + DoubleRead->Read1->offset_ar[i]);
//             }
//             else
//             {
//                 // printf("%s\n", DoubleRead->Read1->mem + DoubleRead->Read1->offset_ar[i]);
//                 fprintf(OutputFileArr[DoubleRead->w1->ar[(int) i/4]], "%s\n", DoubleRead->Read1->mem + DoubleRead->Read1->offset_ar[i]);
//             }
//         }
//     }
//     if(DoubleRead->w2 != NULL)
//     {
//         for(int i = 0; i < DoubleRead->Read2->l; i++)
//         {
//             if(BarMode == 2 && Round == 0 && DoubleRead->fb->ar[(int) i/4] != 1)
//             {
//                 fprintf(OutputFileArr[7], "%s\n", DoubleRead->Read2->mem + DoubleRead->Read2->offset_ar[i]);
//             }
//             else
//             {
//                 fprintf(OutputFileArr[DoubleRead->w2->ar[(int) i/4]], "%s\n", DoubleRead->Read2->mem + DoubleRead->Read2->offset_ar[i]);
//             }
//         }
//     }
//     return 0;
// }

int writeOutputToFile(darr_k * DoubleRead, FILE ** OutputFileArr, int BarMode, int Round)
{
    if(DoubleRead->w1 != NULL)
    {
        // printf("Double Read : %ld - %ld\n", DoubleRead->Read1->l, DoubleRead->Read2->l);
        for(int i = 0; i < DoubleRead->Read1->l; i++)
        {
            fprintf(OutputFileArr[DoubleRead->w1->ar[(int) i/4]], "%s\n", DoubleRead->Read1->mem + DoubleRead->Read1->offset_ar[i]);
        }
    }
    if(DoubleRead->w2 != NULL)
    {
        for(int i = 0; i < DoubleRead->Read2->l; i++)
        {
            fprintf(OutputFileArr[DoubleRead->w2->ar[(int) i/4]], "%s\n", DoubleRead->Read2->mem + DoubleRead->Read2->offset_ar[i]);
        }
    }
    if(DoubleRead->w3 != NULL)
    {
        for(int i = 0; i < DoubleRead->Read3->l; i++)
        {
            fprintf(OutputFileArr[DoubleRead->w3->ar[(int) i/4]], "%s\n", DoubleRead->Read3->mem + DoubleRead->Read3->offset_ar[i]);
        }
    }
    if(DoubleRead->w4 != NULL)
    {
        for(int i = 0; i < DoubleRead->Read4->l; i++)
        {
            fprintf(OutputFileArr[DoubleRead->w4->ar[(int) i/4]], "%s\n", DoubleRead->Read4->mem + DoubleRead->Read4->offset_ar[i]);
        }

    }
    return 0;
}

int parseRead(pthrs_k * ThreadArray, fqhd_k * fqh, int ChunkSize, arwl_k *** ThreadArwlBarstat, barcl_k * barcl, kh_m16_t ** hart, 
                arwl_k * idfcount, int min_qual, int dropN, int BarMode, FILE ** OutputFileArr, int Round, FragTypeStat * ThreadFragTypeStat, int a, 
                darr_k ** ResultStrArr, barc_k ** BarCounterArr, arwl_k * BarCounter, FragTypeStat * TotalTypeStat, int TMPReadWrite)
{
    int MaxThread = 0, ProcessedRead = 0;
    int FinishFlag = 0, RoundCount = 0;
    parseq_k * SeqConfigArr = (parseq_k *) calloc(ThreadArray->l, sizeof(parseq_k));
    for(;;)
    {
        int ThreadCreated = 0;
        for(size_t i = 0; i < ThreadArray->l; i++)
        {
            // printf("i : %ld - state: %ld\n", i, ThreadArray->statear[i]);
            if(ThreadArray->statear[i] == 1)
            {
                
                fqpool_k SeqPool = loadFastq(fqh, ChunkSize);

                if(SeqPool.count > 0)
                {
                    ThreadArwlBarstat[i] = (arwl_k **) calloc(barcl->l, sizeof(arwl_k *));
                    ThreadFragTypeStat[i].DNA = 0; ThreadFragTypeStat[i].RNA = 0; ThreadFragTypeStat[i].NDNR = 0;
                    ThreadFragTypeStat[i].FullDNA = 0; ThreadFragTypeStat[i].FullRNA = 0; ThreadFragTypeStat[i].FullNDNR = 0;
                    SeqConfigArr[i].ftc = &ThreadFragTypeStat[i];
                    SeqConfigArr[i].pread = ProcessedRead;
                    SeqConfigArr[i].torder = i;
                    SeqConfigArr[i].barcl = barcl;
                    SeqConfigArr[i].fqpool = SeqPool;
                    SeqConfigArr[i].har = hart;
                    SeqConfigArr[i].idfcount = idfcount;
                    SeqConfigArr[i].barstat = ThreadArwlBarstat[i];
                    SeqConfigArr[i].fb = a;
                    SeqConfigArr[i].min_qual = min_qual;
                    SeqConfigArr[i].dropN = dropN;
                    SeqConfigArr[i].barmode = BarMode;
                    SeqConfigArr[i].readwrite = TMPReadWrite;
                    // ResultStrArr[i]->w1 = -1;
                    // ResultStrArr[i]->w2 = -1;
                    SeqConfigArr[i].DoubleRead = ResultStrArr[i];
                    SeqConfigArr[i].barcounter = BarCounterArr[i];
                    SeqConfigArr[i].round = Round;
                    pthread_create(&ThreadArray->ptar[i], NULL, parseSeq, &SeqConfigArr[i]);
                    ThreadCreated++;
                    ThreadArray->statear[i] = 0;
                    ProcessedRead = ProcessedRead + SeqPool.count;
                }
                else
                {
                    FinishFlag = 1;
                }
            }
        }
        if(ThreadCreated > MaxThread) MaxThread = ThreadCreated;

        RoundCount++;
        for(size_t i = 0; i < ThreadCreated; i++)
        {
            pthread_join(ThreadArray->ptar[i], NULL);
            ThreadArray->statear[i] = 1;
            for(int j = 0; j < barcl->l; j++) 
            {
                // printf("%ld - %ld\n", ThreadArwlBarstat[i][j]->l, BarCounter[j].l);
                for(int k = 0; k < ThreadArwlBarstat[i][j]->l; k++)
                {
                    // printf("Num : %ld\n", ThreadArwlBarstat[i][j]->ar[k]);
                    BarCounter[j].ar[ThreadArwlBarstat[i][j]->ar[k]]++;
                }
                kaFree(ThreadArwlBarstat[i][j]);
                free(ThreadArwlBarstat[i][j]);
            }
            free(ThreadArwlBarstat[i]);

            if(BarMode == 2 && Round != 1)
            {
                for(int j = 0; j < barcl->l; j++)
                {
                    for(int k = 0; k < barcl->barch[j].l; k++)
                    {
                        for(int m = 0; m < barcl->barch[j].barh[k].l; m++)
                        {
                            if(barcl->barch[j].barh[k].bar[m].hash != NULL) 
                            {
                                for(int p = 0; p < BarCounterArr[i]->l; p++)
                                {
                                    addCountForBarcode(barcl->barch[j].barh[k].bar[m].hash, BarCounterArr[i]->bar[p]);
                                }
                                kbFree(BarCounterArr[i]);
                            }
                        }
                    }
                }
            }
            // Type stat num = 6
            TotalTypeStat->DNA = ThreadFragTypeStat[i].DNA + TotalTypeStat->DNA; 
            TotalTypeStat->RNA = TotalTypeStat->RNA + ThreadFragTypeStat[i].RNA;
            TotalTypeStat->NDNR = TotalTypeStat->NDNR + ThreadFragTypeStat[i].NDNR;
            TotalTypeStat->FullDNA = TotalTypeStat->FullDNA + ThreadFragTypeStat[i].FullDNA;
            TotalTypeStat->FullRNA = TotalTypeStat->FullRNA + ThreadFragTypeStat[i].FullRNA;
            TotalTypeStat->FullNDNR = TotalTypeStat->FullNDNR + ThreadFragTypeStat[i].FullNDNR;
            writeOutputToFile(ResultStrArr[i], OutputFileArr, BarMode, Round);
            arrstrFree(ResultStrArr[i]->Read1);
            arrstrFree(ResultStrArr[i]->Read2);
            arrstrFree(ResultStrArr[i]->Read3);
            arrstrFree(ResultStrArr[i]->Read4);
            kaFree(ResultStrArr[i]->fb);
            kaFree(ResultStrArr[i]->w1);
            kaFree(ResultStrArr[i]->w2);
            kaFree(ResultStrArr[i]->w3);
            kaFree(ResultStrArr[i]->w4);
        }
        processeinfo("%d.", ProcessedRead);
        if(FinishFlag == 1) break;
    }
    free(SeqConfigArr);
    return 0;
}

// convert hash key to barcode string. chain: nth identifier.
int convertKeyToStr(uint64_t key, barcl_k * barcl, kstring_t * str, int chain)
{
    uint64_t VUnit = barcl->barch[0].maxb;
    uint64_t ValueExtract = (1 << VUnit) - 1;
    arrstr_k * StrArr = (arrstr_k *) calloc(1, sizeof(arrstr_k));

    for(int i = 0; i < barcl->l; i++)
    {
        barch_k * barch = &barcl->barch[i];
        for(int j = barch->l - 1; j >= 0; j--)
        {
            uint64_t Value = key & ValueExtract;
            kstring_t TmpStr = {0, 0, NULL};

            if(barch->barh[j].l > 1)
            {
                continue;
            }
            else
            {
                int AddFlag = 0;
                if(barch->barh[j].bar[0].ord != SIZE_MAX && barch->barh[j].bar[0].hash)
                {
                    char * BarName = barch->barh[j].bar[0].rec->mem + barch->barh[j].bar[0].rec->offset_ar[Value];
                    kputs(BarName, &TmpStr);
                    AddFlag = 1;
                    key = key >> VUnit;
                }
                if(AddFlag == 0) addStrToArrstrRepXC(StrArr, "-"); else addStrToArrstrRepXC(StrArr, TmpStr.s);
            }

            ks_free(&TmpStr);

        }
    }

    for(int i = 0; i < StrArr->l; i++)
    {
        if(i > 0) kputc('|', str);
        kputs(StrArr->mem + StrArr->offset_ar[i], str);
    }

    return 1;
}

int writeHashTable(kh_m16_t * har, char * fout, barcl_k * barcl, int chain)
{
    int WriteIdb = 1;
    // if there are random barcode in the chain, won't write idb file
    for(int i = 0; i < barcl->l; i++)
    {
        barch_k * barch = &barcl->barch[i];
        for(int j = 0; j < barch->l; j++)
        {
            for(int k = 0; k < barch->barh[j].l; k++)
            {
                if(barch->barh[j].bar[k].random == 1)
                {
                    WriteIdb = 0;
                    return 0;
                }
            }
        }
    }

    if(WriteIdb == 1)
    {
        FILE * fo = fopen(fout, "w");
        hidv_k * ret;
        khint_t k;
        uint64_t key;
        kstring_t str = {0, 0, NULL};
        for( k = 0; k < kh_end(har); ++k)
        {
            if(kh_exist(har, k))
            {
                key = kh_key(har, k);
                ret = kh_val(har, k);
                convertKeyToStr(key, barcl, &str, chain);
                fprintf(fo, "%lu\t%ld\t%ld\t%s\n", key, ret->idfcount, ret->count, str.s);
                ks_free(&str);
            }
        }
        fclose(fo);
    }
    return 0;
}

int parseBarcode(int argc, char * argv[])
{
    int Nargs, ThreadNum = 1, a = 0, ChunkSize = 10000, dropN = 0;
    char * Prefix = "-", * Config = "-", * IdfStr = "-", * Read1Input = "-", * Read2Input = "-", * Read3Input = "-", * Read4Input = "-";
    // file input handler
    fqhd_k fqh = {0};
    fqh.rd1 = "-";
    fqh.rd2 = "-";
    fqh.rd3 = "-";
    fqh.rd4 = "-";
    // Model: default ; rnat; atact; minq 
    char * Model = "-";
    int BarMode = 0;

    static const struct option lopts[] = {
        {"dropN", no_argument, NULL, 0}, 
        {"min-mistq", required_argument, NULL, 0},
        {"stdout", 0, NULL, 0},
        {"mode", required_argument, NULL, 0},
        {NULL, 0, NULL, 0}
    };
    int min_qual = -1;
    int c;
    int opidx = 0;

    if(argc == 1) 
    {
        Usage();
        exit(1);
    }
    while ((c = getopt_long(argc, argv, "1:2:3:4:c:o:@:ab:m:q:s:", lopts, &opidx)) >= 0)
    {
        if(c == -1) break;

        switch(c)
        {
            case '1' : Read1Input = optarg; break;
            case '2' : Read2Input = optarg; break;
            case '3' : Read3Input = optarg; break;
            case '4' : Read4Input = optarg; break;
            case 'c' : Config = optarg; break;
            case 'o' : Prefix = optarg; break;
            case 'b' : IdfStr = optarg; break;
            case '@' : ThreadNum = atoi(optarg); break;
            case 'a' : a = 1; break;
            case 's' : ChunkSize = atoi(optarg); break;
            // case 'm' : {
            //     char * q;
            //     MaxMem = strtol(optarg, &q, 0);
            //     if( *q == 'k' || *q == 'K') MaxMem <<= 10;
            //     else if(*q == 'm' || *q == 'M') MaxMem <<= 20;
            //     else if(*q == 'g' || *q == 'G') MaxMem <<= 30;
            //     break;
            // }
            case 'q' : min_qual = atoi(optarg); break;
            case 0 : 
            {
                const char * lp = lopts[opidx].name;
                if(strcmp(lp, "dropN") == 0) 
                {
                    dropN = 1; break;
                }
                if(strcmp(lp, "min-miseq") == 0) 
                {
                    min_qual = atoi(optarg); break;
                }
                // todo: output to stdout.
                if(strcmp(lp, "stdout") == 0) 
                {
                    break;
                }
                if(strcmp(lp, "mode") == 0) 
                {
                    Model = optarg; 
                    break;
                }
                break;
            }
            case '?' : Usage(); return -1;
        }
    }

    Nargs = argc - optind;
    if(Nargs > 0 && isatty(STDIN_FILENO))
    {
        warnings("Please specify read 1 and read 2 through option \"-1 Read1\" and \"-2 Read2\" \n");
        Usage();
        exit(1);
    }

    // mode 0 : default : find and stop; 1 : 10x rna seq; 2 10x atac seq; 3 minus quality
    if((strcmp(Model, "default") == 0) || strcmp(Model, "-") == 0 || strcmp(Model, "atact") == 0) 
    {
        BarMode = 2;
    }
    else if(strcmp(Model, "rnat") == 0)
    {
        BarMode = 1;
    }
    else if(strcmp(Model, "fcfg") == 0)
    {
        BarMode = 0;
    }
    else if(strcmp(Model, "minq") == 0)
    {
        BarMode = 3;
    }
    else
    {
        warnings("Barcode identification mode %s is not supported by BARP.\n", Model);
    }

    size_t i, j;
    initMistComArray();
    initPrekQArray();

    if(strcmp(Config, "-") == 0) 
    {
        fprintf(stderr, ANSI_COLOR_RED "[error]  " ANSI_COLOR_RESET ANSI_COLOR_MAGENTA "Config file is necessary!" ANSI_COLOR_RESET" \n");
        return -1;
    }
    
    arrstr_k * Read1List = splitStrByChar(Read1Input, ',');
    arrstr_k * Read2List = splitStrByChar(Read2Input, ',');
    arrstr_k * Read3List = splitStrByChar(Read3Input, ',');
    arrstr_k * Read4List = splitStrByChar(Read4Input, ',');

    size_t FileFlag = 0;
    int FN[4] = {-1, -1, -1, -1};
    size_t MaxFileNum = 0;
    if((strcmp(Read1List->mem + Read1List->offset_ar[0], "-") != 0))
    {
        FN[0] = Read1List->l;
        MaxFileNum = Read1List->l;
        FileFlag = FileFlag | 1;
    }
    if((strcmp(Read2List->mem + Read2List->offset_ar[0], "-") != 0))
    {
        FN[1] = Read2List->l;
        if(Read2List->l > MaxFileNum) MaxFileNum = Read2List->l;
        FileFlag = FileFlag | 2;
    }
    if((strcmp(Read3List->mem + Read3List->offset_ar[0], "-") != 0))
    {
        FN[2] = Read3List->l;
        if(Read3List->l > MaxFileNum) MaxFileNum = Read3List->l;
        FileFlag = FileFlag | 4;
    }    
    if((strcmp(Read4List->mem + Read4List->offset_ar[0], "-") != 0))
    {
        FN[3] = Read4List->l;
        if(Read4List->l > MaxFileNum) MaxFileNum = Read4List->l;
        FileFlag = FileFlag | 8;
    }
    
    fqh.flag = FileFlag;

    for(int i = 0; i < 4; i++)
    {
        if(FN[i] != -1 && FN[i] != MaxFileNum) error("File number dose not match.");
    }


    fqh.Read1List = Read1List;
    fqh.Read2List = Read2List;
    fqh.Read3List = Read3List;
    fqh.Read4List = Read4List;    

    if(openFastq(&fqh)) error("Open FASTQ file failed.");

    barcl_k * barcl = parseConfig(Config);
    int BigHash = 0;
    for(i = 0; i < barcl->l; i++)
    {
        fprintf(stdout, ANSI_COLOR_GREEN "[Barcode Chain] :" ANSI_COLOR_CYAN);
        for(j = 0; j < barcl->barch[i].l; j++)
        {
            if(j == 0) fprintf(stdout, " %s", barcl->barch[i].barh[j].bar[0].na); else fprintf(stdout, " - %s", barcl->barch[i].barh[j].bar[0].na);
            if(j == 0) 
            {
                if(barcl->barch[i].barh[j].bar[0].bighash)
                {
                    BigHash = 1;
                }
            }
            for(int m = 1; m < barcl->barch[i].barh[j].l; m++)
            {
                if(barcl->barch[i].barh[j].bar[m].bighash)
                {
                    BigHash = 1;
                }
                fprintf(stdout, " or %s", barcl->barch[i].barh[j].bar[m].na);
            }
        }
        fprintf(stdout, ANSI_COLOR_RESET "\n");
    }
    if(BigHash == 0) BarMode = 0;
    kh_m16_t ** hart = (kh_m16_t **) calloc(barcl->idl, sizeof(kh_m16_t *));

    // for(i = 0; i < barcl->jump->l; i++)
    // {
    //     printf("Jump[%ld]: %ld\n",i, barcl->jump->ar[i]);
    // }

    for(int i = 0; i < barcl->idl; i++)
    {
        hart[i] = kh_init(m16);
    }

    arwl_k * idfcount;
    if(barcl->idf)
    {
        idfcount = (arwl_k *) calloc(1, sizeof(arwl_k));
        idfcount->l = barcl->idl;
        idfcount->ar = (size_t *) calloc(1, barcl->idl * sizeof(size_t));
    }
    // add identifier from file : CELL:CELL.idb-COMPLEX:COMPLEX.idb
    if(strcmp(IdfStr, "-") != 0)
    {
        fprintf(stdout, "IdfStr: %s\n", IdfStr);
        arrstr_k * IdenFileStrArray = parseIdenFile(IdfStr);
        for(int i = 0; i < barcl->idl; i++)
        {
            for(int j = 0; j < IdenFileStrArray->l; j++)
            {
                if(strcmp(barcl->idf[i].na, IdenFileStrArray->mem + IdenFileStrArray->offset_ar[j]) == 0)
                {
                    fprintf(stdout, ANSI_COLOR_CYAN "Read identifier %s from %s.", barcl->idf[i].na, IdenFileStrArray->mem + IdenFileStrArray->offset_ar[j + 1]);
                    fprintf(stdout, ANSI_COLOR_RESET "\n");
                    if(!(idfcount->ar[i] = pushIdentifierFromFile(hart[i], IdenFileStrArray->mem + IdenFileStrArray->offset_ar[j + 1]))) error("Read Identifier from file corrupted\n");
                }
            }
        }
        arrstrFree(IdenFileStrArray);
        free(IdenFileStrArray);
    }

    if(ThreadNum > 1) initLocker();

    FragTypeStat TotalTypeStat = {0};
    char * OutputFileNamePre[OUTPUTTPYENUM] = {"DNA_1", "DNA_2","DNA_3", "DNA_4", "RNA_1", "RNA_2","RNA_3", "RNA_4", "NDNR_1", "NDNR_2","NDNR_3", "NDNR_4", "ATAC_1_TMP", "ATAC_2_TMP", "ATAC_3_TMP", "ATAC_4_TMP"};
    kstring_t * OutputFileName = (kstring_t *) calloc(OUTPUTTPYENUM, sizeof(kstring_t));
    for(int i = 0; i < OUTPUTTPYENUM; i++) 
    {
        ks_initialize(&OutputFileName[i]);
        if(strcmp(Prefix, "-") != 0) 
        {
            kputs(Prefix, &OutputFileName[i]);
            kputs("_", &OutputFileName[i]);
        }
        kputs(OutputFileNamePre[i], &OutputFileName[i]);
        if(a == 1)
        {
            kputs("_all.fq", &OutputFileName[i]);
        }
        else
        {
            kputs(".fq", &OutputFileName[i]);
        }
    }

    FILE * OutputFileArr[OUTPUTTPYENUM];
    for(int j = 0; j < OUTPUTTPYENUM; j++)
    {
        OutputFileArr[j] = fopen(OutputFileName[j].s, "w");
    }

    pthrs_k ThreadArray;
    initThread(&ThreadArray, ThreadNum);
    arwl_k *** ThreadArwlBarstat = (arwl_k ***) calloc(ThreadNum, sizeof(arwl_k **));
    // statistic for barcode chain
    arwl_k * BarCounter = (arwl_k *) calloc(barcl->l, sizeof(arwl_k));
    FragTypeStat * ThreadFragTypeStat = (FragTypeStat *) calloc(ThreadNum, sizeof(FragTypeStat));
    // barcount count barcode chain situation for each barcode chain.
    for(int i = 0; i < barcl->l; i++)
    {
        barch_k * curbarch = &barcl->barch[i];
        BarCounter[i].l = (1 << curbarch->l);
        BarCounter[i].ar = (size_t *) calloc(BarCounter[i].l, sizeof(size_t));
    }
    FILE * BarStatFile = fopen("BarcodeIdentification_statistic.txt", "w");
    if(BarStatFile == NULL)
    {
        warnings("Failed to open BarcodeIdentification_statistic.txt\n");
    }
    darr_k ** ResultStrArr = (darr_k **) calloc(ThreadNum, sizeof(darr_k *));
    for(int i = 0; i < ThreadNum; i++)
    {
        ResultStrArr[i] = (darr_k *) calloc(1, sizeof(darr_k));
        ResultStrArr[i]->Read1 = (arrstr_k *) calloc(1, sizeof(arrstr_k));
        ResultStrArr[i]->w1 = (arwl_k *) calloc(1, sizeof(arwl_k));
        ResultStrArr[i]->Read2 = (arrstr_k *) calloc(1, sizeof(arrstr_k));
        ResultStrArr[i]->w2 = (arwl_k *) calloc(1, sizeof(arwl_k));
        ResultStrArr[i]->Read3 = (arrstr_k *) calloc(1, sizeof(arrstr_k));
        ResultStrArr[i]->w3 = (arwl_k *) calloc(1, sizeof(arwl_k));
        ResultStrArr[i]->Read4 = (arrstr_k *) calloc(1, sizeof(arrstr_k));
        ResultStrArr[i]->w4 = (arwl_k *) calloc(1, sizeof(arwl_k));
        ResultStrArr[i]->fb = (arwl_k *) calloc(1, sizeof(arwl_k));
    }

    barc_k ** BarCounterArr = (barc_k **) calloc(ThreadNum, sizeof(barc_k *));
    for(int i = 0; i < ThreadNum; i++)
    {
        BarCounterArr[i] = (barc_k *) calloc(1, sizeof(barc_k));
    }
    // first round of barcode identification
    parseRead(&ThreadArray, &fqh, ChunkSize, ThreadArwlBarstat, barcl, hart, idfcount, min_qual, dropN, BarMode, OutputFileArr, 0, ThreadFragTypeStat, a, ResultStrArr, BarCounterArr, BarCounter, &TotalTypeStat, FileFlag);
    // write exact barcode count
    // if(BarMode == 2) writeBarcodeHash(barcl->barch[0].barh[0].bar[0].hash, "BC-1");

    if(BarMode == 2)
    {
        gzclose(fqh.r1);
        gzclose(fqh.r2);
        gzclose(fqh.r3);
        gzclose(fqh.r4);
        fclose(OutputFileArr[12]);
        fclose(OutputFileArr[13]);
        fclose(OutputFileArr[14]);
        fclose(OutputFileArr[15]);
        fqh.r1 = gzopen(OutputFileName[12].s, "r");
        fqh.r2 = gzopen(OutputFileName[13].s, "r");
        fqh.r3 = gzopen(OutputFileName[14].s, "r");
        fqh.r4 = gzopen(OutputFileName[15].s, "r");
        if(fqh.r1) 
        {
            if(fqh.ks1) kseq_destroy(fqh.ks1);
            if((FileFlag & 1) > 0) fqh.ks1 = kseq_init(fqh.r1); 
        }
        else 
        {
            fprintf(stderr,ANSI_COLOR_RED "[error]  " ANSI_COLOR_RESET ANSI_COLOR_MAGENTA  "No such file or file error: %s."ANSI_COLOR_RESET "\n", fqh.rd1);
            exit(1);
        }

        if(fqh.r2)
        {
            if(fqh.ks2) kseq_destroy(fqh.ks2);
            if((FileFlag & 2) > 0)fqh.ks2 = kseq_init(fqh.r2); 
        }
        else 
        {
            fprintf(stderr,ANSI_COLOR_RED "[error]  " ANSI_COLOR_RESET ANSI_COLOR_MAGENTA  "No such file or file error: %s."ANSI_COLOR_RESET "\n", fqh.rd2);
            exit(1);
        }
        if(fqh.r3)
        {
            if(fqh.ks3) kseq_destroy(fqh.ks3);
            if((FileFlag & 4) > 0)fqh.ks3 = kseq_init(fqh.r3); 
        }
        else 
        {
            fprintf(stderr,ANSI_COLOR_RED "[error]  " ANSI_COLOR_RESET ANSI_COLOR_MAGENTA  "No such file or file error: %s."ANSI_COLOR_RESET "\n", fqh.rd3);
            exit(1);
        }
        if(fqh.r4)
        {
            if(fqh.ks4) kseq_destroy(fqh.ks4);
            if((FileFlag & 8) > 0) fqh.ks4 = kseq_init(fqh.r4); 
        }
        else 
        {
            fprintf(stderr,ANSI_COLOR_RED "[error]  " ANSI_COLOR_RESET ANSI_COLOR_MAGENTA  "No such file or file error: %s."ANSI_COLOR_RESET "\n", fqh.rd4);
            exit(1);
        }

        parseRead(&ThreadArray, &fqh, ChunkSize, ThreadArwlBarstat, barcl, hart, idfcount, min_qual, dropN, BarMode, OutputFileArr, 1, ThreadFragTypeStat, a, ResultStrArr, BarCounterArr, BarCounter,  &TotalTypeStat, FileFlag);
    }

    if(ThreadNum > 1)
    {
        destroyLocker();
    }

    int fields = 0;
    for(int i = 0; i < OUTPUTTPYENUM - 4; i++)
    {
        fclose(OutputFileArr[i]);
    }

    for(int j = 0; j < OUTPUTTPYENUM - 4; j++)
    {
        if((fields = open(OutputFileName[j].s, O_RDONLY)))
        {
            struct stat buf;
            fstat(fields, &buf);
            if( buf.st_size == 0)
            {
                close(fields);
                remove(OutputFileName[j].s);
            }
        }
    }
    for(int j = OUTPUTTPYENUM - 4; j < OUTPUTTPYENUM; j++)
    {
        remove(OutputFileName[j].s);
    }
    char NotFound[10] = "Not_Found";

    if(BarStatFile)
    {
        if(barcl->idf)
        {
            for(int j = 0; j < idfcount->l; j++)
            {
                if((strcmp(barcl->idf[j].na, "DNA") != 0) && (strcmp(barcl->idf[j].na, "RNA") != 0))
                {
                    fprintf(BarStatFile, "%s : %ld\n", barcl->idf[j].na ,idfcount->ar[j]);
                }
            }
        }
        for(int i = 0; i < barcl->l; i++)
        {
            barch_k * curbarch = &barcl->barch[i];
            arrstr_k * ResultAtrStatArray = (arrstr_k *) calloc(BarCounter[i].l, sizeof(arrstr_k));
            for(int k = 0; k < BarCounter[i].l; k++)
            {
                int QFlag = 1 << (curbarch->l - 1);
                for(int j = 0; j < curbarch->l; j++)
                {
                    char * CurName;
                    if((k & QFlag) != QFlag) 
                    {
                        CurName = NotFound;
                    }
                    else
                    {
                        char * PreName = strConcate(curbarch->barh[j].bar[0].na, "", "");
                        for(int m = 1; m < curbarch->barh[j].l; m++)
                        {
                            char * TmpName = strConcate(PreName, curbarch->barh[j].bar[m].na, "/");
                            free(PreName);
                            PreName = TmpName;
                        }
                        CurName = PreName;
                    } 
                    addStrToArrstrRepXC(&ResultAtrStatArray[k], CurName);
                    QFlag = QFlag >> 1;
                }
            }
            size_t * SortedResultStat = selectSortOrder(&BarCounter[i], 1);
            for(int k = 0; k < BarCounter[i].l; k++)
            {
                for(int j = 0; j < ResultAtrStatArray[k].l; j++)
                {
                    if(j == 0) 
                    {
                        fprintf(BarStatFile,"[%s]", ResultAtrStatArray[SortedResultStat[k]].mem + ResultAtrStatArray[SortedResultStat[k]].offset_ar[j]);
                    }
                    else
                    {
                        fprintf(BarStatFile, " [%s]", ResultAtrStatArray[SortedResultStat[k]].mem + ResultAtrStatArray[SortedResultStat[k]].offset_ar[j]);
                    } 

                }
                fprintf(BarStatFile, " | Read Count : %ld\n", BarCounter[i].ar[SortedResultStat[k]]);
            }
            free(SortedResultStat);
            for(int f = 0; f < BarCounter[i].l; f++) arrstrFree(&ResultAtrStatArray[f]);
            free(ResultAtrStatArray);
        }
        fclose(BarStatFile);
    }

    for(int i = 0; i < barcl->l; i++) kaFree(&BarCounter[i]);
    free(BarCounter);
    char * StatFile;
    if(strcmp(Prefix, "-") == 0) StatFile = "stat"; else StatFile = strConcateByChar(Prefix, "stat", '.');
    FILE * fps = fopen(StatFile, "w");
    fprintf(fps, "Total reads:%ld\nTotal DNA reads:%ld\nFully barcoded DNA reads:%ld\n"
    "Not fully barcoded DNA reads:%ld\nTotal RNA reads:%ld\nFully barcoded RNA reads:%ld\n"
    "Not fully barcoded RNA reads:%ld\nTotal non-defined reads:%ld\n"
    "Fully barcoded non-defined reads:%ld\nNot fully barcoded non-defined reads:%ld\n"
    , (TotalTypeStat.DNA + TotalTypeStat.RNA + TotalTypeStat.NDNR), TotalTypeStat.DNA, TotalTypeStat.FullDNA, 
    (TotalTypeStat.DNA-TotalTypeStat.FullDNA), TotalTypeStat.RNA, TotalTypeStat.FullRNA, 
    (TotalTypeStat.RNA-TotalTypeStat.FullRNA), TotalTypeStat.NDNR, TotalTypeStat.FullNDNR, 
    (TotalTypeStat.NDNR-TotalTypeStat.FullNDNR));
    fclose(fps);
    if(strcmp(Prefix, "-") != 0) free(StatFile);
    for(int i = 0; i < barcl->idl; i++)
    {
        if(strcmp(barcl->idf[i].na, "DNA") == 0 || strcmp(barcl->idf[i].na, "RNA") == 0)
        {
            continue;
        }
        else
        {
            kstring_t HashFile = {0, 0, NULL};
            if(strcmp(Prefix, "-") != 0)
            {
                kputs(Prefix, &HashFile);
                kputs(".", &HashFile);
            }
            kputs(barcl->idf[i].na, &HashFile);
            kputs(".idb", &HashFile);
            writeHashTable(hart[i], HashFile.s, barcl, i); 
            free(HashFile.s);
        }
    }
    if(fqh.ks1) kseq_destroy(fqh.ks1);
    if(fqh.ks2) kseq_destroy(fqh.ks2);
    if(fqh.ks3) kseq_destroy(fqh.ks3);
    if(fqh.ks4) kseq_destroy(fqh.ks4);
    gzclose(fqh.r1);
    gzclose(fqh.r2);
    gzclose(fqh.r3);
    gzclose(fqh.r4);
    barclFree(barcl);
    free(barcl);
    destroyMistComArray();
    destroyPrekQArray();
    return 0;
}
