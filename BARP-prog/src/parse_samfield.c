#include <stdbool.h>
#include <stdlib.h>
#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <sys/stat.h>
#include <unistd.h>
#include <getopt.h>
#include <pthread.h>
#include "hts.h"
#include "khash.h"
#include "sam.h"
#include "kutils.h"
#include "sam_opts.h"
#include "parse_samfield.h"

#define BARPVERSION "0.1.1"

typedef struct kvpair
{
    char * key;
    char * val;
} kvpa_k;

typedef struct bam1_pool{
    size_t count;
    int mem_full;
    bam1_t ** bam_record;
    uint8_t * bam_mem;
} bam1_pool;


const size_t SORT_MIN_MEGS_PER_THREAD = 1;
const size_t SORT_DEFAULT_MEGS_PER_THREAD = 768;


static void complain_about_memory_setting(size_t max_mem) {
    char  *suffix = "";
    const size_t nine_k = 9<<10;
    if (max_mem > nine_k) { max_mem >>= 10; suffix = "K"; }
    if (max_mem > nine_k) { max_mem >>= 10; suffix = "M"; }

    fprintf(stderr,
"[bam_sort] -m setting (%zu%s bytes) is less than the minimum required (%zuM).\n\n"
"Trying to run with -m too small can lead to the creation of a very large number\n"
"of temporary files.  This may make sort fail due to it exceeding limits on the\n"
"number of files it can have open at the same time.\n\n"
"Please check your -m parameter.  It should be an integer followed by one of the\n"
"letters K (for kilobytes), M (megabytes) or G (gigabytes).  You should ensure it\n"
"is at least the minimum above, and much higher if you are sorting a large file.\n",
            max_mem, suffix, SORT_MIN_MEGS_PER_THREAD);
}

static inline bam1_pool loadSam(samFile * fp, sam_hdr_t * header, size_t max_mem)
{
    bam1_pool ret = {0};
    ret.bam_record = (bam1_t **) calloc(1, sizeof(bam1_t *));
    uint8_t * bam_mem;
    ret.bam_mem = malloc(max_mem);
    bam_mem = ret.bam_mem;
    bam1_t * b = bam_init1();
    if(bam_mem == NULL) error("Could not allocate memory for sam pool.\n");
    int res;
    size_t count = 0, bam_mem_offset = 0, max_k = 0;
    while((bam_mem_offset + sizeof(*b) + b->l_data < max_mem) && ((res = sam_read1(fp, header, b)) >= 0))
    {
        if(count == max_k)
        {
            bam1_t ** newbuf;
            max_k = max_k? max_k<<1 : 0x10000;
            newbuf = realloc(ret.bam_record, max_k * sizeof(bam1_t *));
            if(!newbuf) error("Wronged at allocate memory\n");
            ret.bam_record = newbuf;
        }
        ret.bam_record[count] = (bam1_t *) (bam_mem + bam_mem_offset);
        * ret.bam_record[count]  = * b;
        //             .--------------------.
        //    |   | *data  |       |                   |
        //    {       bam1_t       }{     data         }
        ret.bam_record[count]->data = (uint8_t *) ((char *) ret.bam_record[count] + sizeof(bam1_t));
        memcpy(ret.bam_record[count]->data, b->data, b->l_data);
        // store next BAM record in next 8-byte-aligned address after current one
        bam_mem_offset = (bam_mem_offset + sizeof(*b) + b->l_data + 8 - 1) & ~((size_t) 8 - 1);
        count++;
    }
    if(bam_mem_offset + sizeof(*b) + b->l_data < max_mem) ret.mem_full = 1;
    bam_destroy1(b);
    if(count < 1) error("No bam record found in input file\n");
    ret.count = count;
    return ret;
}

typedef struct m2value
{
    arrstr_k * add;
    optar_k * opts;
    arrstr_k * fathers;
    char * name;
}m2v_k;

typedef struct ArrstrWithFlag
{
    arrstr_k * arr;
    arwl_k * flag;
}arrsf_k;

KHASH_MAP_INIT_INT64(m2, m2v_k *)
typedef struct parsename
{
    bam1_pool bampool;
    size_t index;
    arrsf_k * FieldArray;
    arrstr_k * SonArray;
    arrsf_k * AddFieldArray;
    sam_hdr_t * header;
    kh_m2_t * har;
    samFile * BamOut;
    int rr;
    int writebam;
    int maxlen;
    int strict;
    int atac;
    int extpos;
    int extneg;
}parnam_k;

// a tmp struct to store identifier information
typedef struct IdentifierStruct
{
    char name[128];
}idfs_k;

pthread_mutex_t slock;

int addOptToOpts(optar_k * opts, opts_k opt)
{
    opts_k * tmp = (opts_k *) realloc(opts->opts, (opts->l + 1) * sizeof(opts_k));
    opts->opts = tmp;
    memcpy(opts->opts + opts->l, &opt, sizeof(opts_k));
    opts->l++;
    return 0;
}

m2v_k * queryLastSon(kh_m2_t * har, size_t key, opts_k opt, char * name)
{
    int absent;
    m2v_k * ret;
    khint_t k = kh_get(m2, har, key);
    if(k != kh_end(har))
    {
        m2v_k * m2k = kh_value(har, k);
        addOptToOpts(m2k->opts, opt);
        kh_val(har, k) = m2k;
        ret = m2k;
    }
    else
    {
        k = kh_put(m2, har, key, &absent);
        optar_k * optars = (optar_k *) calloc(1, sizeof(optar_k));
        optars->l = 1;
        optars->opts = calloc(1, sizeof(opts_k));
        memcpy(optars->opts, &opt, sizeof(opts_k));
        m2v_k * m2kv = (m2v_k *) calloc(1, sizeof(m2v_k));
        m2kv->opts = optars;
        kh_val(har, k) = m2kv;
        m2kv->fathers = (arrstr_k *) calloc(1, sizeof(arrstr_k));
        m2kv->name = strCopy(name);
        m2kv->add = (arrstr_k *) calloc(1, sizeof(arrstr_k));
        ret = m2kv;
    }
    return ret;
}

    // if(curbam->core.isize >= -600 && curbam->core.isize <= 600 && ((curbam->core.flag & 0x01) > 0))
    // {
    //     // samfile annotation 14
    //     opt.start = curbam->core.pos;
    //     opt.end = curbam->core.pos + abs(curbam->core.isize);
    // }
    // else
    // {
    //     opt.start = curbam->core.pos;
    //     opt.end = MaxUint64(curbam->core.pos, curbam->core.pos + bam_cigar2qlen(curbam->core.n_cigar, bam_get_cigar(curbam)));
    // }
    // opt.chr = (char *) sam_hdr_tid2name(header, curbam->core.tid);
    // if(opt.chr != NULL)
    // {
    //     m2v_k * Result = queryLastSon(har, CurIdenOrder, opt, IdenArray->mem + IdenArray->offset_ar[i]);
    //     FatherArray = Result->fathers;
    //     AddArray = Result->add;
    // }

// take : 0: default,  1: read 1; 2: read 2;
opts_k getBamPosition(bam1_t ** bampair, int pairflag, size_t maxlen)
{  
    opts_k opt = {0, 0, -1};
    uint64_t s1, s2, l1, l2;
    int32_t c1, c2;
    if(pairflag)
    {
        c1 = bampair[0]->core.tid; c2 = bampair[1]->core.tid;
        if(c1 == -1 && c2 != -1)
        {
            opt.start = bampair[1]->core.pos;
            opt.end = bampair[1]->core.pos + bam_cigar2rlen(bampair[1]->core.n_cigar, bam_get_cigar(bampair[1]));
            opt.tid = c2;
        }
        else if(c1 != -1 && c2 == -1)
        {
            opt.start = bampair[0]->core.pos;
            opt.end = bampair[0]->core.pos + bam_cigar2rlen(bampair[0]->core.n_cigar, bam_get_cigar(bampair[0]));
            opt.tid = c1;
        }
        else
        {
            s1 = bampair[0]->core.pos;
            s2 = bampair[1]->core.pos;
            l1 = bam_cigar2rlen(bampair[0]->core.n_cigar, bam_get_cigar(bampair[0]));
            l2 = bam_cigar2rlen(bampair[1]->core.n_cigar, bam_get_cigar(bampair[1]));
            if(c1 == c2)
            {
                // if(bampair[0]->core.isize > 0)
                // {
                //     if(s2 + l2 - s1 <= maxlen) 
                //     {
                //         opt.start = s1;
                //         opt.end = s2 + l2;
                //         opt.tid = c1;
                //         return opt;
                //     }
                // }
                // else
                // {
                //     if(s1 + l1 - s2 <= maxlen)
                //     {
                //         opt.start = s2;
                //         opt.end = s1 + l1;
                //         opt.tid = c2;
                //         return opt;
                //     }
                // }
                size_t e1 = s1 + l1;
                size_t e2 = s2 + l2;
                opt.end = Max(e1, e2);
                opt.start = Min(s1, s2);
                opt.tid = c1;
                if(opt.end - opt.start <= maxlen)
                {
                    return opt;
                }
            }
            else
            {
                if(l1 >= l2)
                {
                    opt.start = s1;
                    opt.end = s1 + l1;
                    opt.tid = c1;
                }
                else
                {
                    opt.start = s2;
                    opt.end = s2 + l2;
                    opt.tid = c2;
                }
            }
            return opt;
        }
    }
    else
    {
        opt.start = bampair[0]->core.pos + 1;
        opt.end = bampair[0]->core.pos + 1 + bam_cigar2rlen(bampair[0]->core.n_cigar, bam_get_cigar(bampair[0]));
        opt.tid = bampair[0]->core.tid;
    }
    return opt;
}  

void * parseSam(void * arg)
{
    parnam_k SamConfig = * ((parnam_k *) arg);
    arrsf_k * FieldArray = SamConfig.FieldArray;
    arrstr_k * SonArray = SamConfig.SonArray;
    arrsf_k * AddFieldArray = SamConfig.AddFieldArray;
    bam1_pool bampool = SamConfig.bampool;
    kh_m2_t * har = SamConfig.har;
    sam_hdr_t * header = SamConfig.header;
    int writebam = SamConfig.writebam;
    samFile * fp = SamConfig.BamOut;
    int maxlen = SamConfig.maxlen;
    int strict = SamConfig.strict;
    int atac = SamConfig.atac;
    int extpos = SamConfig.extpos;
    int extneg = SamConfig.extneg;
    // 0 : based on read 1 and read 2; 1: based on read 1, 2 : based on read 2; rr : retain read.
    int rr = SamConfig.rr;

    int * IdenFindFlag;
    if(FieldArray)
    {
        IdenFindFlag = (int *) calloc(FieldArray->arr->l, sizeof(int));
    }
    idfs_k * WaitQueue = NULL;
    if(SonArray)
    {
        WaitQueue = (idfs_k *) calloc(SonArray->l, sizeof(idfs_k));
    }
    char TagName[2] = {0};
    opts_k opt = {0, 0, -1};
    bam1_t * bampair[2];
    for(int i = 0; i < 2; i++) 
    {
        bampair[i] = bam_init1();
    }

    int PairFlag = 0;
    char * NextName;
    for(int p = 0; p < bampool.count; p++)
    {
        bam1_t * oribam = bampool.bam_record[p];
        char * SeqName = bam_get_qname(oribam);
        // printf("SeqName : %s\n", SeqName);
        if((oribam->core.flag & 0x01) > 0)
        {
            if(rr == 0)
            {
                bampair[0] = bam_copy1(bampair[0], oribam);
                if((p + 1) < bampool.count)
                {
                    NextName = bam_get_qname(bampool.bam_record[p + 1]);
                    if(strcmp(SeqName, NextName) == 0)
                    {
                        PairFlag = 1;
                        bampair[1] = bam_copy1(bampair[1], bampool.bam_record[++p]);
                        goto ProcessBam;
                    }
                    else
                    {
                        if(strict)
                        {
                            warnings("%s is marked as paired, but its mate does not occur next to it in you BAM file. Will be processed as a not qualified fragment.", SeqName);
                            continue;
                        }
                        else
                        {
                            warnings("%s is marked as paired, but its mate does not occur next to it in you BAM file. Will be processed as a individual fragment.", SeqName);
                        }
                    }
                }
            }
            else if(rr == 1)
            {
                if((oribam->core.flag & 0x40) > 0)
                {
                    bampair[0] = bam_copy1(bampair[0], oribam);
                    PairFlag = 0;
                    goto ProcessBam;
                }
                else
                {
                    continue;
                }
            }
            else if(rr == 2)
            {
                if((oribam->core.flag & 0x80) > 0)
                {
                    bampair[0] = bam_copy1(bampair[0], oribam);
                    PairFlag = 0;
                    goto ProcessBam;
                }
                else
                {
                    continue;
                }
            }
        }
        else
        {
            bampair[0] = bam_copy1(bampair[0], oribam);
        }

        ProcessBam:
        {
            // split name to [Name][Barcode][Iden]
            arrstr_k * AllInfoArray = splitStrByStr(SeqName, "|||");
            arrstr_k * FatherArr, * AddArr;
            int AddFlag = 0;
            if(FieldArray != NULL)
            {
                if(AllInfoArray->l < 3)
                {
                    warnings("%s do not have identifier in its name. Skipping...", SeqName);
                }
                else
                {
                    memset(IdenFindFlag, 0, FieldArray->arr->l);
                    int j, FatherCount = 0;
                    // cell_1|complex_1
                    char * IdenName = AllInfoArray->mem + AllInfoArray->offset_ar[2];
                    arrstr_k * IdenArray = splitStrByChar(IdenName, '|');
                    for(int i = 0; i < IdenArray->l; i++)
                    {
                        arrstr_k * CurIden = splitStrByChar(IdenArray->mem + IdenArray->offset_ar[i], '_');
                        // 1 2 2 1 2 2
                        char * CurIdenName = CurIden->mem + CurIden->offset_ar[0];
                        // parse to bam field
                        size_t CurIdenOrder = atoi(CurIden->mem + CurIden->offset_ar[1]);
                        for(j = 0; j < FieldArray->arr->l; j++)
                        {
                            if(strcmp(CurIdenName, FieldArray->arr->mem + FieldArray->arr->offset_ar[j]) == 0)
                            {
                                IdenFindFlag[j] = 1;
                                // printf("%ld\n", FieldArray->flag->ar[j]);
                                // FieldArray->flag->ar[j] can only be 1 or 2
                                if(FieldArray->flag->ar[j] == 2 || writebam > 0)
                                {
                                    IdenFindFlag[j + 1] = 1;
                                    strncpy(TagName, FieldArray->arr->mem + FieldArray->arr->offset_ar[j + 1], 2);
                                    bam_aux_update_int(bampair[0], TagName, CurIdenOrder);
                                    if(PairFlag)
                                    {
                                        bam_aux_update_int(bampair[1], TagName, CurIdenOrder);
                                    }
                                    writebam = 1;
                                    j++;
                                }
                                break;
                            }
                        }
                        opt = getBamPosition(bampair, PairFlag, maxlen);
                        if(atac)
                        {
                            opt.start = opt.start + extpos;
                            opt.end = opt.end + extneg;
                        }
                        if(SonArray)
                        {
                            for(j = 0; j < SonArray->l; j++)
                            {
                                if(strcmp(CurIdenName, SonArray->mem + SonArray->offset_ar[j]) == 0)
                                {
                                    // not the smallest identifier
                                    if(j != SonArray->l - 1) 
                                    {
                                        strcpy(WaitQueue[j].name, IdenArray->mem + IdenArray->offset_ar[i]);
                                        WaitQueue[j].name[strlen(IdenArray->mem + IdenArray->offset_ar[i])] = '\0';
                                        FatherCount++;
                                    }
                                    else
                                    {
                                        if(opt.tid != -1)
                                        {
                                            m2v_k * Result = queryLastSon(har, CurIdenOrder, opt, IdenArray->mem + IdenArray->offset_ar[i]);
                                            FatherArr = Result->fathers;
                                            AddArr = Result->add;
                                            AddFlag = 1;
                                        }
                                        if(FatherArr)
                                        {
                                            if(FatherArr->l == 0)
                                            {
                                                for(j = 0; j < FatherCount; j++)
                                                {
                                                    addStrToArrstrRepXC(FatherArr, WaitQueue[j].name);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        else
                        {
                            m2v_k * Result = queryLastSon(har, CurIdenOrder, opt, IdenArray->mem + IdenArray->offset_ar[i]);
                            AddArr = Result->add;
                            AddFlag = 1;
                        }
                        arrstrFree(CurIden);
                        free(CurIden);
                    }
                    for(j = 0; j < FieldArray->arr->l; j++)
                    {
                        if(IdenFindFlag[j] == 0) warnings("Field \"%s\" not found at %s", FieldArray->arr->mem + FieldArray->arr->offset_ar[j], SeqName);
                    }
                    arrstrFree(IdenArray);
                    free(IdenArray);
                }
            }

            if(AddFieldArray != NULL)
            {
                char * BarChainName = AllInfoArray->mem + AllInfoArray->offset_ar[1];
                arrstr_k * BarChain = splitStrByChar(BarChainName, '|');
                if(AddFieldArray->flag->ar[AddFieldArray->flag->l - 1] == 1)
                {
                    // loop through all barcode
                    for(int i = 0; i < BarChain->l; i++)
                    {
                        arrstr_k * TagContent = splitStrByChar(BarChain->mem + BarChain->offset_ar[i], ':');
                        char * CodeName = TagContent->mem + TagContent->offset_ar[0];
                        char * CodeContent = TagContent->mem + TagContent->offset_ar[1];
                        for(int j = 0; j < AddFieldArray->arr->l; j++)
                        {
                            switch(AddFieldArray->flag->ar[j])
                            {
                                case 1:
                                case 3:
                                {
                                    if(strcmp(AddFieldArray->arr->mem + AddFieldArray->arr->offset_ar[j], CodeName) == 0)
                                    {
                                        strncpy(TagName, AddFieldArray->arr->mem + AddFieldArray->arr->offset_ar[j + 1], 2);
                                        bam_aux_update_str(bampair[0], TagName, strlen(CodeContent) + 1, CodeContent);
                                        if(PairFlag)
                                        {
                                            bam_aux_update_str(bampair[1], TagName, strlen(CodeContent) + 1, CodeContent);
                                        }
                                        j++;
                                        if(AddFlag)
                                        {
                                            addStrToArrstrRepXC(AddArr, CodeContent);
                                        }
                                        break;
                                    }
                                }
                                default:
                                    break;
                            }
                        }
                        arrstrFree(TagContent);
                        free(TagContent);
                    }
                }
                else
                {
                    for(int i = 0; i < AddFieldArray->arr->l; i++)
                    {
                        switch(AddFieldArray->flag->ar[i])
                        {
                            case 0:
                            {
                                int BarcodeOrder = atoi(AddFieldArray->arr->mem + AddFieldArray->arr->offset_ar[i]) - 1;
                                arrstr_k * TagContent = splitStrByChar(BarChain->mem + BarChain->offset_ar[BarcodeOrder], ':');
                                strncpy(TagName, AddFieldArray->arr->mem + AddFieldArray->arr->offset_ar[i + 1], 2);
                                bam_aux_update_str(bampair[0], TagName, strlen(TagContent->mem + TagContent->offset_ar[1]) + 1, TagContent->mem + TagContent->offset_ar[1]);
                                if(PairFlag)
                                {
                                    bam_aux_update_str(bampair[1], TagName, strlen(TagContent->mem + TagContent->offset_ar[1]) + 1, TagContent->mem + TagContent->offset_ar[1]);
                                }
                                i++;
                                if(AddFlag)
                                {
                                    addStrToArrstrRepXC(AddArr, TagContent->mem + TagContent->offset_ar[1]);
                                }
                                arrstrFree(TagContent);
                                free(TagContent);
                                break;                                
                            }
                            case 2:
                            {
                                int BarcodeOrder = atoi(AddFieldArray->arr->mem + AddFieldArray->arr->offset_ar[i]);
                                arrstr_k * TagContent = splitStrByChar(BarChain->mem + BarChain->offset_ar[BarcodeOrder - 1], ':');
                                if(strcmp(TagContent->mem + TagContent->offset_ar[0], AddFieldArray->arr->mem + AddFieldArray->arr->offset_ar[i+1]) == 0)
                                {
                                    strncpy(TagName, AddFieldArray->arr->mem + AddFieldArray->arr->offset_ar[i + 2], 2);
                                    bam_aux_update_str(bampair[0], TagName, strlen(TagContent->mem + TagContent->offset_ar[1]) + 1, TagContent->mem + TagContent->offset_ar[1]);
                                    if(PairFlag)
                                    {
                                        bam_aux_update_str(bampair[1], TagName, strlen(TagContent->mem + TagContent->offset_ar[1]) + 1, TagContent->mem + TagContent->offset_ar[1]);
                                    }
                                }
                                i = i + 2;
                                if(AddFlag)
                                {
                                    addStrToArrstrRepXC(AddArr, TagContent->mem + TagContent->offset_ar[1]);
                                }
                                arrstrFree(TagContent);
                                free(TagContent);
                                break;
                            }
                            default:
                                break;
                        }
                        
                    }
                }
                arrstrFree(BarChain);
                free(BarChain);
            }

            arrstrFree(AllInfoArray);
            free(AllInfoArray);

            if(writebam > 0) 
            {
                if(sam_write1(fp, header, bampair[0]) < 0) warnings("Failed to write %s\n", SeqName);
                if(PairFlag) 
                {
                    if(sam_write1(fp, header, bampair[0]) < 0) warnings("Failed to write %s\n", SeqName);
                }
            }
        }
        PairFlag = 0;
    }
    for(int i = 0; i < 2; i++) bam_destroy1(bampair[i]);
    free(bampool.bam_record);
    free(bampool.bam_mem);
    free(WaitQueue);
    return NULL;
}

int bam_p2s_core_ext(arrsf_k * FieldArray, arrstr_k * SonArray, arrsf_k * AddFieldArray, const char * fn, const char * fnout, const char * modeout, size_t _max_mem, int n_threads, const htsFormat * in_fmt, 
                        const htsFormat * out_fmt, char * arg_list, int writebam, int rr, int maxlen, samFile * BamFileOut, int strict, int atac, int extpos, int extneg)
{
    sam_hdr_t * header = NULL;
    samFile * fp;
    kh_m2_t * har = kh_init(m2);

    if(n_threads < 2) n_threads = 1;
    int ThreadNum = n_threads;
    // 1 GB
    _max_mem = (1 << 30);
    fp = sam_open_format(fn, "r", in_fmt);
    if(fp == NULL) error("Can not open \"%s\"\n", fn); 
    header = sam_hdr_read(fp);
    if(header == NULL) error("Failed to read header from \"%s\"\n", fn);

    if(sam_hdr_write(BamFileOut, header) == -1) warnings("Write header failed");

    pthrs_k ThreadArray;
    initThread(&ThreadArray, ThreadNum);
    pthread_mutex_init(&slock, NULL);
    int FinishFlag = 0, RoundCount = 0;
    int MaxThread = 0, ProcessedRead = 0; 
    for(;;)
    {
        int ThreadCreated = 0;
        for(size_t i = 0; i < ThreadArray.l; i++)
        {
            if(ThreadArray.statear[i] == 1)
            {
                bam1_pool BamPool = loadSam(fp, header, _max_mem);
                if(BamPool.mem_full == 1) FinishFlag = 1;
                if(BamPool.count > 0)
                {
                    parnam_k SamConfig;
                    SamConfig.bampool = BamPool;
                    SamConfig.FieldArray = FieldArray;
                    SamConfig.SonArray = SonArray;
                    SamConfig.index = i;
                    SamConfig.AddFieldArray = AddFieldArray;
                    SamConfig.har = har;
                    SamConfig.header = header;
                    SamConfig.writebam = writebam;
                    SamConfig.rr = rr;
                    SamConfig.maxlen = maxlen;
                    SamConfig.BamOut = BamFileOut;
                    SamConfig.strict = strict;
                    SamConfig.atac = atac;
                    SamConfig.extneg = extneg;
                    SamConfig.extpos = extpos;
                    pthread_create(&ThreadArray.ptar[i], NULL, parseSam, &SamConfig);
                    ThreadCreated++;
                    ThreadArray.statear[i] = 0;
                    ProcessedRead = ProcessedRead + BamPool.count;
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
            pthread_join(ThreadArray.ptar[i], NULL);
            ThreadArray.statear[i] = 1;
        }
        processeinfo("%d.", ProcessedRead);
        if(FinishFlag == 1) break;
    }

    char * Out;
    if(strcmp(fnout, "-") == 0) Out = "cluster"; else Out = strConcateByChar((char *) fnout, "cluster", '.');
    FILE * fo = fopen(Out, "w");
    khint_t kk;

    fprintf(fo, "#This file is produced by barp p2s.\n");
    fprintf(fo, "@CN");
    if(SonArray)
    {
        for(int i = 0; i < SonArray->l; i++)
        {
            fprintf(fo, "\t%s", SonArray->mem + SonArray->offset_ar[i]);
        }
    }
    else
    {
        if(FieldArray)
        {
            for(int i = 0; i < FieldArray->flag->l; i++)
            {
                fprintf(fo, "\t%s", FieldArray->arr->mem + FieldArray->arr->offset_ar[i]);
                if(FieldArray->flag->ar[i] == 2) i++;
            }
        }
    }
    fprintf(fo, "\tCount");
    if(AddFieldArray != NULL)
    {
        for(int i = 0; i < AddFieldArray->arr->l/2; i++)
        {
            fprintf(fo, "\t%s", AddFieldArray->arr->mem + AddFieldArray->arr->offset_ar[i * 2 + 1]);
        }
    }
    fprintf(fo, "\n");
    for(kk = 0; kk < kh_end(har); ++kk)
    {
        if(kh_exist(har, kk))
        {
            m2v_k * m2v = kh_val(har, kk);
            for(int i = 0; i < m2v->fathers->l; i++)
            {
                fprintf(fo, "%s\t", m2v->fathers->mem + m2v->fathers->offset_ar[i]);
            }
            fprintf(fo, "%s\t%ld\t", m2v->name, m2v->opts->l);
            for(int i = 0; i < m2v->opts->l; i++)
            {
                fprintf(fo, "%s\t%li\t%li\t", sam_hdr_tid2name(header, m2v->opts->opts[i].tid), m2v->opts->opts[i].start, m2v->opts->opts[i].end);
            }
            for(int i = 0; i < m2v->add->l; i++)
            {
                fprintf(fo, "%s\t", m2v->add->mem + m2v->add->offset_ar[i]);
            }
            fprintf(fo, "\n");
        }
    }
    fclose(fo);
    if(strcmp(fnout, "-") != 0) free(Out);
    sam_hdr_destroy(header);
    return 0;
}

arrsf_k * parseFieldStr(char * fad)
{
    arrstr_k * Tmp = splitStrByChar(fad, '-');
    arrsf_k * Ret = (arrsf_k *) calloc(1, sizeof(arrsf_k));
    Ret->arr = (arrstr_k *) calloc(1, sizeof(arrstr_k));
    Ret->flag = (arwl_k *) calloc(1, sizeof(arwl_k));

    // flag = 0: ; flag = 1: CELL; flag = 2: CELL:CO
    for(int i = 0; i < Tmp->l; i++)
    {
        arrstr_k * Field = splitStrByChar(Tmp->mem + Tmp->offset_ar[i], ':');
        switch(Field->l)
        {
            case 1 :
            {
                addEleToArwlRep(Ret->flag, 1);
                addStrToArrstrXC(Ret->arr, Field->mem + Field->offset_ar[0]);
                break;
            }
            case 2 :
            {
                addEleToArwlRep(Ret->flag, 2); addEleToArwlRep(Ret->flag, 2);
                if(!addStrToArrstrXC(Ret->arr, Field->mem + Field->offset_ar[0])) warnings("%s added twice, will be considered as the last time its specified.\n", Field->mem + Field->offset_ar[0]);
                if(strlen(Field->mem + Field->offset_ar[1]) != 2) error("[p2s] Field add format error at \"%s\", field TAG should be a two characters string!\n", Tmp->mem + Tmp->offset_ar[i]);
                if(!addStrToArrstrXC(Ret->arr, Field->mem + Field->offset_ar[1])) error("[p2s] field-add format error at \"%s\", identifier or field TAG can not duplicate!\n", Tmp->mem + Tmp->offset_ar[i]);
                break;
            }
            default :
            {
                error("[p2s] Field add format error at \"%s\",if you want to write to bam, the parameter should be like CELL:CO-COMPLEX:CP\n", Tmp->mem + Tmp->offset_ar[i]);
            }
        }
        arrstrFree(Field);
        free(Field);
    }
    arrstrFree(Tmp);
    free(Tmp);
    return Ret;
}

arrstr_k * parseSonStr(char * son, arrstr_k * FieldArray)
{
    // son will be like CELL:COMPLEX
    arrstr_k * Tmp = (arrstr_k *) splitStrByChar(son, ':');
    arrstr_k * Tmp2 = arrstrCopy(FieldArray);
    arrstr_k * Ret = (arrstr_k *) calloc(1, sizeof(arrstr_k));
    if(Tmp->l < 1) error("[p2s] Son format error at \"%s\" should be like CELL:COMPLEX!\n", son);
    for(int i = 0; i < Tmp->l; i++)
    {
        if(!addStrToArrstrXC(Tmp2, Tmp->mem + Tmp->offset_ar[i]))
        {
            addStrToArrstrXC(Ret, Tmp->mem + Tmp->offset_ar[i]);
        }
        else
        {
            error("[p2s] %s in son not found in --field-add! should be like CELL:COMPLEX!\n", son);
        }
    }
    arrstrFree(Tmp);
    free(Tmp);
    arrstrFree(Tmp2);
    free(Tmp2);
    return Ret;
}

arrsf_k * parseFieldAddStr(char * adf)
{
    // The last flag shows if loop through all the barcodes.
    arrsf_k * Ret = (arrsf_k *) calloc(1, sizeof(arrsf_k));
    Ret->arr = (arrstr_k *) calloc(1, sizeof(arrstr_k));
    Ret->flag = (arwl_k *) calloc(1, sizeof(arwl_k));

    arrstr_k * AddFieldArray = splitStrByChar(adf, '-');
    char * CurADF;
    size_t LoopFlag = 0;
    // flag shows what this position represent : 0 : 2 field and order of barcode; 1 : 2 field and name of barcode; 2 : 3 field and order of barcode; 3 : 3 field and name of barcode; 4 : tagname.
    for(int i = 0; i < AddFieldArray->l; i++)
    {
        CurADF = AddFieldArray->mem + AddFieldArray->offset_ar[i];
        arrstr_k * PosArray = splitStrByChar(CurADF, ':');
        if(PosArray->l == 2)
        {
            if(strlen(PosArray->mem + PosArray->offset_ar[1]) != 2)
            {
                error("Extra field : \"%s\" format error, should be like : 3:UMI:UM-DPM:DP-5:CP", CurADF);
            }
            else
            {
                if(atoi(PosArray->mem + PosArray->offset_ar[0]) != 0)
                {
                    addEleToArwlRep(Ret->flag, 0);
                    addEleToArwlRep(Ret->flag, 4);
                }
                else
                {
                    addEleToArwlRep(Ret->flag, 1);
                    addEleToArwlRep(Ret->flag, 4);
                    LoopFlag = 1;
                }
                addStrToArrstrRepXC(Ret->arr, PosArray->mem + PosArray->offset_ar[0]);
                addStrToArrstrRepXC(Ret->arr, PosArray->mem + PosArray->offset_ar[1]);
            }
        }
        else if(PosArray->l == 3)
        {
            if(strlen(PosArray->mem + PosArray->offset_ar[2]) != 2)
            {
                error("Extra field : \"%s\" format error, should be like : 3:UMI:UM-DPM:DP-5:CP", CurADF);
            }
            if(atoi(PosArray->mem + PosArray->offset_ar[0]) != 0)
            {
                addEleToArwlRep(Ret->flag, 2);
                addEleToArwlRep(Ret->flag, 3);
                addEleToArwlRep(Ret->flag, 4);
                addStrToArrstrRepXC(Ret->arr, PosArray->mem + PosArray->offset_ar[0]);
                addStrToArrstrRepXC(Ret->arr, PosArray->mem + PosArray->offset_ar[1]);
                addStrToArrstrRepXC(Ret->arr, PosArray->mem + PosArray->offset_ar[2]);
            }
            else
            {
                error("Extra field : \"%s\" format error, should be like : 3:UMI:UM-DPM:DP-5:CP", CurADF);
            }


        }
        else
        {
            error("Extra field : \"%s\" format error, should be like : 3:UMI:UM-DPM:DP-5:CP", CurADF);
        }
        arrstrFree(PosArray);
        free(PosArray);
    }
    arrstrFree(AddFieldArray);
    free(AddFieldArray);
    addEleToArwlRep(Ret->flag, LoopFlag);
    return Ret;
}

static void p2s_usage(FILE *fp)
{
    fprintf(fp, "\n"
"Usage: barp p2s [options...] [in.bam]\n"
"Options:\n"
"  -b STR     Identifier name parse string. Example: UMI:MI-CELL:CO-COMPLEX:CP or CELL-COMPLEX. If separate with : will write to bam tag.\n"
"  -d STR     Additional field to output besides identifiers : should be like 3:UMI\n"
"  -s STR     Identifier belongings. Example: CELL:COMPLEX\n"
"  -w         Write output bam file\n"
"  -i         Write output bam file and complex file\n"
"  -m INT     Set maximum memory per thread; suffix K/M/G recognized [768M]\n"
"  --read1    Write output with read 1 position\n"
"  --read2    Write output with read 2 position\n"
"  --strict   Write coordinate if both read appear in the BAM file\n\n"
);

}

int parseToSamField(int argc, char * argv[])
{
    size_t max_mem = SORT_DEFAULT_MEGS_PER_THREAD << 20;
    int c, level = -1, nargs, strict = 0;
    char * fnout = "-", * fad = "-", * son = "-", *arg_list = NULL, * adf = "-";
    char  modeout[12];
    int writebam = 0, writeall = 0;
    int rr = 0, maxlen = 600, atac = 0;
    int extpos = 0, extneg = 0;
    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
    // field-add : CELL:CO-COMPLEX:CP output-by : CELL:COMPLEX-UMI
    kstring_t tmpprefix = {0, 0, NULL};
    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, 'O', 0, 0, '@'),
        {"field-add", required_argument, NULL, 'b'},
        {"output-by", required_argument, NULL, 's'},
        {"read1", 0, NULL, '1'},
        {"read2", 0, NULL, '2'},
        {"strict", 0, NULL, '3'},
        {"atac", 0, NULL, '4'},
        {NULL, 0, NULL, 0}
    };
    while((c = getopt_long(argc, argv, "l:m:o:O:T:@:ub:s:d:w3", lopts, NULL)) >= 0)
    {
        switch(c)
        {
            case 'o': fnout = optarg; break;
            case 'm': {
                char * q;
                max_mem = strtol(optarg, &q, 0);
                if( *q == 'k' || * q == 'K') max_mem <<= 10;
                else if(*q == 'm' || *q == 'M') max_mem <<= 20;
                else if(*q == 'g' || *q == 'G') max_mem <<= 30;
                break;
            }
            case '1': rr = 1; break;
            case '2': rr = 2; break;
            case '4': extpos = 4; extneg = -5; atac = 1; break;
            case 'd': adf = optarg; break;
            case 'T': kputs(optarg, &tmpprefix); break;
            case 'b': fad = optarg; break;
            case 's': son = optarg; break;
            case 'r': rr = atoi(optarg); break;
            case 'w': writebam = 2; break;
            case 'a': writeall = 1; break;
            case '3': strict = 1; break;
            case 't': maxlen = atoi(optarg); break;
            default: if(parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
            case '?': p2s_usage(stderr); exit(1);
        }
    }
    if(writeall) writebam = 1;
    // 0与进程的标 准输入相关联，文件描述符1与标准输出相关联，文件描述符2与标准出错相关联
    nargs = argc - optind;
    if(nargs == 0 && isatty(STDIN_FILENO))
    {
        p2s_usage(stdout);
        exit(1);
    }
    else if(nargs >= 2)
    {
        // If exactly two, user probably tried to specify legacy <out.prefix>
        if (nargs == 2)
            fprintf(stderr, "[p2s] Use -T PREFIX / -o FILE to specify temporary and final output files\n");
        p2s_usage(stderr);
        exit(1);
    }

    if (max_mem < (SORT_MIN_MEGS_PER_THREAD << 20))
    {
        complain_about_memory_setting(max_mem);
        exit(1);
    }
    strcpy(modeout, "wb");
    // open mode detected by file suffix;
    sam_open_mode(modeout+1, fnout, NULL);
    if(level >= 0) sprintf(strchr(modeout, '\0'), "%d", level < 9 ? level : 9);
    if(tmpprefix.l == 0)
    {
        if(strcmp(fnout, "-") != 0)
        {
            kputs(".tmp", &tmpprefix);
        }
        else
        {
            kputc('.', &tmpprefix);
        }
    }

    char * FileOut;
    if(strcmp(fnout, "-") == 0)
    {
        FileOut = "P2S.bam";
    }
    else
    {
        FileOut = strConcateByChar(fnout, "P2S.bam", '_');
    }
    samFile * BamOut = sam_open(FileOut, "w");
    if(strcmp(fnout, "-") != 0) free(FileOut);

    if(!(arg_list = stringify_argv(argc + 1, argv - 1)))
    {
        error("Failed to create arg_list\n");
    }
    // todo multi thread
    // if(stat(tmpprefix.s, &st) == 0 && S_ISDIR(st.st_mode))
    // {
    //     unsigned t = ((unsigned) time(NULL)) ^ ((unsigned) clock());
    //     if(tmpprefix.s[tmpprefix.l-1] != '/') kputc('/', &tmpprefix);
    //     ksprintf(&tmpprefix, "barp.%d.%u.tmp", (int) getpid(), t%10000);
    // }
    // writebam = 2; only parse to bam; = 1, write bam and complex line; = 0, only write complex line.
    arrsf_k * FieldArray = NULL;
    arrstr_k * SonArray = NULL;

    if(strcmp(fad, "-") == 0)
    {
        writebam = 2;
        warnings("Field of identifier not specified, will parse barcode only.");
    }
    else 
    {
        FieldArray = parseFieldStr(fad);
        if(FieldArray->arr->l > 1)
        {
            if(strcmp(son, "-") == 0) warnings("--output-by not specified! Will output each identifier separately!\n"); else SonArray = parseSonStr(son, FieldArray->arr);
        }
    }
    arrsf_k * AddFieldArray = NULL;
    if(strcmp(adf, "-") != 0) AddFieldArray =  parseFieldAddStr(adf);
    bam_p2s_core_ext(FieldArray, SonArray, AddFieldArray, (nargs > 0)? argv[optind] : "-", fnout, modeout, max_mem, ga.nthreads, &ga.in, &ga.out, arg_list, writebam, rr, maxlen, BamOut, strict, atac, extpos, extneg);
    sam_close(BamOut);
    return 0;
}