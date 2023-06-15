#include <stdio.h>
#include <stdlib.h>
#include "parse_samfield.h"
#include "kutils.h"
#include "sctable.h"
#include <getopt.h>
#include <sys/stat.h>
// offset 0 : chr offset; 1 start; 2 end
typedef struct FragChrCoord
{
    uint64_t s, e;
    char * chr;
}FragChrCoord;

typedef struct FragChrCoordArr
{
    FragChrCoord * frag;
    // additional field
    arrstr_k * add;
    // additional field num per fragment, ml : merged length.
    size_t l, adsize, ml;
}FragChrCoordArr;

typedef struct Interval
{
    uint64_t s, e;
} Interval;

void FragChrCoordArrFree(FragChrCoordArr * s)
{
    if(s)
    {
        if(s->frag)
        {
            for(int i = 0; i < s->l; i++)
            {
                if(s->frag[i].chr) 
                {
                    free(s->frag[i].chr);
                    s->frag[i].chr = NULL;
                }
            }
        }
        if(s->add)
        {
            arrstrFree(s->add);
            free(s->add);
        }
        free(s->frag);
    }
}

static void mergefragment_usage(FILE *fp)
{
    fprintf(fp, "\n"
    "Usage: barp rdp [options...] [in.cluster]\n"
    "Options:\n"
    "  -d INT     Merge distance.\n"
    "  -e STR     Extend 3' end\n"
    "  -o FILE    Write final output to FILE rather than standard output\n\n");
}

int cmpfunc(const void * a, const void * b)
{
    FragChrCoord * DataOne = (FragChrCoord *) a;
    FragChrCoord * DataTwo = (FragChrCoord *) b;
    return strcmp(DataOne->chr, DataTwo->chr) >= 0 ? strcmp(DataOne->chr, DataTwo->chr) > 0 ? 1 : DataOne->s - DataTwo->s : -1; 
}

// return the length of merged intervals
void mergeIntervals(FragChrCoordArr * FragArr, char * delim, int64_t ext)
{
    FragChrCoord * arr = FragArr->frag;
    // arrstr_k * AddArr = (arrstr_k *) calloc(1, sizeof(arrstr_k));
    // sort interval in increasing order of start position
    qsort(arr, FragArr->l, sizeof(FragChrCoord), cmpfunc);
    // store the last element of output array
    int index = 0;
    arrstr_k * NewAdd = (arrstr_k *) calloc(1, sizeof(arrstr_k));
    // traverse all input Intervals
    kstring_t * MergedAdd = (kstring_t *) calloc(FragArr->adsize, sizeof(kstring_t));
    // for( int i = 0; i < FragArr->l; i++)
    // {
    //     printf("arr[%d]:%lu\n", i, arr[i].e);
    // }
    for(int i = 0; i < FragArr->l; i++)
    {
        // int64_t tmp = (int64_t) (arr[i].s - ext) > 0 ? arr[i].s - ext : 0;
        int64_t tmp = (int64_t) arr[i].s;
        // printf("Start : %?ld\n", tmp);
        if((arr[index].e + ext >= tmp) && (strcmp(arr[index].chr, arr[i].chr) == 0))
        {
            // printf("Ith End : %lu - Index End %lu\n", arr[i].e, arr[index].e);
            arr[index].e = MaxUint64(arr[index].e, arr[i].e);
            // printf("Index End : %lu\n", arr[index].e);
            for(int j = 0; j < FragArr->adsize; j++)
            {
                if(i != 0) kputs(delim, MergedAdd);
                kputs(FragArr->add->mem + FragArr->add->offset_ar[i * 2 + j] , MergedAdd);
            }
        }
        else
        {
            index++;
            arr[index] = arr[i];
            // printf("Index : %d - End : %ld\n",index, arr[index].s);
            addStrToArrstrRepXC(NewAdd, MergedAdd->s);
            ks_free(MergedAdd);
        }
    }
    arrstrFree(FragArr->add);
    free(FragArr->add);
    FragArr->add = NewAdd;
    FragArr->ml = index + 1;
}

int mergeCore(char * fi, arwl_k * col, FILE * fp, char * delim, int64_t ext, int extend, int minlen)
{
    FILE * fn = fopen(fi, "r");
    size_t skip = col->ar[col->l - 1];
    char * line;
    size_t LineCount = 0;
    size_t FragCount;
    char * over;
    while (1)
    {
        if((line = getLine(fn)) != NULL)
        {
            if(*line == '#' || *line == '@')
            {
                continue;
            }
            else
            {
                LineCount++;
                if(LineCount%10000 == 0) processeinfo("%ld.", LineCount);
                arrstr_k * InfoArray = splitStrByChar(line, '\t');
                int QuaFragCount = 0;
                if(InfoArray->l > skip)
                {
                    // printf("Count : %s\n", InfoArray->mem + InfoArray->offset_ar[skip - 1]);
                    if((FragCount = atoi(InfoArray->mem + InfoArray->offset_ar[skip - 1])) > 0)
                    {
                        FragChrCoordArr * FragArr = (FragChrCoordArr *) calloc(1, sizeof(FragChrCoordArr));
                        FragArr->add = (arrstr_k *) calloc(1, sizeof(arrstr_k));
                        FragArr->frag = (FragChrCoord *) calloc(FragCount, sizeof(FragChrCoord));
                        FragArr->adsize = col->l - 4;
                        for(int i = 0; i < FragCount; i++)
                        {
                            // printf("%s - %s - %s\n", InfoArray->mem + InfoArray->offset_ar[col->ar[0] + i * 3 + skip], InfoArray->mem + InfoArray->offset_ar[col->ar[1] + i * 3 + skip], InfoArray->mem + InfoArray->offset_ar[col->ar[2] + i * 3 + skip]);
                            char * Chr = strCopy(InfoArray->mem + InfoArray->offset_ar[col->ar[0] + i * 3 + skip]);
                            uint64_t Start = strtoul(InfoArray->mem + InfoArray->offset_ar[col->ar[1] + i * 3 + skip], &over, 10);
                            uint64_t End = strtoul(InfoArray->mem + InfoArray->offset_ar[col->ar[2] + i * 3 + skip], &over, 10);
                            if(End - Start < minlen)
                            {
                                // printf("Length Less Than %d\n", minlen);    
                                // free(Chr);
                                continue;
                            } 
                            else
                            {
                                // printf("Extend : %d\n", extend);
                                FragArr->frag[QuaFragCount].chr = Chr;
                                FragArr->frag[QuaFragCount].s = Start;
                                FragArr->frag[QuaFragCount].e = End + extend;
                                QuaFragCount++;
                            }
                            // printf("Chr : %s\n", InfoArray->mem + InfoArray->offset_ar[col->ar[0] + i * 3 + skip]);
                            // printf("Start : %lu\n", FragArr->frag[i].s);
                            // printf("End : %lu\n", FragArr->frag[i].e);
                            for(int j = 3; j < col->l - 1; j++)
                            {
                                // printf("Have Add\n");
                                addStrToArrstrRepXC(FragArr->add, InfoArray->mem + InfoArray->offset_ar[col->ar[j] + i * 3 + skip]);
                            }
                        }
                        FragArr->l = QuaFragCount;
                        // printf("FragArr : %d\n", QuaFragCount);
                        if(QuaFragCount > 0)
                        {
                            mergeIntervals(FragArr, delim, ext);
                            // printf("Merged Length: %ld\n", FragArr->ml);
                            for(int i = 0; i < skip - 1; i++)
                            {
                                fprintf(fp, "%s\t", InfoArray->mem + InfoArray->offset_ar[i]);
                            }
                            fprintf(fp, "%ld\t", FragArr->ml);
                            for(int i = 0; i < FragArr->ml; i++)
                            {
                                fprintf(fp, "%s\t%lu\t%lu\t", FragArr->frag[i].chr, FragArr->frag[i].s, FragArr->frag[i].e);
                            }
                            for(int i = 0; i < FragArr->ml; i++)
                            {
                                for(int j = 0; j < FragArr->adsize; j++)
                                {
                                    fprintf(fp, "%s\t", FragArr->add->mem + FragArr->add->offset_ar[2 * i + j]);
                                }
                            }
                            fprintf(fp,"\n");
                        }
                    }
                    else
                    {
                        warnings("%s", line);
                    }
                }
                arrstrFree(InfoArray);
                free(InfoArray);
            }
            free(line); 
        }
        else
        {
            break;
        }
    }
    fclose(fp);
    return 1;
}

int mergeFragment(int argc, char * argv[])
{
    int64_t ext = 0;
    // extend to 5' end .
    int nargs, extend = 0, minlen = 0;
    char * adf = "-", *fieldn = "-", * over;
    arwl_k * col;
    char * delim = ";";
    static const struct option lopts[] = {
        {"delim", required_argument, NULL, 'm'},
        {NULL, 0, NULL, 0}
    };
    rpst_k * rpst = initRpst();
    char c;
    while((c = getopt_long(argc, argv, "d:c:f:m:e:l:", lopts, NULL)) >= 0)
    {
        switch(c)
        {
            case 'd': ext = strtol(optarg, &over, 10); break;
            case 'c': adf = optarg; break;
            case 'f': fieldn = optarg; break;
            case 'm': delim = optarg; break;
            case 'e': extend = atoi(optarg); break;
            case 'l': minlen = atoi(optarg); break;
            case '?': mergefragment_usage(stderr); exit(1);
        }
    }
    if(extend < 0) extend = 0;
    if(minlen < 0) minlen = 0;
    nargs = argc - optind;
    if(nargs == 0 && isatty(STDIN_FILENO))
    {
        mergefragment_usage(stdout);
        return 1;
    }
    char * fi = argv[optind], split = '\t';
    if(strcmp(fieldn, "-") == 0) fieldn = readHeader(fi); else split = '-';
    // if(fieldn == NULL) error("Field not found in input file header!\n");
    char * ChrCoord = "CHR-START-END";
    char * ColStr = strConcateByChar(ChrCoord, adf, '-');
    if(!parseNewField(fieldn, rpst, split)) col = parseCol(ColStr, rpst);
    for(int i = 0; i < rpst->fathers->l; i++)
    {
        printf("%s\n", rpst->fathers->mem + rpst->fathers->offset_ar[i]);
    }
    for(int i = 0; i < rpst->addar->l; i++)
    {
        printf("%s\n", rpst->addar->mem + rpst->addar->offset_ar[i]);
    }
    for(int i = 0; i < col->l; i++)
    {
        printf("%ld\n", col->ar[i]);
    }
    printf("Col : %ld\n", col->l);
    char * result = "Merged_Cluster.txt";
    FILE * fp = fopen(result, "w");
    fprintf(fp, "#This file is produced by barp merge from %s.\n", fi);
    fprintf(fp, "@CN\n");
    mergeCore(fi, col, fp, delim, ext, extend, minlen);
    return 0;
}