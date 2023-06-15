#include "kutils.h"
#include "parse_samfield.h"
#include <getopt.h>
#include <sys/stat.h>
#include "sctable.h"

static void rdp_usage(FILE *fp)
{
    fprintf(fp,
"Usage: barp rdp [options...] [in.cluster]\n"
"Options:\n"
"  -r STR     Remove duplicates by field.\n"
"  -f STR     Re-specify the representation of field.\n"
"  -o FILE    Write final output to FILE rather than standard output\n");
}

void addFragBefore(fna_k * fna, char * fraginfo, size_t fragmem, size_t pos)
{
    // printf("Add Frag Before %ld\n", pos);
    fna->mem = (char *) realloc(fna->mem, fna->m + fragmem);
    memmove(fna->mem + fna->m, fraginfo, fragmem);
    size_t * tmp = (size_t *) calloc(fna->l, sizeof(size_t));
    memmove(tmp, fna->offset_ar, fna->l * sizeof(size_t));
    fna->offset_ar = (size_t *) realloc(fna->offset_ar, (fna->l + 1) * sizeof(size_t));
    memmove(fna->offset_ar, tmp, pos * sizeof(size_t));
    fna->offset_ar[pos] = fna->m;
    memmove(fna->offset_ar + pos + 1, tmp + pos, sizeof(size_t) * (fna->l - pos));
    free(tmp);
    fna->m = fna->m + fragmem;
    fna->l++;
}
// pos is start from 0
void addFragAfter(fna_k * fna, char * fraginfo, size_t fragmem, size_t pos)
{
    // printf("Add Frag After %ld\n", pos);
    fna->mem = (char *) realloc(fna->mem, fna->m + fragmem);
    memmove(fna->mem + fna->m, fraginfo, fragmem);
    size_t * tmp = (size_t *) calloc(fna->l, sizeof(size_t));
    memmove(tmp, fna->offset_ar, fna->l * sizeof(size_t));
    fna->offset_ar = (size_t *) realloc(fna->offset_ar, (fna->l + 1) * sizeof(size_t));
    memmove(fna->offset_ar, tmp, (pos + 1) * sizeof(size_t));
    fna->offset_ar[pos + 1] = fna->m;
    memmove(fna->offset_ar + pos + 2, tmp + pos + 1, sizeof(size_t) * (fna->l - pos - 1));
    free(tmp);
    fna->m = fna->m + fragmem;
    fna->l++;
}

int findFragInFna(fna_k * fna, char * fraginfo, size_t fragmem)
{
    if(fna->l == 0) 
    {
        fna->mem = realloc(fna->mem, fragmem);
        fna->m = fragmem;
        memmove(fna->mem, fraginfo, fragmem);
        fna->l = 1;
        fna->offset_ar = (size_t *) realloc(fna->offset_ar, 1 * sizeof(size_t));
        fna->offset_ar[0] = 0;
        return 1;
    }
    else
    {
        // judge the memory big or small; small 
        size_t cmpbyte = CPUBYTENUM;
        int startpoint = 0, endpoint = fna->l - 1, midpoint;
        char * fnasptr, * fnamptr, * fnaeptr;
        char * fragptr = fraginfo;
        int cmps, cmpe, cmpm, added = 0;
        int startflag = 1, endflag = 1;
        while (cmpbyte <= fragmem)
        {
            // printf("1 - Start : %d - End : %d - Mid : %d\n", startpoint, endpoint, midpoint);
            // printf("Frag : %s\t%s\n", fragptr, fragptr + strlen(fragptr) + 1);
            if(startflag)
            {
                fnasptr = fna->mem + fna->offset_ar[startpoint];
                // printf("Start Str : %s\t%s\n", fnasptr, fnasptr + strlen(fnasptr) + 1);
                cmps = memcmp(fragptr, fnasptr, CPUBYTENUM);
                if(cmps < 0)
                {
                    added = 1;
                    addFragBefore(fna, fraginfo, fragmem, startpoint);
                    return 1;
                }
            }
            if(endflag)
            {
                fnaeptr = fna->mem + fna->offset_ar[endpoint];
                // printf("End Str : %s\t%s\n", fnaeptr, fnaeptr + strlen(fnaeptr) + 1);
                cmpe = memcmp(fnaeptr, fragptr, CPUBYTENUM);
                if(cmpe < 0)
                {
                    added = 1;
                    addFragAfter(fna, fraginfo, fragmem, endpoint);
                    return 1;
                }
            }
            midpoint = (startpoint + endpoint)/2;
            fnamptr = fna->mem + fna->offset_ar[midpoint];
            // printf("Mid Str : %s\t%s\n", fnamptr, fnamptr + strlen(fnamptr) + 1);
            cmpm = memcmp(fragptr, fnamptr, CPUBYTENUM);
            // printf("2 - Start : %d - End : %d - Mid : %d - cmpm : %d\n", startpoint, endpoint, midpoint, cmpm);
            if(cmpm < 0)
            {
                endpoint = midpoint - 1;
                // no need to compare start point
                startflag = 0;
            }
            else if(cmpm > 0)
            {
                startpoint = midpoint + 1;
                // no need to compare end point
                endflag = 0;
            }
            else
            {
                while(1)
                {
                    // printf("Deeping Down\n");
                    if(cmpbyte == fragmem) return 0;
                    cmpbyte += CPUBYTENUM;
                    fragptr = fragptr + CPUBYTENUM;
                    fnamptr = fnamptr + CPUBYTENUM;
                    cmpm = memcmp(fragptr, fnamptr, CPUBYTENUM);
                    if(cmpm != 0) break;
                }
                if(cmpm < 0)
                {
                    endpoint = midpoint - 1;
                    startflag = 0;
                    fragptr = fraginfo;
                    cmpbyte = CPUBYTENUM;
                }
                else
                {
                    startpoint = midpoint + 1;
                    endflag = 0;
                    fragptr = fraginfo;
                    cmpbyte = CPUBYTENUM;
                }
            }
            // printf("3 - Start : %d - End : %d - Mid : %d - cmpm : %d\n", startpoint, endpoint, midpoint, cmpm);
            if(startpoint > endpoint)
            {
                if(cmpm < 0)
                {
                    addFragBefore(fna, fraginfo, fragmem, midpoint);
                    return 1;
                }
                else
                {
                    addFragAfter(fna, fraginfo, fragmem, midpoint);
                    return 1;
                }
            }
        }
    }
}

int removeCore(char * fi, arwl_k * col, char * fo)
{
    FILE * fn = fopen(fi, "r");
    FILE * fp = fopen(fo, "a");
    size_t skip = col->ar[col->l - 1];
    char * line;
    size_t Duplicates = 0;
    size_t TotalReads = 0;
    size_t LineCount = 0;

    while(1)
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
                if(LineCount%10000 == 0) processeinfo("%ld.", LineCount);;
                arrstr_k * InfoArray = splitStrByChar(line, '\t');
                arrstr_k * AddArray = (arrstr_k *) calloc(1, sizeof(arrstr_k));
                arwl_k * Reserve = (arwl_k *) calloc(1, sizeof(arwl_k));
                Reserve->ar = (size_t *)calloc(1, sizeof(size_t));
                size_t Count = 0;
                if(InfoArray->l > skip)
                {
                    for(int i = 0; i < skip - 1; i++)
                    {
                        fprintf(fp, "%s\t", InfoArray->mem + InfoArray->offset_ar[i]);
                    }
                    // printf("\n");
                    size_t FragCount;
                    if((FragCount = atoi(InfoArray->mem + InfoArray->offset_ar[skip - 1])) != 0)
                    {
                        fna_k FragArray = {NULL, NULL, 0, 0};
                        char * FragInfo = NULL;
                        for(int i = 0; i < FragCount; i++)
                        {
                            size_t InitMem = 0;
                            size_t FragMem = 0;
                            for(int j = 0; j < col->l - 1; j++)
                            {
                                FragInfo = realloc(FragInfo, FragMem + strlen(InfoArray->mem + InfoArray->offset_ar[col->ar[j] + i * 3 + skip]) + 1);
                                strcpy(FragInfo + FragMem, InfoArray->mem + InfoArray->offset_ar[col->ar[j] + i * 3 + skip]);
                                FragMem = FragMem + strlen(InfoArray->mem + InfoArray->offset_ar[col->ar[j] + i * 3 + skip]);
                                // printf(" %s\t", InfoArray->mem + InfoArray->offset_ar[col->ar[j] + i * 3 + skip]);
                                FragInfo[FragMem++] = '\0';
                            }
                            while(FragMem > InitMem) InitMem += CPUBYTENUM;
                            FragMem = InitMem;
                            TotalReads++;
                            if(findFragInFna(&FragArray, FragInfo, FragMem) == 1)
                            {
                                Count++;
                                addEleToArwlRep(Reserve, i);
                            }
                            else
                            {
                                Duplicates++;
                            }
                        }
                        if(FragArray.offset_ar) free(FragArray.offset_ar);
                        if(FragArray.mem) free(FragArray.mem);
                    }
                    else
                    {
                        continue;
                    }
                }
                fprintf(fp, "%ld\t", Count);
                if(Reserve->l > 0)
                {
                    for(int i = 0; i < Reserve->l - 1; i++)
                    {
                        fprintf(fp, "%s\t", InfoArray->mem + InfoArray->offset_ar[Reserve->ar[i] * 3 + skip]);
                        fprintf(fp, "%s\t", InfoArray->mem + InfoArray->offset_ar[Reserve->ar[i] * 3 + 1 + skip]);
                        fprintf(fp, "%s\t", InfoArray->mem + InfoArray->offset_ar[Reserve->ar[i] * 3 + 2 + skip]);
                    }
                    fprintf(fp, "%s\t", InfoArray->mem + InfoArray->offset_ar[Reserve->ar[Reserve->l - 1] * 3 + skip]);
                    fprintf(fp, "%s\t", InfoArray->mem + InfoArray->offset_ar[Reserve->ar[Reserve->l - 1] * 3 + 1 + skip]);
                    fprintf(fp, "%s\n", InfoArray->mem + InfoArray->offset_ar[Reserve->ar[Reserve->l - 1] * 3 + 2 + skip]);
                }
                arrstrFree(InfoArray);
                free(InfoArray);
                arrstrFree(AddArray);
                free(AddArray);
                kaFree(Reserve);
                free(Reserve);
            }
            free(line);
        }
        else
        {
            break;
        }
    }
    fclose(fp);
    warnings("Duplicate : %ld - Total Fragment : %ld\n", Duplicates, TotalReads);
}

int removeDuplicate(int argc, char * argv[])
{
    int ret, nargs;
    char c;
    rpst_k * rpst = initRpst();
    struct stat st;
    arwl_k * col;
    char * rmstr = "START", * fnout = "-", * fieldn = "-"; 
    static const struct option lopts[] = {
        {"remove-by", required_argument, NULL, 'r'},
        {"field-order", required_argument, NULL, 'f'},
        {NULL, 0, NULL, 0}
    };
    int count = 0;
    while((c = getopt_long(argc, argv, "r:f:o:", lopts, NULL)) >= 0)
    {
        switch(c)
        {
            case 'o': fnout = optarg; break;
            case 'r': rmstr = optarg; break;
            case 'f': fieldn = optarg; break;
            case '?': rdp_usage(stderr); ret = EXIT_FAILURE; goto RDP_END;
        }
    }
    nargs = argc - optind;
    if(nargs == 0 && isatty(STDIN_FILENO))
    {
        rdp_usage(stdout);
        ret = EXIT_SUCCESS;
        goto RDP_END;
    }
    char * fi = argv[optind], split = '\t';
    if(strcmp(fieldn, "-") == 0) fieldn = readHeader(fi); else split = '-';
    if(!parseNewField(fieldn, rpst, split)) col = parseCol(rmstr, rpst);
    char * result;
    if(strcmp(fnout, "-") == 0) result = "Dedup.cluster"; else result = strConcateByChar(fnout, "Dedup.cluster", '.');
    FILE * fo = fopen(result, "w");
    fprintf(fo, "#This file is produced by barp rdp.\n");
    fprintf(fo, "@CN");
    for(int i = 0; i < rpst->fathers->l; i++)
    {
        fprintf(fo, "\t%s", rpst->fathers->mem + rpst->fathers->offset_ar[i]);
    }
    fprintf(fo, "\tCount");
    for(int i = 0; i < rpst->addar->l; i++)
    {
        fprintf(fo, "\t%s", rpst->addar->mem + rpst->addar->offset_ar[i]);
    }
    fprintf(fo, "\n");
    fclose(fo);
    if(col->l < 2) error("No column specified to remove duplicates"); else removeCore(fi, col, result);
    if(strcmp(fnout, "-") != 0) free(result);
    RpstFree(rpst);
    return 0;
RDP_END:
    return 1;
}