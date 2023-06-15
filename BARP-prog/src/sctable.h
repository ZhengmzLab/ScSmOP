#ifndef _SCTABLE_H_
#define _SCTABLE_H_

#include "kutils.h"

typedef struct rdpst
{
    // field record number
    size_t frn;
    arrstr_k * fathers;
    arrstr_k * addar;
    optar_k * opts;
}rpst_k;

static inline rpst_k * initRpst()
{
    rpst_k * rpst = (rpst_k *) calloc(1, sizeof(rpst_k));
    rpst->addar = (arrstr_k *) calloc(1, sizeof(arrstr_k));
    rpst->fathers = (arrstr_k *) calloc(1, sizeof(arrstr_k));
    return rpst;
}

static inline void RpstFree(rpst_k * s)
{
    if(s)
    {
        if(s->addar)
        {
            arrstrFree(s->addar);
            free(s->addar);
        }
        if(s->fathers)
        {
            arrstrFree(s->fathers);
            free(s->fathers);
        }
    }
}

typedef struct FragNodeArray
{
    size_t * offset_ar;
    char * mem;
    size_t m;
    size_t l;
}fna_k;

static inline char * convertOpt(char * s)
{
    int i;
    char * t = s;
    for(i = 0; i < strlen(s); i++)
    {
        if(! isspace(s[i]))
        {
            if(s[i] >= 65 && s[i] <= 90)
            {
                * t = s[i] + 32;
            }
            else if(s[i] >= 97 && s[i] <= 122)
            {
                * t = s[i];
            }
            t++;
        }
    }
    * t = '\0';
    return s;
}

// start from start to additional field
static inline arwl_k * parseCol(char * rmstr, rpst_k * rpst)
{
    arwl_k * Ret = (arwl_k *) calloc(1, sizeof(arwl_k));
    size_t FieldBefore = rpst->fathers->l + 1;
    arrstr_k * Tmp = splitStrByChar(rmstr, '-');
    for(int i = 0; i < Tmp->l; i++)
    {
        convertOpt(Tmp->mem + Tmp->offset_ar[i]);
        // chr 0 start 1 end 2 add 3... 
        for(int j = 0; j < rpst->addar->l; j++)
        {
            if(strcmp(Tmp->mem + Tmp->offset_ar[i], rpst->addar->mem + rpst->addar->offset_ar[j]) == 0)
            {
                if(addEleToArwlXCSorted(Ret, j + 3) == -1) warnings("%s existed more than once, will treated as once\n", Tmp->mem + Tmp->offset_ar[i]);
            }
        }
        if(strcmp(Tmp->mem + Tmp->offset_ar[i], "start") == 0) 
        {
            addEleToArwlXCSorted(Ret, 0);
            addEleToArwlXCSorted(Ret, 1);
        }
        if(strcmp(Tmp->mem + Tmp->offset_ar[i], "end") == 0)
        {
            addEleToArwlXCSorted(Ret, 0);
            addEleToArwlXCSorted(Ret, 2);
        }
    }
    arrstrFree(Tmp);
    free(Tmp);
    arwl_k * Final = selectSortArray(Ret, 0);
    addEleToArwlRep(Final, FieldBefore);
    kaFree(Ret);
    free(Ret);
    return Final;
}

static inline char * getLine(FILE * fp)
{
    char c;
    kstring_t Tmp = {0, 0, NULL};
    while((c = fgetc(fp)) != EOF)
    {
        if(c != '\n') kputc(c, &Tmp); else break;
    }
    if(Tmp.l > 0) return Tmp.s;
    return NULL;
}

static inline char * readHeader(char * fn)
{
    // printf("Read Header\n");
    FILE * fi = fopen(fn, "r");
    if(fi == NULL) error("Failed to read %s\n", fn);
    while(1)
    {
        char * line = getLine(fi);
        if(strncmp(line, "@CN", 3) == 0)
        {
            fclose(fi);
            return line; 
        } 
        else
        {
            free(line);
        } 
    }
    fclose(fi);
    return "-";
}

static inline void writeHeader(char * fn, char * fo)
{
    FILE * fp = fopen(fo, "w");
    FILE * fi = fopen(fn, "r");
    if(fi == NULL) error("Failed to read %s\n", fn);
    while(1)
    {
        char * line = getLine(fi);
        if(strncmp(line, "@CN", 3) == 0)
        {
            fclose(fi);
            fprintf(fp, "%s\n", line);
        } 
        else
        {
            free(line);
        } 
    }
    fclose(fi);
    fclose(fp);
}

static inline int parseNewField(char * fdn, rpst_k * rpst, char split)
{
    arrstr_k * Tmp = splitStrByChar(fdn, split);
    int CountFlag = 1;
    for(int i = 1; i < Tmp->l; i++)
    {
        if(CountFlag)
        {
            if(strcmp(convertOpt(Tmp->mem + Tmp->offset_ar[i]), "count") == 0)
            {
                CountFlag = 0;
                continue;
            }
            addStrToArrstrRepXC(rpst->fathers, convertOpt(Tmp->mem + Tmp->offset_ar[i]));
        }
        else
        {
            addStrToArrstrRepXC(rpst->addar, convertOpt(Tmp->mem + Tmp->offset_ar[i]));
        }

    }
    arrstrFree(Tmp);
    free(Tmp);
    return 0;
}

#endif