#include "read_json.h"
#include <unistd.h>
#include "kstring.h"
#include "kson.h"
#include <ctype.h>
#include <stdio.h>
#include "parse_config.h"


#define BARCODECHAIN 0
#define BARCODETYPE 1
#define IDENTIFIER 2
#define SON 3
#define BUNIT 4

void barcodeFree(bar_k * curbar)
{
    if(curbar)
    {
        if(curbar->na) 
        {
            free(curbar->na);
            curbar->na = NULL;
        }
        if(curbar->rec)
        {
            arrstrFree(curbar->rec);
            free(curbar->rec);
            curbar->rec = NULL;
        }
        if(curbar->hash) 
        {
            hashDestroyB(curbar->hash);
            curbar->hash = NULL;
        }
    }
}

void barcodeChainFree(barch_k * barch, arrstr_k * freedBar, arrstr_k * allBar)
{
    if(barch)
    {
        barh_k * curbarh = barch->barh;
        for(int j = 0; j < barch->l; j++)
        {
            bar_k * curbar = curbarh->bar;
            for(int i = 0; i < curbarh->l; i++)
            {
                int FreeFlag = 0;
                for(int q = 0; q < allBar->l; q++)
                {
                    if(strcmp(allBar->mem + allBar->offset_ar[q], curbar->na) == 0) FreeFlag = 1;
                }
                if(FreeFlag)
                {
                    if(addStrToArrstrXC(freedBar, curbar->na) != 0)
                    {
                        barcodeFree(curbar);
                    }
                }
                curbar++;
            }
            free(curbarh->bar);
            curbarh++;
        }
        free(barch->barh);
    }
}

void idfFree(id_k * s)
{
    if(s)
    {
        if(s->v)
        {
            kaFree(s->v);
            free(s->v);
            s->v = NULL;
        }
        if(s->na) 
        {
            // printf("IDF NAME: %s\n", s->na);
            free(s->na);
            s->na = NULL;
        }
    }
}

void barclFree(barcl_k * s)
{
    if(s)
    {   
        arrstr_k * AllBarcode = (arrstr_k *) calloc(1, sizeof(arrstr_k));
        for(int i = 0; i < s->l; i++)
        {
            for(int j = 0; j < s->barch[i].l; j++)
            {
                for(int m = 0; m < s->barch[i].barh[j].l; m++)
                {
                    addStrToArrstrXC(AllBarcode, s->barch[i].barh[j].bar[m].na);
                }
            }
        }
        arrstr_k * FreedBarcode = (arrstr_k *) calloc(1, sizeof(arrstr_k));
        if(s->barch)
        {
            barch_k * curbarch = s->barch;
            for(int i = 0; i < s->l; i++)
            {
                barcodeChainFree(curbarch, FreedBarcode, AllBarcode);
                curbarch++;
            }
            free(s->barch);
        }
        arrstrFree(FreedBarcode);
        free(FreedBarcode);

        arrstrFree(AllBarcode);
        free(AllBarcode);

        if(s->jump){
            kaFree(s->jump);
            free(s->jump);
        }

        if(s->idf)
        {
            id_k * curidf = s->idf;
            for(int i = 0; i < s->idl; i++)
            {
                idfFree(curidf);
                curidf++;
            }
            free(s->idf);
        }
        if(s->sonch)
        {
            arrstrFree(s->sonch);
            free(s->sonch);
        }
    }
}
// remove the empty words in json option and convert them to lower case.
int optConvert(char * s)
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
    char * opt[] = {
        "barcodechain",
        "barcodetype",
        "identifier",
        "son",
        NULL,
    };
    char ** ptr;
    size_t num = 0;
    for(ptr = opt; * ptr != NULL; ptr++)
    {
        if(strcmp(*ptr, s) == 0) return num; else num++;
    }
    return -1;
}

barch_k parseBarcodeChain(const kson_node_t * p)
{
    barch_k ret = {0};
    arrstr_k * BarcodeArray = splitStrByChar(p->v.child[0]->key, '-');
    ret.l = BarcodeArray->l;
    ret.barh = (barh_k *) calloc(ret.l, sizeof(barh_k));
    for(int i = 0; i < ret.l; i++)
    {
        arrstr_k * BarHole = splitStrByChar(BarcodeArray->mem + BarcodeArray->offset_ar[i], '/');
        ret.barh[i].l = BarHole->l;
        ret.barh[i].bar = (bar_k *) calloc(BarHole->l, sizeof(bar_k));
        for(int j = 0; j < BarHole->l; j++) ret.barh[i].bar[j].na = strCopy(BarHole->mem + BarHole->offset_ar[j]);
    }

    arrstrFree(BarcodeArray);
    free(BarcodeArray);

    arrstr_k * rdps = splitStrByChar(p->v.child[0]->v.str, ':');
    if(rdps->l != 2) error("Barcode chain position format error! Should be like \"R1:1\"!");
    arrstr_k * rd = splitStrByChar(rdps->mem + rdps->offset_ar[0], 'R');
    if(!rd)
    {
        error("Barchain start position format error! Should be started with R and followed by 1 or 2 or 3 or 4!");
    }
    else
    {
        switch ((rd->mem + rd->offset_ar[0])[0])
        {
            case '1':
                ret.rdp.rd = 1;
                break;
            case '2':
                ret.rdp.rd = 2;
                break;
            case '3':
                ret.rdp.rd = 3;
                break;
            case '4':
                ret.rdp.rd = 4;
                break;
            default:
                error("Currently SCSM only support for at most 4 read file!\n");
                break;
        }
    }
    arrstrFree(rd);
    free(rd);
    if(atoi(rdps->mem + rdps->offset_ar[1])) ret.rdp.pos = atoi(rdps->mem + rdps->offset_ar[1]); else error("Barcode chain position format error: read position not specified! Should be like \"R1:1\"!\n");
    arrstrFree(rdps);
    free(rdps);
    return ret;
}

void rearrangeArwlByReadPos(arwl_k * s, barch_k * barch)
{
    if(s)
    {
        arwl_k * tmp1 = (arwl_k *) calloc(1, sizeof(arwl_k));
        arwl_k * tmp2 = (arwl_k *) calloc(1, sizeof(arwl_k));
        arwl_k * tmp3 = (arwl_k *) calloc(1, sizeof(arwl_k));
        arwl_k * tmp4 = (arwl_k *) calloc(1, sizeof(arwl_k));

        for(int i = 0; i < s->l; i++)
        {
            switch(barch[s->ar[i]].rdp.rd)
            {
                case 1:
                    addEleToArwlXC(tmp1, s->ar[i]);
                    break;
                case 2:
                    addEleToArwlXC(tmp2, s->ar[i]);
                    break;
                case 3:
                    addEleToArwlRep(tmp3, s->ar[i]);
                    break;
                case 4:
                    addEleToArwlRep(tmp4, s->ar[i]);
                    break;
                default:
                    error("Not supported read number!");
                    break;
            }
        }

        for(int i = 0; i < tmp1->l; i++)
        {
            s->ar[i] = tmp1->ar[i];
        }
        for(int i = tmp1->l; i < tmp1->l + tmp2->l; i++)
        {
            s->ar[i] = tmp2->ar[i - tmp1->l];
        }
        for(int i = tmp1->l + tmp2->l; i < tmp1->l + tmp2->l + tmp3->l; i++ )
        {
            s->ar[i] = tmp3->ar[i - (tmp1->l + tmp2->l)];
        }
        for(int i = tmp1->l + tmp2->l + tmp3->l; i < tmp1->l + tmp2->l + tmp3->l + tmp4->l; i++)
        {
            s->ar[i] = tmp4->ar[i - ( tmp1->l + tmp2->l + tmp3->l)];
        }
        kaFree(tmp1);
        kaFree(tmp2);
        kaFree(tmp3);
        kaFree(tmp4);
        free(tmp1);
        free(tmp2);
        free(tmp3);
        free(tmp4);
    }
}

int parseBarcodeChainList(barcl_k * barcl, const kson_node_t * p)
{
    // unsorted barch 
    barch_k * tmp = (barch_k *) calloc(p->n, sizeof(barch_k));
    // temporary barcode chain order based on the GENOME TYPE and read position
    arwl_k * BarcodeChainOrderTMP = (arwl_k *) calloc(1, sizeof(arwl_k));
    BarcodeChainOrderTMP->ar = (size_t *) calloc(p->n, sizeof(size_t));
    int k;

    arrstr_k * AllGenomeTypeArray = (arrstr_k *) calloc(1, sizeof(arrstr_k));
    // process each barcode chain
    for(int i = 0; i < p->n; i++)
    {
        // flag show that if the barcode chain have genome (DNA/RNA) already.
        size_t GFlag = 1;
        tmp[i] = parseBarcodeChain(p->v.child[i]);
        // barcode chain i
        for(int j = 0; j < tmp[i].l; j++)
        {
            if(GFlag)
            {
                // get the genome type of current barcode chain
                if(tmp[i].barh[j].l == 1)
                {
                    arrstr_k * GenomeTypeArray = splitStrByStr(tmp[i].barh[j].bar[0].na, "GENOME");
                    if(GenomeTypeArray->l > 0 && strncmp(tmp[i].barh[j].bar[0].na, "GENOME", 6) == 0)
                    {
                        switch (tmp[i].rdp.rd)
                        {
                            case 1:
                                tmp[i].rw = 1;
                                break;
                            case 2:
                                // printf("Write 2 Ready\n");
                                tmp[i].rw = 2;
                                break;
                            case 3:
                                tmp[i].rw = 4;
                                break;
                            case 4:
                                tmp[i].rw = 8;
                                break;
                            default:
                                break;
                        }
                        arrstr_k * GenomeSubArray = splitStrByChar(GenomeTypeArray->mem + GenomeTypeArray->offset_ar[0], '|');
                        tmp[i].gtype = (GenomeSubArray->mem + GenomeSubArray->offset_ar[0])[0];
                        for(k = 0; k < AllGenomeTypeArray->l; k++)
                        {
                            if(strcmp(GenomeSubArray->mem + GenomeSubArray->offset_ar[0], AllGenomeTypeArray->mem + AllGenomeTypeArray->offset_ar[k]) == 0)
                            {                
                                BarcodeChainOrderTMP->ar[i] = k + 1;
                                GFlag = 0;
                                break;
                            }
                        }
                        if(k == AllGenomeTypeArray->l)
                        {
                            BarcodeChainOrderTMP->l++;
                            BarcodeChainOrderTMP->ar[i] = AllGenomeTypeArray->l + 1;
                            addStrToArrstrXC(AllGenomeTypeArray, GenomeSubArray->mem + GenomeSubArray->offset_ar[0]);
                            GFlag = 0;
                        }
                        arrstrFree(GenomeSubArray);
                        free(GenomeSubArray);
                    }
                    arrstrFree(GenomeTypeArray);
                    free(GenomeTypeArray);
                }
            }
        }
    }
    arrstrFree(AllGenomeTypeArray);
    free(AllGenomeTypeArray);


    arwl_k * BarcodeChainOrderArray = (arwl_k *) calloc(BarcodeChainOrderTMP->l + 1, sizeof(arwl_k));
    for(int i = 0; i < p->n; i++) addEleToArwlXC(&BarcodeChainOrderArray[BarcodeChainOrderTMP->ar[i]], i);
    // rearrange the order of barcode chain of GENOMEA GENOMEB GENOMEC separately by their read pos, right now the barcode chains are arranged based on the order GENOME: GENOMEA|2: R2:1, GENOMEA|1: R1:1, this function only rearrange inside the GENOME.
    for(int i = 0; i < BarcodeChainOrderTMP->l + 1; i++)
    {
        // This function is very dangerous cause the element of ar should be checked if exceeded barcl length.
        rearrangeArwlByReadPos(&BarcodeChainOrderArray[i], tmp);
    }

    
    size_t * FinalArray = (size_t *) calloc(p->n, sizeof(size_t));
    size_t count = 0, pos = 0;

    for(int i = 0; i < BarcodeChainOrderTMP->l + 1; i++)
    {
        // i == 0 : common barcode.
        if(i == 0)
        {
            for(int j = 0; j < BarcodeChainOrderArray[0].l; j++)
            {
                FinalArray[p->n - BarcodeChainOrderArray[0].l + j] = BarcodeChainOrderArray[0].ar[j];
                barcl->jump->ar[p->n - BarcodeChainOrderArray[0].l + j] = p->n - BarcodeChainOrderArray[0].l + j;
            }
        }
        // GENOME
        else
        {
            // printf("Length : %ld\n", BarcodeChainOrderArray[i].l);
            for(int k = 0; k < BarcodeChainOrderArray[i].l; k++)
            {
                barcl->jump->ar[pos + k] = pos + BarcodeChainOrderArray[i].l;
                FinalArray[count] = BarcodeChainOrderArray[i].ar[k];
                count++;
            }
            pos = pos + BarcodeChainOrderArray[i].l;
            // size_t InnerCount = 0, preRead = tmp[BarcodeChainOrderArray[i].ar[0]].rdp.rd;
            // for(j = 0; j < BarcodeChainOrderArray[i].l; j++)
            // {
            //     FinalArray[count] = BarcodeChainOrderArray[i].ar[j];
            //     if(preRead == tmp[BarcodeChainOrderArray[i].ar[j]].rdp.rd)
            //     {
            //         InnerCount++;
            //     }
            //     else
            //     {
            //         for(k = pos; k < pos + j; k++)
            //         {
            //             barcl->jump->ar[k] = OutCount + InnerCount;
            //         }
            //         pos = k++;
            //         preRead = tmp[BarcodeChainOrderArray[i].ar[j]].rdp.rd;
            //     }
            //     count++;
            // }
            // OutCount = OutCount + BarcodeChainOrderArray[i].l;
            // for(int k = pos; k < OutCount - BarcodeChainOrderArray[i].l + j; k++)
            // {
            //     barcl->jump->ar[k] = OutCount;
            // }
            // pos = k++;
        }
    }
    for(int j = 0; j < BarcodeChainOrderTMP->l+1; j++) kaFree(&BarcodeChainOrderArray[j]);
    kaFree(BarcodeChainOrderTMP);
    free(BarcodeChainOrderTMP);
    free(BarcodeChainOrderArray);
    for(int i = 0; i < p->n; i++) memmove(&barcl->barch[i], &tmp[FinalArray[i]], sizeof(barch_k));
    free(FinalArray);
    free(tmp);
    return 0;
}

int parseBarcodeHash(barch_k * barch, const kson_node_t * p)
{
    const kson_node_t * barp;
    for(int i = 0; i < barch->l; i++)
    {
        size_t Count = 0, Marker = 1, TmpMaxB = 0;
        arrstr_k * BarcodeNameArray = splitStrByStr(barch->barh[i].bar[0].na, "GENOME");
        // Barcode is Genome
        if(BarcodeNameArray->l > 0 && strncmp(barch->barh[i].bar[0].na, "GENOME", 6) == 0)
        {
            barch->barh[i].bar->ord = SIZE_MAX;
            barp = kson_by_key(p, barch->barh[i].bar[0].na);
            if(barp)
            {
                if(!kson_by_key(barp, "SPACE"))
                {
                    barch->barh[i].bar[0].sp = 0;
                }
                else
                {
                    if(atoi(kson_by_key(barp, "SPACE")->v.str) >= 0) barch->barh[i].bar[0].sp = atoi(kson_by_key(barp, "SPACE")->v.str); else barch->barh[i].bar[0].sp = 0;
                }
                if(!kson_by_key(barp, "LENGTH"))
                {
                    warnings("Genome \"Length\" %s not provided, will output sequence just after the barcode before it!", barch->barh[i].bar[0].na);
                    barch->barh[i].bar[0].minlen = barch->barh[i].bar[0].maxlen = SIZE_MAX;
                }
                else
                {
                    arrstr_k * bl = splitStrByChar(kson_by_key(barp, "LENGTH")->v.str, '-');
                    if(bl->l == 1)
                    {
                        if(atoi(bl->mem + bl->offset_ar[0]) >= 0)
                        {
                            barch->barh[i].bar[0].minlen = barch->barh[i].bar[0].maxlen = atoi(bl->mem + bl->offset_ar[0]);
                        }
                        else
                        {
                            barch->barh[i].bar[0].minlen = barch->barh[i].bar[0].maxlen = SIZE_MAX;
                        }
                    }
                    else if( bl->l == 2)
                    {
                        warnings("Genome %s \"LENGTH\" format is variable, will be set to the smaller value!", barch->barh[i].bar[0].na);
                        barch->barh[i].bar[0].minlen = barch->barh[i].bar[0].maxlen = Min(atoi(bl->mem + bl->offset_ar[0]), atoi(bl->mem + bl->offset_ar[1]));         
                    }
                    else
                    {
                        error("Genome %s \"LENGTH\" format is not correct, please check!", barch->barh[i].bar[0].na);
                    }
                    arrstrFree(bl);
                    free(bl);
                }
            }
            else
            {
                warnings("Genome %s not provided,  will output sequence just after the barcode before it!", barch->barh[i].bar[0].na);
                barch->barh[i].bar[0].sp = 0;
                barch->barh[i].bar[0].minlen = barch->barh[i].bar[0].maxlen = SIZE_MAX;
            }
            barch->barh[i].bar[0].hash = NULL;
            barch->barh[i].bar[0].rec = NULL;
            barch->barh[i].bar[0].la = 0;
            barch->barh[i].bar[0].mist = 0;
        }
        else
        {
            barp = kson_by_key(p, barch->barh[i].bar[0].na);
            if(!barp)
            {
                error("Barcode %s do not exist in the config file!", barch->barh[i].bar[0].na);
                return -1;
            }
            if(!kson_by_key(barp, "SPACE"))
            {
                warnings("\"SPACE\" do not exist in barcode %s, setting to default: 0.", barch->barh[i].bar[0].na);
                barch->barh[i].bar[0].sp = 0;
            }
            else
            {
                if(atoi(kson_by_key(barp, "SPACE")->v.str) >= 0)
                {
                    barch->barh[i].bar[0].sp = atoi(kson_by_key(barp, "SPACE")->v.str);
                }
                else
                {
                    warnings("\"SPACE\" in barcode %s provided a negative value, setting to default: 0.", barch->barh[i].bar[0].na);
                    barch->barh[i].bar[0].sp = 0;
                }
            }

            if(!kson_by_key(barp, "DENSE"))
            {
                warnings("\"DENSE\" do not exist in barcode %s, setting to default: SPARSE.", barch->barh[i].bar[0].na);
                barch->barh[i].bar[0].bighash = 0;
            }
            else
            {
                if(atoi(kson_by_key(barp, "DENSE")->v.str) >= 0)
                {
                    barch->barh[i].bar[0].bighash = atoi(kson_by_key(barp, "DENSE")->v.str);
                }
                else
                {
                    warnings("\"DENSE\" in barcode %s provided a negative value, setting to default: SPARSE.", barch->barh[i].bar[0].na);
                    barch->barh[i].bar[0].bighash = 0;
                }
            }

            if(!kson_by_key(barp, "RANDOM"))
            {
                warnings("\"RANDOM\" do not exist in barcode %s, setting to default: WHITE LIST needed.", barch->barh[i].bar[0].na);
                barch->barh[i].bar[0].random = 0;
            }
            else
            {
                if(atoi(kson_by_key(barp, "RANDOM")->v.str) >= 0)
                {
                    barch->barh[i].bar[0].random = atoi(kson_by_key(barp, "RANDOM")->v.str);
                }
                else
                {
                    warnings("\"RANDOM\" in barcode %s provided a negative value, setting to default: WHITE LIST needed.", barch->barh[i].bar[0].na);
                    barch->barh[i].bar[0].random = 0;
                }
            }

            if(!kson_by_key(barp, "LAXITY"))
            {
                warnings("\"LAXITY\" do not exist in barcode %s, setting to default: 0.", barch->barh[i].bar[0].na);
                barch->barh[i].bar[0].la = 0;
            }
            else
            {
                if(atoi(kson_by_key(barp, "LAXITY")->v.str) >= 0)
                {
                    barch->barh[i].bar[0].la = atoi(kson_by_key(barp, "LAXITY")->v.str);
                }
                else
                {
                    warnings("\"LAXITY\" in barcode %s provided a negative value, setting to default: 0.", barch->barh[i].bar[0].na);
                    barch->barh[i].bar[0].la = 0;
                }
            }

            if(kson_by_key(barp, "LENGTH") == NULL)
            {
                warnings("\"LENGTH\" do not exist in barcode %s, setting to default: 0.", barch->barh[i].bar[0].na);
                barch->barh[i].bar[0].maxlen = barch->barh[i].bar[0].minlen = 0;
            }
            else
            {
                arrstr_k * bl = splitStrByChar(kson_by_key(barp, "LENGTH")->v.str, '-');
                switch(bl->l)
                {
                    case 1: 
                        barch->barh[i].bar[0].minlen = barch->barh[i].bar[0].maxlen = atoi(bl->mem + bl->offset_ar[0]);
                        break;
                    case 2:
                        barch->barh[i].bar[0].minlen = Min(atoi(bl->mem + bl->offset_ar[0]), atoi(bl->mem + bl->offset_ar[1]));
                        barch->barh[i].bar[0].maxlen = Max(atoi(bl->mem + bl->offset_ar[0]), atoi(bl->mem + bl->offset_ar[1]));
                        break;
                    default:
                        warnings("Barcode %s \"LENGTH\" format is not supported, please check!", barch->barh[i].bar[0].na);
                }
                arrstrFree(bl);
                free(bl);
            }

            if(kson_by_key(barp, "MISMATCH") == NULL)
            {
                warnings("\"MISMATCH\" do not exist in barcode %s, setting to default: 0.", barch->barh[i].bar[0].na);
                barch->barh[i].bar[0].mist = 0;
            }
            else
            {
                if(atoi(kson_by_key(barp, "MISMATCH")->v.str) >= 0)
                {
                    barch->barh[i].bar[0].mist = atoi(kson_by_key(barp, "MISMATCH")->v.str);
                }
                else
                {
                    barch->barh[i].bar[0].mist = 0;
                }
            }

            // harh_k * har = initParseHashTable(barch->barh[i].bar[0].mist);
            kh_m8_t * har = kh_init(m8);
            const kson_node_t * wl = kson_by_key(barp, "WHITE LIST");
            // whitelsit did not provided, either construct hash table by random or no hash table at all. 
            if(barch->barh[i].bar[0].random == 0)
            {
                if(wl == NULL || barch->barh[i].bar[0].minlen == 0)
                {
                    warnings("%s do not have white list or \"LENGTH\" set to 0, will not be idenfied!", barch->barh[i].bar[0].na);
                    hashDestroyB(har);
                    // free(har);
                    har = NULL;
                    arrstrFree(BarcodeNameArray);
                    free(BarcodeNameArray);
                    continue;
                }
                // whitelist provided, construct hash table based on whitelist.
                else
                {
                    arrstr_k * TmpArray = (arrstr_k *) calloc(1, sizeof(arrstr_k));
                    for( int j = 0; j < wl->n; j++)
                    {
                        char * seq, * name;
                        switch(wl->v.child[j]->n)
                        {
                            case 0:
                            {
                                seq = wl->v.child[j]->v.str;
                                name = strConcate(barch->barh[i].bar[0].na, seq, ":");
                                goto ADDELETOHASH;
                            }
                            case 1:
                            {
                                seq = wl->v.child[j]->v.child[0]->v.str;
                                name = (char *) calloc(strlen(wl->v.child[j]->v.child[0]->key) + 1, sizeof(char));
                                memmove(name, wl->v.child[j]->v.child[0]->key, strlen(wl->v.child[j]->v.child[0]->key) * sizeof(char));
                                name[strlen(wl->v.child[j]->v.child[0]->key)] = '\0';
                                goto ADDELETOHASH;
                            }
                            default:
                            {
                                error("\"WHITE LIST\" format do not supported at %s!\n", barch->barh[i].bar[0].na);
                                goto END;
                            }
                        }
                        ADDELETOHASH:
                        {
                            
                            arrstr_k * Tmp = splitStrByChar(seq, 'N');
                            if(Tmp->l == 0)
                            {
                                hashDestroyB(har);
                                har = NULL;
                                arrstrFree(TmpArray);
                                free(TmpArray);
                                TmpArray = NULL;
                                free(name);
                                name = NULL;
                                arrstrFree(Tmp);
                                free(Tmp);
                                break;
                            }
                            else
                            {
                                if(Count > SIZE_MAX)
                                {
                                    error("Barcode %s number exceeded limits: %zu", barch->barh[i].bar[0].na, SIZE_MAX);
                                }
                                else
                                {
                                    arwl_k * Conflict;
                                    if(barch->barh[i].bar[0].bighash)
                                    {
                                        Conflict = addSeqToHashTableXCB(seq, har, Count, 0);
                                    }
                                    else
                                    {
                                        Conflict = addSeqToHashTableXCB(seq, har, Count, barch->barh[i].bar[0].mist);
                                    }
                                    Count++;
                                    if(Count > Marker) 
                                    {
                                        Marker = Marker << 1;
                                        TmpMaxB++;
                                    }
                                    addStrToArrstrRepXC(TmpArray, name);
                                    char * OldflictSeq = (char *) calloc(1, sizeof(char));
                                    for(int i = 0; i < Conflict->l; i++)
                                    {
                                        char * NewConflictSeq = strConcateByChar(OldflictSeq, TmpArray->mem + TmpArray->offset_ar[Conflict->ar[i]], ' ');
                                        free(OldflictSeq);
                                        OldflictSeq = NewConflictSeq;
                                    }
                                    if(Conflict->l > 1) warnings("%s have conflict with: %s, result will output as: %s\n", name, OldflictSeq, TmpArray->mem + TmpArray->offset_ar[Conflict->ar[0]]); 
                                    free(OldflictSeq);
                                    kaFree(Conflict);
                                    free(Conflict);
                                    arrstrFree(Tmp);
                                    free(Tmp);
                                    continue;
                                }
                            }
                        }
                        END: break;
                    }
                    barch->barh[i].bar[0].rec = TmpArray;
                }
            }
            else
            {
                size_t RANDOM_B = (size_t) 64/barch->l - 1;
                if(barch->maxb < RANDOM_B) barch->maxb = RANDOM_B;
                barch->barh[i].bar[0].rec = (arrstr_k *) calloc(1, sizeof(arrstr_k));
            }
            barch->barh[i].bar[0].hash = har;
        }
        if(TmpMaxB > barch->maxb) barch->maxb = TmpMaxB;
        arrstrFree(BarcodeNameArray);
        free(BarcodeNameArray);
        printf( ANSI_COLOR_CYAN "Added %ld Barcode of %s." ANSI_COLOR_RESET "\n", Count, barch->barh[i].bar[0].na);  
    }
    return 0;
}       

id_k parseBarcodeIdentifier(const kson_node_t * p, arrstr_k * strar)
{
    id_k ret = {0};
    ret.na = strCopy(p->v.child[0]->key);

    arrstr_k * IdentifierArray = splitStrByChar(p->v.child[0]->v.str, '-');
    // mark the barcode need or not to be find in the barcode chain -1 no; 1 yes
    arwl_k * Tmp = (arwl_k *) calloc(1, sizeof(arwl_k));
    // barcode hole array;
    arwl_k * Var = (arwl_k *) calloc(IdentifierArray->l, sizeof(arwl_k));
    // printf("Identifier Length : %ld\n", IdentifierArray->l);
    size_t curpos = 1;

    for(int i = 0; i < IdentifierArray->l; i++)
    {
        char * curbar = IdentifierArray->mem + IdentifierArray->offset_ar[i];
        arrstr_k * cbar = splitStrByChar(curbar, ':');
        char * barna;
        switch(cbar->l)
        {
            case 1:
            {
                barna = cbar->mem + cbar->offset_ar[0];
                goto PARSE;
            }
            case 2:
            {
                if(atoi(cbar->mem + cbar->offset_ar[0]) > 0 && atoi(cbar->mem + cbar->offset_ar[0]) > curpos)
                {
                    int MoveLen = (atoi(cbar->mem + cbar->offset_ar[0]) - curpos);
                    ret.anchor = 1;
                    for(int k = 0; k < MoveLen; k++)
                    {
                        addEleToArwlRep(Tmp, SIZE_MAX);
                        curpos++;
                    }
                }
                else
                {
                    error("Identifier format error at %s:%s! Barcode position should be assigned in order!\n",  p->v.child[0]->v.str, curbar);
                }
                barna = cbar->mem + cbar->offset_ar[1];
                goto PARSE;
            }
            default:
            {
                error("Identifier format error at %s:%s!\n", p->v.child[0]->v.str, curbar);
                break;
            }
            PARSE:
            {
                if(strcmp(barna, "ANY") == 0)
                {
                    addEleToArwlRep(Tmp, SIZE_MAX);
                    curpos++;
                }
                else
                {
                    addEleToArwlRep(Tmp, 1);
                    arrstr_k * Hole = splitStrByChar(barna, '/');
                    for(int k = 0; k < Hole->l; k++)
                    {
                        for(int o = 0; o < strar->l; o++)
                        {
                            // o + 1 is the barcode ord.
                            if(strcmp(Hole->mem + Hole->offset_ar[k], strar->mem + strar->offset_ar[o]) == 0)
                            {
                                // printf("Var[%d] : %ld\n", i, (size_t) o + 1);
                                addEleToArwlRep(&Var[i], (size_t) o + 1);
                            }
                        }
                    }
                    curpos++;
                    arrstrFree(Hole);
                    free(Hole);
                }
            }
        }
        arrstrFree(cbar);
        free(cbar);
    }
    
    ret.l = Tmp->l;
    size_t q = 0;
    // identifier q; 10111; do not need the 2 one to be identified.
    for(int i = 0; i < Tmp->l; i++)
    {
        // printf("Tmp[%d] : %ld\n", i, Tmp->ar[i]);
        if(Tmp->ar[i] != SIZE_MAX) q = (q << 1) | 1; else q = (q << 1) | 0;
    }

    ret.q = q;
    kaFree(Tmp);
    free(Tmp);

    // number of combination of barcode.
    size_t vcount = 1;
    for(int i = 0; i < IdentifierArray->l; i++)
    {
        if(Var[i].l != 0) vcount = vcount * Var[i].l; 
    }

    ret.v = (arwl_k *) calloc(1, sizeof(arwl_k));
    ret.v->l = vcount;
    ret.v->ar = (size_t *) calloc(ret.v->l, sizeof(size_t));

    size_t count = 0, split = ret.v->l;
    arwl_k * Tmp2 = (arwl_k *) calloc(ret.v->l, sizeof(arwl_k));
    for(int i = 0; i < ret.v->l; i++)
    {
        Tmp2[i].ar = (size_t *) calloc(ret.l, sizeof(size_t));
        Tmp2[i].l = ret.l;
    }

    for(int i = 0; i < IdentifierArray->l; i++)
    {
        count = 0;
        if(Var[i].l > 0)
        {
            split = (size_t) split/Var[i].l;
            while(count < ret.v->l)
            {
                int k = 0;
                for(int j = 0; j < (size_t) ret.v->l/split; j++)
                {
                    for(int o = split * j; o < split * (j + 1); o++)
                    {
                        // printf("Var[%d].ar[%d] : %ld\n", i, k, Var[i].ar[k]);
                        Tmp2[o].ar[i] = Var[i].ar[k];
                        // printf("Tmp2[%d].ar[%d] : %ld\n", o, i, Tmp2[o].ar[i]);
                        count++;
                    }
                    k++;
                    if(k == Var[i].l) k = 0;
                }
            }
        }
    }

    for(int i = 0; i < ret.v->l; i++)
    {
        size_t vq = ret.q, tmpv = 0, flag = 1 << (ret.l - 1), Counter = 0;
        for(int j = 0; j < ret.l; j++)
        {
            // printf("vq : %ld - flag : %ld\n", vq, flag);
            if((vq & flag) > 0)
            {
                // printf("Tmp2[%d].ar[%ld] : %ld\n",i, Counter, Tmp2[i].ar[Counter]);
                tmpv  = (tmpv << BUNIT) | Tmp2[i].ar[Counter];
                Counter++;
            }
            flag = flag >> 1;
        }
        ret.v->ar[i] = tmpv;
        // printf("ret.v : %ld\n", ret.v->ar[i]);
    }

    for(int i = 0; i < ret.v->l; i++) kaFree(&Tmp2[i]);
    free(Tmp2);

    for(int i = 0; i < IdentifierArray->l; i++)  kaFree(&Var[i]);
    free(Var);

    arrstrFree(IdentifierArray);
    free(IdentifierArray);
    return ret;
}

int parseBarcodeIdentifierList(arrstr_k * StrArray, const kson_node_t * p, barcl_k * barcl)
{
    barcl->idf = (id_k *) calloc(p->n, sizeof(id_k));
    barcl->idl = p->n;
    for(int i = 0; i < p->n; i++)
    {
        barcl->idf[i] = parseBarcodeIdentifier(p->v.child[i], StrArray);
    }
    return 0;
}

arrstr_k parseSon(const kson_node_t * p)
{
    arrstr_k ret = {0};
    ret.l = 2;
    if(p->v.child[0]->n == 0)
    {
        addStrToArrstrRepXC(&ret, p->v.child[0]->key);
        addStrToArrstrRepXC(&ret, p->v.child[0]->v.str);
    }
    else
    {
        error("SCSM only support 1 level of son!\n");
    }
    return ret;
}

int parseSonList(const kson_node_t * p, barcl_k * barcl)
{
    barcl->sonch = (arrstr_k *) calloc(p->n, sizeof(arrstr_k));
    barcl->sonl = p->n;
    for(int i = 0; i < p->n; i++)
    {
        barcl->sonch[i] = parseSon(p->v.child[i]);
    }
    return 0;
}

barcl_k * parseConfig(char * fname)
{
    barcl_k * barcl = (barcl_k *) calloc(1, sizeof(barcl_k));
    char * ConfigStr = readJson(fname);
    if(!ConfigStr) 
    {
        error("There are some error occur in config file! Please check!\n");
    }

    kson_t * ConfigJson = kson_parse(ConfigStr);
    free(ConfigStr); 

    const kson_node_t * p;
    int FlagChain = -1, FlagHash = -1 , FlagIden = -1, FlagSon = -1;
    // collect options not supported occur in config file.
    arrstr_k * NonOpt = (arrstr_k *) calloc(1, sizeof(arrstr_k));
    char * NonOptStr = (char *) calloc(1, sizeof(char));

    for(long i = 0; i < ConfigJson->root->n; i++)
    {
        p = ConfigJson->root->v.child[i];
        switch(optConvert(p->key))
        {
            case BARCODECHAIN:
            {
                if(FlagChain != 1) goto PARSECHAIN; else continue;
            }
            case BARCODETYPE:
            {
                if(FlagHash != 1 && FlagChain == 1) goto CREATEHASH; else continue;
            }
            case IDENTIFIER:
            {
                if(FlagIden != 1 && FlagHash == 1) goto PARSEIDENTIFY; else continue;
            }
            case SON:
            {
                if(FlagSon != 1 && FlagIden == 1) goto PARSESON; else continue;
            }
            default:
            {
                addStrToArrstrXC(NonOpt, p->key);
                continue;
            }
        }

        PARSECHAIN:
        {
            // fprintf(stdout, "Parsing barcode chain\n");
            barcl->l = p->n;
            barcl->jump = (arwl_k *) calloc(1, sizeof(arwl_k));
            barcl->jump->ar = (size_t *) calloc(p->n, p->n * sizeof(size_t));
            barcl->jump->l = p->n;
            barcl->barch = (barch_k *) calloc(p->n, sizeof(barch_k));
            if(!parseBarcodeChainList(barcl, p)) fprintf(stderr, ANSI_COLOR_LIGHTBLACK "Finished parsing barcode chain list!" ANSI_COLOR_RESET "\n");
            FlagChain = 1;
            i = -1;
            continue;
        }
        CREATEHASH:
        {
            arrstr_k * TmpStrArray = (arrstr_k *) calloc(1, sizeof(arrstr_k));
            for(int j = 0; j < barcl->l; j++)
            {
                for(int k = 0; k < barcl->barch[j].l; k++)
                {
                    for(int m = 0; m < barcl->barch[j].barh[k].l; m++) addStrToArrstrXC(TmpStrArray, barcl->barch[j].barh[k].bar[m].na);
                }
            }
            barch_k * BarclHash = (barch_k *) calloc(1, sizeof(barch_k));
            BarclHash->l = TmpStrArray->l;
            BarclHash->barh = (barh_k *) calloc(BarclHash->l, sizeof(barh_k));
            for(int j = 0; j < TmpStrArray->l; j++)
            {
                BarclHash->barh[j].l = 1;
                BarclHash->barh[j].bar = (bar_k *) calloc(1, sizeof(bar_k));
                BarclHash->barh[j].bar[0].na = strCopy(TmpStrArray->mem + TmpStrArray->offset_ar[j]);
                
            }
            parseBarcodeHash(BarclHash, p);
            for(int j = 0; j < barcl->l; j++)
            {
                for(int k = 0; k < barcl->barch[j].l; k++)
                {
                    for(int m = 0; m < barcl->barch[j].barh[k].l; m++)
                    {
                        for(int o = 0; o < BarclHash->l; o++)
                        {
                            if(strcmp(barcl->barch[j].barh[k].bar[m].na, BarclHash->barh[o].bar[0].na) == 0)
                            {
                                memmove(&barcl->barch[j].barh[k].bar[m], &BarclHash->barh[o].bar[0], sizeof(bar_k));
                                break;
                            }
                        }
                    }
                }
                barcl->barch[j].maxb = BarclHash->maxb;
            }
            arrstrFree(TmpStrArray);
            free(TmpStrArray);
            fprintf(stderr,ANSI_COLOR_LIGHTBLACK "Finished parsing barcode white list!" ANSI_COLOR_RESET "\n");
            FlagHash = 1;
            i = -1;
            continue;
        }
        PARSEIDENTIFY:
        {
            arrstr_k * TmpArray = (arrstr_k *) calloc(1, sizeof(arrstr_k));
            // get unique barcode
            for(int j = 0; j < barcl->l; j++)
            {
                for(int k = 0; k < barcl->barch[j].l; k++)
                {
                    for(int m = 0; m < barcl->barch[j].barh[k].l; m++)
                    {
                        arrstr_k * Array = splitStrByStr(barcl->barch[j].barh[k].bar[m].na, "GENOME");
                        if(Array->l > 0 && strncmp(barcl->barch[j].barh[k].bar[m].na, "GENOME", 6) == 0)
                        {
                            arrstrFree(Array);
                            free(Array);
                            continue;
                        }
                        else
                        {
                            addStrToArrstrXC(TmpArray, barcl->barch[j].barh[k].bar[m].na);
                        }
                        arrstrFree(Array);
                        free(Array);
                    }
                }
            }
            // assign a unique code for each barcode
            for(int j = 0; j < barcl->l; j++)
            {
                for(int k = 0; k < barcl->barch[j].l; k++)
                {
                    for(int m = 0; m < barcl->barch[j].barh[k].l; m++)
                    {
                        for(int o = 0; o < TmpArray->l; o++)
                        {
                            if(strcmp(barcl->barch[j].barh[k].bar[m].na, TmpArray->mem + TmpArray->offset_ar[o]) == 0)
                            {
                                barcl->barch[j].barh[k].bar[m].ord = o + 1;
                                break;
                            }
                        }
                    }
                }
            }
            parseBarcodeIdentifierList(TmpArray, p, barcl);
            arrstrFree(TmpArray);
            free(TmpArray);
            fprintf(stderr,ANSI_COLOR_LIGHTBLACK "Finished parsing identifier!"ANSI_COLOR_RESET "\n");
            FlagIden = 1;
            i = -1;
            continue;
        }

        PARSESON:
        {
            parseSonList(p, barcl);
            fprintf(stderr, ANSI_COLOR_LIGHTBLACK "Finished parsing Son!" "\n");
            continue;
        }
    }
    for(int j = 0; j < NonOpt->l; j++)
    {
        char * OldNonOptStr = NonOptStr;
        if(j == 0) strConcate(NonOptStr, NonOpt->mem + NonOpt->offset_ar[j], ""); else strConcate(NonOptStr, NonOpt->mem + NonOpt->offset_ar[j], ", ");
        free(OldNonOptStr);
    }
    if(NonOpt->l > 0)
    {
        if(NonOpt->l > 1) 
        {
            warnings("Options %s are not supported yet, ignoring...\n", NonOptStr);
        }
        else 
        {
            warnings("Option %s is not supported yet, ignoring...\n", NonOptStr);
        }
    }
    free(NonOptStr);
    kson_destroy(ConfigJson);
    arrstrFree(NonOpt);
    free(NonOpt);
    NonOpt = NULL;
    return barcl;
}
