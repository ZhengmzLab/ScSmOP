#include "kutils.h"
#include "config_hash.h"
#include "unistd.h"
#include "parse_config.h"
#include <getopt.h>

static int ConvertUsage()
{
    fprintf(stdout, "\n");
    fprintf(stdout, "Usage: BARP convert -b [STR] -c [Config file] -f [Field]\n"
    "Options:\n"
    "  -b           STR     Identifier string.\n"
    "  -c           FILE    Config file include barcode chain and white list\n"
    "  -f           INT     TenX barcode field.\n"
    );
    fprintf(stdout, "\n");
    return -1;
}

static inline arrstr_k * parseIdenFile(char * idfstr)
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

int ConVert(int argc, char * argv[])
{

    char * IDFStr = "-", * Config = "-", * Name = "BC";
    int BarcodeColumn = 0;
    char c;
    int nargs;
    int opidx = 0;
    // ConvertTo: 1 : TenX to BARP; ConvertTo: 2 : BARP to TenX.
    int ConvertTo = 1;
    static const struct option lopts[] = {
        {NULL, 0, NULL, 0}
    };
    while ((c = getopt_long(argc, argv, "b:f:c:", lopts, &opidx)) >= 0)
    {
        if(c == -1) break;

        switch(c)
        {
            case 'f' : BarcodeColumn = atoi(optarg); break;
            case 'c' : Config = optarg; break;
            case 'b' : IDFStr = optarg; break;
            case 'n' : Name = optarg; break;
            case 'v' : 
            case '?' : ConvertUsage(); return -1;
        }
    }
    nargs = argc - optind;
    if(nargs == 0)
    {
        ConvertUsage();
        return 1;
    }

    if(strcmp(IDFStr, "-") == 0)
    {
        error("No identifier string provided.");
    }

    if(strcmp(Config, "-") == 0)
    {
        error("No config file found.");
    }

    if(BarcodeColumn < 0) error("Barcode column error\n");


    char * fi = argv[optind];
    char * fo = strConcateByChar(fi, "BARPID.txt", '_');

    FILE * CHDPFile = fopen(fi, "r");
    FILE * FOUT = fopen(fo, "w");
    


    initMistComArray();
    initPrekQArray();
    barcl_k * barcl = parseConfig(Config);

    kh_m16_t ** hart = (kh_m16_t **) calloc(barcl->idl, sizeof(kh_m16_t *));
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
    arrstr_k * IdenFileStrArray = parseIdenFile(IDFStr);
    int IdenNum = -1;
    for(int i = 0; i < barcl->idl; i++)
    {
        for(int j = 0; j < IdenFileStrArray->l; j++)
        {
            if(strcmp(barcl->idf[i].na, IdenFileStrArray->mem + IdenFileStrArray->offset_ar[j]) == 0)
            {
                fprintf(stdout, ANSI_COLOR_CYAN "Read identifier %s from %s.", barcl->idf[i].na, IdenFileStrArray->mem + IdenFileStrArray->offset_ar[j + 1]);
                fprintf(stdout, ANSI_COLOR_RESET "\n");
                IdenNum = i;
                if(!(idfcount->ar[i] = pushIdentifierFromFile(hart[i], IdenFileStrArray->mem + IdenFileStrArray->offset_ar[j + 1]))) error("Read Identifier from file corrupted\n");
            }
        }
    }
    arrstrFree(IdenFileStrArray);
    free(IdenFileStrArray);

    int i, j;
    int Chain = -1, Hole = 0, Bar = 0;
    for( i = 0; i < barcl->l; i++)
    {
        fprintf(stdout, ANSI_COLOR_GREEN "[Barcode Chain] : " ANSI_COLOR_CYAN);
        for(j = 0; j < barcl->barch[i].l; j++)
        {
            if(j == 0) fprintf(stdout, "%s", barcl->barch[i].barh[j].bar[0].na); else fprintf(stdout, " - %s", barcl->barch[i].barh[j].bar[0].na);
            if(strcmp(barcl->barch[i].barh[j].bar[0].na, Name) == 0)
            {
                Chain = i;
                Hole = j;
                Bar = 0;
            }
            for(int m = 1; m < barcl->barch[i].barh[j].l; m++)
            {
                fprintf(stdout, " or %s", barcl->barch[i].barh[j].bar[m].na);
                if(strcmp(barcl->barch[i].barh[j].bar[m].na, Name) == 0)
                {
                    Chain = i;
                    Hole = j;
                    Bar = m;
                }
            }
        }
        fprintf(stdout, ANSI_COLOR_RESET "\n");
    }

    if(Chain == -1) error("Not barcode found in the config file.");

    kstring_t line = {0, 0, NULL};
    int Count = 0;
    while (kgetline(&line, (kgets_func *) fgets, (FILE *) CHDPFile) != EOF)
    {
        arrstr_k * Content = splitStrByChar(line.s, '\t');
        char * Barcode = Content->mem + Content->offset_ar[BarcodeColumn];
        bandrec_k QueryRet;
        QueryRet = querySeqBestC(Barcode, Barcode, barcl->barch[Chain].barh[Hole].bar[Bar].hash, 0, -1, 0, 0);
        size_t RecResult = QueryRet.rec;
        idvp_k v = {0, 0};
        v.cv = RecResult;
        v.sv = 0;

        size_t COMPLEX = queryIdentifier(v, hart[0], &idfcount->ar[i], 1);

        for(int i = 0; i < BarcodeColumn; i++)
        {
            fprintf(FOUT, "%s\t", Content->mem + Content->offset_ar[i]);
        }
        fprintf(FOUT, "%s_%ld", barcl->idf[IdenNum].na, COMPLEX);
        for(int i = 0; i < Content->l - BarcodeColumn - 1; i++)
        {
            fprintf(FOUT, "\t%s", Content->mem + Content->offset_ar[i]);
        }
        fprintf(FOUT, "\n");

        if(Count%10000 == 0 && Count > 0) processeinfo("%d.", Count);
        Count++;
        arrstrFree(Content);
        free(Content);
        ks_free(&line);
    }
    
    fclose(CHDPFile);
    fclose(FOUT);
}
