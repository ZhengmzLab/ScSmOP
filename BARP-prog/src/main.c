#include <stdio.h>
#include <string.h>
#include "kstring.h"
#include "kutils.h"

#ifndef PACKAGE_VERSION 
#define PACKAGE_VERSION "0.0.1"
#endif


int parseBarcode(int argc, char * argv[]);
int parseToSamField(int argc, char * argv[]);
// int removeDuplicate(int argc, char * argv[]);
int mergeFragment(int argc, char * argv[]);
// int ConVert(int argc, char * argv[]);


static int Usage()
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Program: barp (barcode identification and cluster via space seed)\n");
    fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
    fprintf(stderr, "Contact: Jing Kai <12032179@mail.sustech.edu.cn>\n\n");
    fprintf(stderr, "Usage:   barp <command> [options]\n\n");
    fprintf(stderr, "Command: idb\t\tidentify barcode and group according to config\n");
    fprintf(stderr, "         p2s\t\tparse information from/to sam field\n");
    fprintf(stderr, "         merge\t\tmerge fragments for inline BED file\n");
    fprintf(stderr, "\n");
    return 1;
}


int main(int argc, char * argv[])
{
    int i, ret;
    double t_real;
    t_real = realtime();
    if (argc < 2) return Usage();
    if (strcmp(argv[1], "idb") == 0) ret = parseBarcode(argc-1, argv+1);
    else if (strcmp(argv[1], "p2s") == 0) ret = parseToSamField(argc-1, argv+1);
    // else if (strcmp(argv[1], "rdp") == 0) ret = removeDuplicate(argc-1, argv+1);
    else if (strcmp(argv[1], "merge") == 0) ret = mergeFragment(argc-1, argv+1);
    // else if (strcmp(argv[1], "convert") == 0) ret = ConVert(argc-1, argv+1);
    else {
        fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
        return 1;
	}
    if (ret == 0) {
		fprintf(stderr, "[%s] Version: %s\n", __func__, PACKAGE_VERSION);
		fprintf(stderr, "[%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - t_real, cputime());
	}
    return ret;
}