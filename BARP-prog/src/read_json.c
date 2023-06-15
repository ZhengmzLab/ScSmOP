#include <stdio.h>
#include "kstring.h"
#include <ctype.h>
#include <errno.h>


#define MAXSTACKSIZE 20

// Modify the line to valuable json context; if will be added to json string, return 1; else 0
int parseTxt(kstring_t * s)
{
    if(s->l  == 0)
    {
        return 0;
    }

    char * ss = s->s;
    char * se = s->s + s->l - 1;
    while(ss < se && isspace(*ss)) ss++;
    // skip comment line start with # or //
    if(* ss == '#') return 0;
    if(*ss == '/' && *(ss+1) == '/') return 0;
    while(se > ss && isspace(*se)) se--;

    if(ss != s->s || se != s->s + s->l - 1)
    {
        s->l = se - ss + 1;
        memmove(s->s, ss, s->l);
        s->s[s->l] = '\0';
    }

    return 1;
}

// Check the if the json format have paired {} or []
int checkComplete(kstring_t * s)
{
    size_t i, line = 1;

    char stack[MAXSTACKSIZE];
    int top = -1;

    for(i = 0; i < s->l; i++)
    {
        char curc = s->s[i];
        switch(curc)
        {
            case '\n':
                line++;
                goto END;
            case '{':
            case '[':
                goto StackPush;
            case '"':
            case '\'':
                if(curc != stack[top]) goto StackPush; else goto StackPop;
            case '}':
            case ']':
                if(curc == stack[top] + 2) goto StackPop; else goto PrintError;
            default:
                goto END;
        }
        StackPush:
            if(top < MAXSTACKSIZE) stack[++top] = curc; else fprintf(stderr, "Level of JSON exceeded\n!");
            continue;
        StackPop:
            if(top >= 0) top--; else return 0;
            continue;
        PrintError:
            fprintf(stderr, "Line %ld enclose: %c is not closed!", line, curc);
            return 1;
        END:
            continue;
    }

    return 0;
}

// read json file return string 
char * readJson(const char * fname)
{
    char * ret;
    FILE * fi;
    fi = fopen(fname, "r");

    if(!fi)
    {
        fprintf(stderr, "%s : %s\n", fname, strerror(errno));
        return NULL;
    }

    kstring_t * JsonString = (kstring_t *) calloc(1, sizeof(kstring_t));
    kstring_t * Tmp = (kstring_t *) calloc(1, sizeof(kstring_t));

    while(kgetline(Tmp, (kgets_func *) fgets, fi) >= 0){
        if(parseTxt(Tmp))
        {
            kputs(Tmp->s, JsonString);
            kputc('\n', JsonString);
        }
        Tmp->l = 0;
    }

    fclose(fi);
    if(Tmp->m) ks_free(Tmp);
    free(Tmp);
    if(checkComplete(JsonString)) return NULL;

    ret = JsonString->s;
    free(JsonString);
    return ret;
}

// int main(int argc, char const *argv[])
// {
//     char * fname = "../test/config_file/bgi_config.json";
//     while (1)
//     {
//         char * ConfigStr = readJson(fname);
//         printf("%s\n", ConfigStr);
//         free(ConfigStr);
//     }
//     return 0;
// }

