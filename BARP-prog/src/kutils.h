#ifndef _KUTILS_H_
#define _KUTILS_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <fcntl.h>
#include <sys/time.h>
#include <stdint.h>

#define CPUBYTENUM  sizeof(int)
#define MAX_SIZE 512
#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_GRAY    "\x1b[37m"
#define ANSI_COLOR_RESET   "\x1b[0m"
#define ANSI_COLOR_LIGHTBLACK "\e[1;30m"

#define error(line, ...) do \
{ \
    fprintf(stderr, ANSI_COLOR_RED "[error]  " ANSI_COLOR_RESET ANSI_COLOR_MAGENTA line ANSI_COLOR_RESET "\n", ##__VA_ARGS__);\
    exit(1); \
} while (0)

#define warnings(line, ...) do \
{ \
    fprintf(stderr, ANSI_COLOR_YELLOW "[warnings]  " line ANSI_COLOR_RESET "\n", ##__VA_ARGS__);\
} while (0)

#define processeinfo(line, ...) do \
{ \
    fprintf(stderr , ANSI_COLOR_LIGHTBLACK "[processed]  " line ANSI_COLOR_RESET "\n", ##__VA_ARGS__); \
} while (0)

// copied from bwa utils.h
static inline double cputime(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

static inline double realtime(void)
{
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}


//size_t array with length
typedef struct 
{
    size_t l;
    size_t * ar;
} arwl_k;
// short array
typedef struct arrstr_k
{
    size_t l, m;
    size_t * offset_ar;
    char * mem;
} arrstr_k; 

// pstring and key pair key used as the key for hash table
typedef struct prek_k
{
    uint64_t v;
    size_t q;
} prek_k;
// pk array all information: qar:[1, n, n(n-1), n(n-1)(n-mist)]; qlar: n, n-1, n-1, n-1 ...; ar: pk pair array; 
typedef struct pkarwl_k
{
    size_t l, * qlar;
    arwl_k * qar;
    prek_k * ar;
} pkarwl_k;

// one mismatch array; 
typedef struct omar_k
{
    size_t ql, l;
    prek_k * ar;
} omar_k;

typedef struct idfvpair
{
    // uint64_t cv, ov;
    uint64_t cv;
    size_t sv;
}idvp_k;

typedef struct idfvcountpair
{
    size_t od, count;
}idvcp_k;

typedef struct numberstack_k
{
    size_t element[MAX_SIZE];
    int top;
    int length;
}sstack_k;

typedef struct PthreadState
{
    pthread_t * ptar;
    size_t l;
    size_t * statear;
} pthrs_k;

static inline void initThread(pthrs_k * s, size_t n)
{
    s->l = n;
    s->ptar = (pthread_t *) calloc(n, sizeof(pthread_t));
    s->statear = (size_t *) calloc(n, sizeof(size_t));
    for(int i = 0; i < n; i++)
    {
        s->statear[i] = 1;
    }
}

static inline void freeThread(pthrs_k * s)
{
    if(s->ptar) free(s->ptar);
    if(s->statear) free(s->statear);
}
typedef struct ndimensionarray
{
    arwl_k * subarray;
    size_t l;
    size_t d;
    size_t top;
}ndarr_k;

typedef struct BarCounter
{
    size_t l;
    uint64_t * bar;
}barc_k;

// static inline void initNdArray(ndarr_k * s, size_t dimension)
// {
//     s->l = 0;
//     s->d = dimension;
//     s->top = -1;
// }

// unfinished. used to add element to ndarr
// static inline void addArrayToNdArray(ndarr_k * s, size_t num)
// {
//     if(s->top == -1)
//     {
//         s->subarray = (arwl_k *) calloc(1, sizeof(arwl_k));
//         s->subarray->ar = (size_t *) calloc(s->d, sizeof(size_t));
//         s->top = 0;
//         s->subarray[s->top].ar[0] = num;
//         s->l = 1;
//     }
//     else
//     {
//         if(s->subarray[s->top].l + 1 == s->d)
//         {
//             s->subarray[s->top].ar[s->subarray[s->top].l] = num;
//         }
//     }
// }

static inline void initNumStack(sstack_k * Stack)
{
    Stack->top = -1;
    Stack->length = 0;
}

static inline void arrstrFree(arrstr_k * s)
{
    if(s->mem) free(s->mem);
    if(s->offset_ar) free(s->offset_ar);
    s->mem = NULL;
    s->offset_ar = NULL;
    s->l = 0;
    s->m = 0;
}

static inline int pushNumStack(sstack_k * Stack, size_t element)
{
    if(Stack->top == MAX_SIZE - 1)
    {
        return 0;
    }
    Stack->top++;
    Stack->element[Stack->top] = element;
    Stack->length++;
    return 1;
}

static inline int popNumStack(sstack_k * Stack, size_t * element)
{
    if(Stack->top == -1)
    {
        return 0;
    }
    * element = Stack->element[Stack->top];
    Stack->top--;
    Stack->length--;
    return 1;
}

static inline void kaFree(arwl_k * s)
{
    if(s)
    {
        if(s->ar)
        {
            free(s->ar);
            s->ar = NULL;
        }
        s->l = 0;
    }
}


static inline void overLayArwl(arwl_k * s, arwl_k * m)
{
    if(s && m)
    {
        s->ar = (size_t *) realloc(s->ar, (s->l + m->l) * sizeof(size_t));
        memmove(s->ar + s->l, m->ar, m->l * sizeof(size_t));
        s->l = s->l + m->l;
        kaFree(m);
        free(m);
        m = NULL;
    }
    else if (!s && m)
    {
        if(s->ar) free(s->ar);
        s->ar = m->ar;
        s->l = m->l;
        free(m);
    }
}

static inline void overLayArrstr(arrstr_k * s, arrstr_k * m)
{
    if(m->m > 0)
    {
        s->mem = realloc(s->mem, (s->m + m->m) * sizeof(char));
        memmove(s->mem + s->m, m->mem, m->m * sizeof(char));
        for(int i = 0; i < m->l; i++) m->offset_ar[i] += s->m;
        s->offset_ar = realloc(s->offset_ar, (s->l + m->l) * sizeof(size_t));
        memmove(s->offset_ar + s->l, m->offset_ar, m->l * sizeof(size_t));
        s->l = s->l + m->l;
        s->m = s->m + m->m;
        arrstrFree(m);
        free(m);
    }
}

static inline size_t Min(size_t a, size_t b)
{
    if(a <= b) return a; else return b;
}

static inline size_t Max(size_t a, size_t b)
{
    if(a >= b) return a; else return b;
}

static inline uint64_t MinUint64(uint64_t a, uint64_t b)
{
    if(a <= b) return a; else return b;
}

static inline uint64_t MaxUint64(uint64_t a, uint64_t b)
{
    if(a >= b) return a; else return b;
}

static inline void kpFree(pkarwl_k * s)
{
    if(s)
    {
        if(s->ar)
        {
            free(s->ar);
            s->ar = NULL;
        }
        if(s->qlar)
        {
            free(s->qlar);
            s->qlar = NULL;
        }
        s->l = 0;
    }
}

static inline int addEleToArwlXC(arwl_k * s, size_t m)
{
    int i;
    for(i = 0; i < s->l; i++)
    {
        if(m == s->ar[i]) return -1;
    }
    if(i == s->l)
    {
        s->ar = realloc(s->ar, (s->l + 1) * sizeof(size_t));
        s->ar[s->l] = m;
        s->l = s->l + 1;
        return 1;
    }
    return 1;
}

static inline void addBarToBackRep(barc_k * s, uint64_t m)
{
    s->bar = realloc(s->bar, (s->l + 1) * sizeof(uint64_t));
    s->bar[s->l] = m;
    s->l = s->l + 1;
}

static inline void kbFree(barc_k * s)
{  
    if(s)
    {
        if(s->bar)
        {
            free(s->bar);
            s->bar = NULL;
        }
        s->l = 0;
    }
}

static inline void addEleToArwlRep(arwl_k * s, size_t m)
{
    s->ar = realloc(s->ar, (s->l + 1) * sizeof(size_t));
    s->ar[s->l] = m;
    s->l = s->l + 1;
}

static inline int addEleToArwlXCSorted(arwl_k * s, size_t m)
{
    int p;

    for(p = 0; p < s->l; p++)
    {
        if(m > s->ar[p])
        {
            continue;
        }
        else if(m < s->ar[p])
        {
            break;
        }
        else
        {
            return -1;
        }
    }

    size_t * tmp = (size_t *) calloc(s->l + 1, sizeof(size_t));
    memmove(tmp, s->ar, (p) * sizeof(size_t));
    tmp[p] = m;
    memmove(tmp + p + 1, s->ar + p, (s->l - p) * sizeof(size_t));
    if(s->ar) free(s->ar);
    s->ar = tmp;
    s->l = s->l + 1;
    return p;
}

static inline size_t addEleToArwlRepSorted(arwl_k * s, size_t m)
{
    size_t * tmp = (size_t *) calloc(s->l + 1, sizeof(size_t));
    size_t i = 0;
    if(s->ar)
    {
        for(i = 0; i < s->l; i++)
        {
            if(m <= s->ar[i]) break;
        }
        memmove(tmp, s->ar, (i) * sizeof(size_t));
        tmp[i] = m;
        memmove(tmp + i + 1, s->ar + i, (s->l - i) * sizeof(size_t));
    }
    else
    {
        s->ar = tmp;
        tmp[0] = m;
    }
    s->ar = tmp;
    s->l = s->l + 1;
    return i;
}

static inline void insertEleToArwl(arwl_k * s, size_t m, size_t pos)
{
    if(pos >= s->l) 
    {
        addEleToArwlRep(s, m);
    }
    else
    {
        size_t * tmp = (size_t *) calloc(s->l + 1, sizeof(size_t));
        memmove(tmp, s->ar, (pos) * sizeof(size_t));
        tmp[pos] = m;
        memmove(tmp + pos + 1, s->ar + pos, (s->l - pos) * sizeof(size_t));
        if(s->ar) free(s->ar);
        s->ar = tmp;
        s->l = s->l + 1;
    }

}

static inline int addStrToArrstrXC(arrstr_k * s, char * str)
{
    for(int i = 0; i < s->l; i++)
    {
        if(strcmp(s->mem + s->offset_ar[i], str) == 0) return 0;
    }
    size_t orim = s->m;
    s->m = s->m + strlen(str);
    s->mem = realloc(s->mem, s->m * sizeof(char) + 1);
    s->offset_ar = realloc(s->offset_ar, (s->l + 1) * sizeof(size_t));
    s->offset_ar[s->l] = orim;
    strcpy(s->mem + orim, str);
    s->mem[s->m++] = '\0';
    s->l = s->l + 1;
    return 1;
}

static inline int addStrToArrstrRepXC(arrstr_k * s, char * str)
{
    if(str == NULL)
    {
        return 0;
    }
    else
    {
        size_t orim = s->m;
        s->m = s->m + strlen(str);
        s->mem = realloc(s->mem, s->m * sizeof(char) + 1);
        s->offset_ar = realloc(s->offset_ar, (s->l + 1) * sizeof(size_t));
        s->offset_ar[s->l] = orim;
        strcpy(s->mem + orim, str);
        s->mem[s->m++] = '\0';
        s->l = s->l + 1;
        return 0;
    }
}

// if splited return length, if not splited return 1 of original str;
static inline arrstr_k * splitStrByChar(char * str, char sp)
{
    arrstr_k * ret = (arrstr_k *) calloc(1, sizeof(arrstr_k));
    ret->mem = (char *) malloc((strlen(str)+1) * sizeof(char));
    ret->mem = strcpy(ret->mem, str);
    ret->m = strlen(str) + 1;

    size_t * offsetar = NULL, mem_offset = 0;
    char * curep = str, * ptr = str;
    size_t curele = 0;
    for(int i = 0; i < strlen(str) + 1; i++)
    {
        if(* ptr == sp)
        {
            if(ptr-curep == 0)
            {
                curep = ++ptr;
                mem_offset++;
                continue;
            }
            offsetar = (size_t *) realloc(offsetar, (curele + 1) * sizeof(size_t));
            offsetar[curele] = mem_offset;
            mem_offset = mem_offset + ptr - curep;
            ret->mem[mem_offset++] = '\0';
            curele++;
            curep = ++ptr;
        }
        else if(* ptr == '\0')
        {
            if(ptr-curep != 0)
            {
                offsetar = (size_t *) realloc(offsetar, (curele + 1) * sizeof(size_t));
                offsetar[curele] = mem_offset;
                mem_offset = mem_offset + ptr - curep;
                ret->mem[mem_offset++] = '\0';
                curele++;
            }
            break;
        }
        else
        {
            ptr++;
        }
    }
    ret->l = curele;
    ret->offset_ar = offsetar;
    return ret;
}

// if can be splited, return length >= 1; if split length less than string length l = 0;
static inline arrstr_k * splitStrByStr(char * str, char * sp)
{
    arrstr_k * ret = (arrstr_k *) calloc(1, sizeof(arrstr_k));
    char * curep = str, * ptr = str;
    size_t curele = 0;
    size_t spl = strlen(sp);
    size_t strl = strlen(str);
    size_t * offsetar = NULL, mem_offset = 0;
    if(strl < spl)
    {
        ret->l = 0;
    }
    else
    {
        ret->mem = (char *) malloc((strl + 1) * sizeof(char));
        for(int i = 0; i < (strl - spl) + 1; i++)
        {
            if(strlen(ptr) >= strlen(sp))
            {
                if(strncmp(ptr, sp, spl) == 0)
                {
                    if(ptr - curep == 0)
                    {
                        ptr = ptr + spl;
                        curep = ptr;
                        continue;
                    }
                    else
                    {
                        offsetar = (size_t *) realloc(offsetar, (curele + 1) * sizeof(size_t));
                        offsetar[curele] = mem_offset;
                        strncpy(ret->mem + mem_offset, curep, ptr-curep);
                        mem_offset = mem_offset + ptr - curep;
                        ret->mem[mem_offset++] = '\0';
                        curele++;
                        curep = ptr = ptr + spl;
                    }
                }          
                else
                {
                    ptr++;
                }
            }
            else
            {
                offsetar = (size_t *) realloc(offsetar, (curele + 1) * sizeof(size_t));
                offsetar[curele] = mem_offset;
                strncpy(ret->mem + mem_offset, curep, strlen(curep));
                mem_offset = mem_offset + strlen(curep);
                ret->mem[mem_offset++] = '\0';
                curele++;
                break;
            }
        }
    }
    ret->offset_ar = offsetar;
    ret->m = mem_offset;
    ret->mem = realloc(ret->mem, ret->m);
    ret->l = curele;
    return ret;
}

static inline char * strConcate(char * s, char * m, char * c)
{
    char * tmp = (char *) calloc(strlen(s) + strlen(m) + strlen(c) + 1, sizeof(char));
    memmove(tmp, s, strlen(s) * sizeof(char));
    memmove(tmp + strlen(s), c, strlen(c) * sizeof(char));
    memmove(tmp + strlen(s) + strlen(c), m, strlen(m) * sizeof(char));
    tmp[strlen(s) + strlen(m) + strlen(c)] = '\0';
    return tmp;
}

static inline char * strConcateByChar(char * s, char * m, char c)
{
    char * tmp = (char *)calloc(strlen(s) + strlen(m) + 2, sizeof(char));
    memmove(tmp, s, strlen(s) * sizeof(char));
    tmp[strlen(s)] = c;
    memmove(tmp + strlen(s) + 1, m, strlen(m) * sizeof(char));
    tmp[strlen(s) + strlen(m) + 1] = '\0';
    return tmp;
}

static inline char * strCopy(char * ori)
{
    char * tmp;
    if(ori)
    {
        tmp = (char *) calloc(1, (strlen(ori) + 1) * sizeof(char));
        strcpy(tmp, ori);
        tmp[strlen(ori)] = '\0';
    }
    else
    {
        tmp = NULL;
    }
    return tmp;
}

static inline arrstr_k * arrstrCopy(arrstr_k * ori)
{
    arrstr_k * target;
    if(ori)
    {
        target = (arrstr_k *) calloc(1, sizeof(arrstr_k));
        target->l = ori->l;
        target->m = ori->m;
        target->mem = (char *) malloc(ori->m);
        target->offset_ar = (size_t *) malloc(ori->l * sizeof(size_t));
        memmove(target->offset_ar, ori->offset_ar, ori->l * sizeof(size_t));
        memmove(target->mem, ori->mem, ori->m);
    } 
    else
    {   
        target = NULL;
    }
    return target;
}

// move from ptr to target
static inline char * movePtrAlongStr(char * ptr, char * target, int step)
{
    if(ptr && target)
    {
        char * tend = target + strlen(target);
        if((tend - ptr) <= (int) step)
        {
            return NULL;
        }
        else
        {
            ptr = ptr + (int) step;
        }
        return ptr;
    }
    else
    {
        return NULL;
    }

}

static inline char * convertSizetToStr(size_t s)
{
    char * str = (char *) calloc(256, sizeof(char));
    snprintf(str, 256, "%zu", s);
    return str;
}
// calculate C m n
static inline size_t cmnCount(size_t m, size_t p)
{
    size_t t = p, r = m;
    for(int i = 1; i < p; i++) 
    {
        r = r * (m - i);
        t = t * i;
    }
    size_t ret = r/t;
    return ret;
}

// sort arwl: sort > 0 max --> min else min --> max
static inline size_t * selectSortOrder(arwl_k * s, int sort)
{
    size_t min, temp, oritemp;
    arwl_k * result = (arwl_k *) calloc(1, sizeof(arwl_k));
    result->l = s->l;
    result->ar = (size_t *) calloc(result->l, sizeof(size_t));
    memmove(result->ar, s->ar, s->l * sizeof(size_t));

    size_t * OriOrder = (size_t *) calloc(s->l, sizeof(size_t));
    for(int i = 0; i < s->l; i++) OriOrder[i] = i;

    for(int i = 0; i < result->l; i++)
    {
        min = i;
        for(int j = i + 1; j < result->l; j++)
        {
            if(sort > 0)
            {
                if(result->ar[j] > result->ar[min]) min = j;
            } 
            else
            {
                if(result->ar[j] < result->ar[min]) min = j;
            }
        }
        if(min != i)
        {
            temp = result->ar[min];
            oritemp = OriOrder[min];

            result->ar[min] = result->ar[i];
            OriOrder[min] = OriOrder[i];

            result->ar[i] = temp;
            OriOrder[i] = oritemp;
        }
    }
    kaFree(result);
    free(result);
    return OriOrder;
}

// return a new awrl of sorted .
static inline arwl_k * selectSortArray(arwl_k * s, int sort)
{
    size_t min, temp, oritemp;
    arwl_k * result = (arwl_k *) calloc(1, sizeof(arwl_k));
    result->l = s->l;
    result->ar = (size_t *) calloc(result->l, sizeof(size_t));
    memmove(result->ar, s->ar, s->l * sizeof(size_t));

    size_t * OriOrder = (size_t *) calloc(s->l, sizeof(size_t));
    for(int i = 0; i < s->l; i++) OriOrder[i] = i;

    for(int i = 0; i < result->l; i++)
    {
        min = i;
        for(int j = i + 1; j < result->l; j++)
        {
            if(sort > 0)
            {
                if(result->ar[j] > result->ar[min]) min = j;
            } 
            else
            {
                if(result->ar[j] < result->ar[min]) min = j;
            }
        }
        if(min != i)
        {
            temp = result->ar[min];
            oritemp = OriOrder[min];

            result->ar[min] = result->ar[i];
            OriOrder[min] = OriOrder[i];

            result->ar[i] = temp;
            OriOrder[i] = oritemp;
        }
    }
    free(OriOrder);
    return result;
}
#endif