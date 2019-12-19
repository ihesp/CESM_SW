#include <slave.h>
#include <stdarg.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include "math_data.h"

#define pe_get(src, dst, size) { \
    volatile unsigned long get_reply; \
    get_reply = 0; \
    athread_get(PE_MODE, (src), (dst), (size), &get_reply, 0, 0, 0); \
    while(get_reply != 1); \
    asm volatile("memb"); \
}
#define pe_put(src, dst, size) { \
    volatile unsigned long put_reply; \
    put_reply = 0; \
    athread_put(PE_MODE, (src), (dst), (size), &put_reply, 0, 0); \
    while(put_reply != 1); \
    asm volatile("memb"); \
}

#define bcast_get(src, dst, size) { \
    volatile unsigned long get_reply; \
    athread_syn(ARRAY_SCOPE, 0xffff); \
    if (_MYID == 0) { \
        get_reply = 0; \
        athread_get(BCAST_MODE, (src), (dst), (size), &get_reply, 0xff, 0, 0);  \
        while(get_reply != 1); \
        asm volatile("memb"); \
    } \
    athread_syn(ARRAY_SCOPE, 0xffff); \
}

// __thread_local int _MPI_RANK = 0;

void set_exp_data_(double *exp_local_data)
{
    exp_data_local_ptr = exp_local_data;
    bcast_get(exp_data, exp_local_data, exp_data_bytes);
}

void reset_exp_data_()
{
    exp_data_local_ptr = exp_data;
}

void athread_get_id_(int *ret)
{
    *ret = _MYID;
}

void athread_get_vid_(int *id)
{
    athread_get_id_(id);
    int r = *id / 8, c = *id % 8;
    *id = r % 2 * 32 + c * 4 + r / 2;
}

// void set_mpi_rank_(int *mpi_rank)
// {
//     _MPI_RANK = *mpi_rank;
// }

int trim_strcmp_(char *a, char *b, int na, int nb)
{
    int i, n;
    for (; a[na - 1] == ' '; --na);
    for (; b[nb - 1] == ' '; --nb);
    n = na < nb ? na : nb;
    for (i = 0; i < n; ++i) {
        if (a[i] != b[i]) {
            return a[i] < b[i] ? -1 : 1;
        }
    }
    return na == nb ? 0 : na < nb ? -1 : 1;
}

// void addr_(char *s, void *addr, long n)
// {
//     char ss[n + 1];
//     strncpy(ss, s, n);
//     ss[n] = 0;
//     if (_MPI_RANK == 0 && _MYID == 0)
//         fprintf(stderr, "addr of %s is 0x%lx.\n", ss, (long) addr);
// }

// void get_sp_(char *s, long n)
// {
//     char ss[n + 1];
//     strncpy(ss, s, n);
//     ss[n] = 0;
//     volatile long sp_addr;
//     asm volatile("mov $sp, %0":"=r"(sp_addr));
//     if (_MPI_RANK == 0 && _MYID == 0)
//         fprintf(stderr, "sp of %s is 0x%lx.\n", ss, (long) sp_addr);
// }

void setftz_()
{
    long fpcr;
    asm volatile("rfpcr %0": "=r"(fpcr));
    fpcr |= 1L << 60 | 1L << 2;
    asm volatile("wfpcr %0": : "r"(fpcr));
}

void unsetftz_()
{
    long fpcr;
    asm volatile("rfpcr %0": "=r"(fpcr));
    fpcr |= 1L << 2;
    fpcr &= ~(1L << 60);
    asm volatile("wfpcr %0": : "r"(fpcr));
}

inline void dmemcpy(void *dst, void *src, uint32_t size)
{
    if (dst < 0x10000 && src < 0x10000) {
        int i;
        for (i = 0; i < size / 4; ++i)
            ((int*) dst)[i] = ((int*) src)[i];
    }
    else if (dst > 0x10000 && src < 0x10000) {
        pe_put(src, dst, size);
    }
    else if (dst < 0x10000 && src > 0x10000) {
        pe_get(src, dst, size);
    }
    else {
        int buf[size / 4];
        pe_get(src, buf, size);
        pe_put(buf, dst, size);
    }
    return 0;
}

inline int dmemset(void *dst, int v, uint32_t size)
{
    int i;
    if (dst < 0x10000) {
        for (i = 0; i < size / 4; ++i)
            ((int*) dst)[i] = v;
    }
    else {
        int buf[size / 4];
        for (i = 0; i < size / 4; ++i)
            buf[i] = v;
        pe_put(buf, dst, size);
    }
    return 0;
}

void dmemcpy_(void *dst, void *src, uint32_t *size)
{
    dmemcpy(dst, src, *size);
}

void dmemset_(void *dst, int *v, uint32_t *size)
{
    dmemset(dst, *v, *size);
}
