void athread_spawn_rrtmg_lw_(void *param)
{
    volatile long addr;
    asm volatile(".globl slave_RRTMG_LW_1PARAM.in.RRTMG_LW_RAD_CPE");
    asm volatile("ldl %0,slave_RRTMG_LW_1PARAM.in.RRTMG_LW_RAD_CPE($gp) !literal":"=r"(addr));
    __real_athread_spawn(addr, param);
}

void athread_spawn_modal_aero_sw_kern_(void *param)
{
    volatile long addr;
    asm volatile(".globl slave_MODAL_AERO_SW_KERN.in.MODAL_AER_OPT_CPE");
    asm volatile("ldl %0,slave_MODAL_AERO_SW_KERN.in.MODAL_AER_OPT_CPE($gp) !literal":"=r"(addr));
    __real_athread_spawn(addr, param);
}

void athread_spawn_modal_aero_lw_kern_(void *param)
{
    volatile long addr;
    asm volatile(".globl slave_MODAL_AERO_LW_KERN.in.MODAL_AER_OPT_CPE");
    asm volatile("ldl %0,slave_MODAL_AERO_LW_KERN.in.MODAL_AER_OPT_CPE($gp) !literal":"=r"(addr));
    __real_athread_spawn(addr, param);
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#define MAX_MALLOC_QUEUE_NUM (9)
#define MAX_MALLOC_SPEC_ADDR_NUM (11)
#define MALLOC_BASE_NUM (70000)
#define MAX_MALLOC_ENTRY_NUM (MALLOC_BASE_NUM * 487)
// #define MAX_MALLOC_ENTRY_NUM (14438400)
#define assert(x) { \
    if (!(x)) { \
        fprintf(stderr, "Assert Failed in %s:%d, \"%s\"\n", __FILE__, __LINE__, #x); \
        MPI_Abort(MPI_COMM_WORLD, 11); \
    } \
}
// #define assert(x) {}

// static const int queue_len[MAX_MALLOC_QUEUE_NUM] = {5000, 2000, 200, 50};
// static const int queue_bits[MAX_MALLOC_QUEUE_NUM] = {4096, 8192, 131072, 1048576};
static const int queue_len[MAX_MALLOC_QUEUE_NUM] = {MALLOC_BASE_NUM, MALLOC_BASE_NUM*3, MALLOC_BASE_NUM, MALLOC_BASE_NUM, MALLOC_BASE_NUM, MALLOC_BASE_NUM, MALLOC_BASE_NUM, 1, 1};
static const int queue_bits[MAX_MALLOC_QUEUE_NUM] = {40, 48, 128, 256, 296, 1024, 1168, 56, 784};

static long *malloc_buffer[MAX_MALLOC_QUEUE_NUM];
// static int malloc_h_buffer[MAX_MALLOC_ENTRY_NUM / (4096 / sizeof(long))];
static int malloc_h_buffer[MAX_MALLOC_ENTRY_NUM / (40 / sizeof(long))];
static long malloc_e_buffer[MAX_MALLOC_ENTRY_NUM];
static int qn[MAX_MALLOC_QUEUE_NUM], l[MAX_MALLOC_QUEUE_NUM], r[MAX_MALLOC_QUEUE_NUM], *h[MAX_MALLOC_QUEUE_NUM];
static int queue_malloc_initialized = 0;
static int spec_maddr_initialized = 0;
static long spec_maddr[MAX_MALLOC_SPEC_ADDR_NUM];

static int queue_malloc_on = 0;
static int queue_malloc_mask = 0;
static int malloc_cnt_on = 0, malloc_cnt = 0;

void print_malloc_cnt_()
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank < 5)
        fprintf(stderr, "rank%d malloc cnt = %d\n", rank, malloc_cnt);
}

void set_malloc_cnt_on_()
{
    malloc_cnt_on = 1;
}

void set_malloc_cnt_off_()
{
    malloc_cnt_on = 0;
}

void reset_malloc_cnt_()
{
    malloc_cnt = 0;
}

void print_queue_malloc(char *caller)
{
    int k;
    char pbuf[128], *p;
    p = pbuf;
    p += sprintf(p, "%s :", caller);
    for (k = 0; k < MAX_MALLOC_QUEUE_NUM; ++k)
        p += sprintf(p, "%d ", qn[k]);
    *p = 0;
    fprintf(stderr, "%s\n", pbuf);
}

void print_queue_malloc_()
{
    print_queue_malloc(__func__);
}

void set_queue_malloc_on_()
{
    // print_queue_malloc(__func__);
    queue_malloc_on = 1;
}

void set_queue_malloc_off_()
{
    // print_queue_malloc(__func__);
    queue_malloc_on = 0;
}

void set_queue_malloc_mask()
{
    queue_malloc_mask = 1;
}
void reset_queue_malloc_mask()
{
    queue_malloc_mask = 0;
}
void set_queue_malloc_mask_()
{
    set_queue_malloc_mask();
}
void reset_queue_malloc_mask_()
{
    reset_queue_malloc_mask();
}

static void queue_malloc_init()
{
    if (getenv("QUEUE_MALLOC_ON") && !strcmp(getenv("QUEUE_MALLOC_ON"), "TRUE"))
        set_queue_malloc_on_();
    int k, i;
    for (k = 0; k < MAX_MALLOC_QUEUE_NUM; ++k) {
        l[k] = 0, r[k] = 0, qn[k] = 0;
        h[k] = (k == 0 ? malloc_h_buffer : h[k-1] + queue_len[k-1]);
        malloc_buffer[k] = (k == 0 ? malloc_e_buffer : malloc_buffer[k-1] + queue_len[k-1] * (queue_bits[k-1] / sizeof(long)));
        assert(h[k] + queue_len[k] <= malloc_h_buffer + MAX_MALLOC_ENTRY_NUM / (queue_bits[0] / sizeof(long)));
        assert(malloc_buffer[k] + queue_len[k] * (queue_bits[k] / sizeof(long)) <= malloc_e_buffer + MAX_MALLOC_ENTRY_NUM);
        memset(malloc_buffer[k], 0, queue_len[k] * queue_bits[k]);
        for (i = 0; i < queue_len[k]; ++i) {
            h[k][i] = i;
        }
    }
    volatile long addr;
    asm volatile(".globl MPIDI_CH3_VC_Init");
    asm volatile("ldl %0,MPIDI_CH3_VC_Init($gp) !literal":"=r"(addr));
    spec_maddr[0] = addr + 84;
    asm volatile(".globl rdma_iba_hca_init");
    asm volatile("ldl %0,rdma_iba_hca_init($gp) !literal":"=r"(addr));
    spec_maddr[1] = addr + 1216;
    spec_maddr[2] = addr + 1272;
    asm volatile(".globl MRAILI_Init_vc");
    asm volatile("ldl %0,MRAILI_Init_vc($gp) !literal":"=r"(addr));
    spec_maddr[3] = addr + 264;
    spec_maddr[4] = addr + 288;
    spec_maddr[5] = addr + 404;
    asm volatile(".globl mlx4_alloc_qp_buf");
    asm volatile("ldl %0,mlx4_alloc_qp_buf($gp) !literal":"=r"(addr));
    spec_maddr[6] = addr + 56;
    spec_maddr[7] = addr + 344;
    // asm volatile(".globl mlx4_create_qp_common");
    // asm volatile("ldl %0,mlx4_create_qp_common($gp) !literal":"=r"(addr));
    // spec_maddr[5] = swlu_find_sym_from_str("mlx4_create_qp_common")->addr + 176;
    asm volatile(".globl mlx4_destroy_ah");
    asm volatile("ldl %0,mlx4_destroy_ah($gp) !literal":"=r"(addr));
    spec_maddr[8] = addr + 224;
    asm volatile(".globl cm_qp_create");
    asm volatile("ldl %0,cm_qp_create($gp) !literal":"=r"(addr));
    spec_maddr[9] = addr + 108;
    spec_maddr[10] = addr + 164;
}

static void *queue_malloc(int k)
{
    ++qn[k];
    if (qn[k] > queue_len[k]) { fprintf(stderr, "k%d\n", k); abort();}
    assert(qn[k] <= queue_len[k]);
    void *ptr = (void*) (malloc_buffer[k] + h[k][l[k]] * (queue_bits[k] / sizeof(long)));
    l[k]  = (l[k] + 1) % queue_len[k];
    return ptr;
}

static void queue_free(int k, void *ptr)
{
    --qn[k];
    assert(qn[k] >= 0);
    h[k][r[k]] = ((long*) ptr - malloc_buffer[k]) / (queue_bits[k] / sizeof(long));
    assert(0 <= h[k][r[k]] && h[k][r[k]] < queue_len[k]);
    r[k] = (r[k] + 1) % queue_len[k];
}

// void* __attribute__((noinline)) __wrap_malloc(size_t size)
// {
//     if (malloc_cnt_on)
//         ++malloc_cnt;
//     return __real_malloc(size);
// }
// 
// void __attribute__((noinline)) __wrap_free(void *ptr)
// {
//     if (malloc_cnt_on)
//         --malloc_cnt;
//     return __real_free(ptr);
// }

// /*
int __attribute__((noinline)) __wrap_posix_memalign(void **memptr, size_t alignment, size_t size)
{
    int ret = __real_posix_memalign(memptr, alignment, size);
    assert(!ret);
    return ret;
}
void* __attribute__((noinline)) __wrap_valloc(size_t size)
{
    int alignment = 8192; // page_size
    void *ptr = __real_malloc(size + alignment - 1);
    // void *ptr = __real_valloc(size);
    assert(ptr);
    // ptr = (void*) (((long) ptr + alignment - 1) / alignment * alignment);
    return ptr;
}
// */

// /*
#define _F90_ALLOCATE_OFFSET 272
#define _ASSIGN_ALLOCATABLE_OFFSET 512
void* __attribute__((noinline)) __wrap_malloc(size_t size)
{
    int k, flag;
    #pragma frequency_hint <never>
    if (!queue_malloc_on || queue_malloc_mask) {
        void *ptr = __real_malloc(size);
        if(!ptr) fprintf(stderr, "Warnning : malloc %d failed\n", size);
        assert(ptr);
        return ptr;
    }
    #pragma frequency_hint <init>
    if (!queue_malloc_initialized) {
        queue_malloc_initialized = 1;
        queue_malloc_init();
    }

    long caller_addr;
    asm volatile("ldl %0, 0($sp)":"=r"(caller_addr));
    flag = 0;
    for (k = 0; k < MAX_MALLOC_SPEC_ADDR_NUM; ++k) {
        if (spec_maddr[k] == caller_addr) {
            flag = 1;
            break;
        }
    }
    if (flag) 
    // void *_F90_ALLOCATE(size_t);
    // void *_ASSIGN_ALLOCATABLE(size_t);
    // void *_f90_alloc_ptr;
    // void *_ass_alloc_ptr;
    // _f90_alloc_ptr = (void*)_F90_ALLOCATE;
    // _ass_alloc_ptr = (void*)_ASSIGN_ALLOCATABLE;
    // if(caller_addr == (long)_f90_alloc_ptr +_F90_ALLOCATE_OFFSET || caller_addr == (long)_ass_alloc_ptr + _ASSIGN_ALLOCATABLE_OFFSET)
        for (k = 0; k < MAX_MALLOC_QUEUE_NUM; ++k)
            // if (size <= queue_bits[k] && qn[k] < queue_len[k])
            if (size == queue_bits[k])
                return queue_malloc(k);
    void *ptr = __real_malloc(size);
    assert(ptr);
    return ptr;
}

// /*
void* __attribute__((noinline)) __wrap_calloc(size_t n, size_t size)
{
    void *ptr = __real_calloc(n, size);
    assert(ptr);
    return ptr;
}

void* __attribute__((noinline)) __wrap_realloc(void *ptr, size_t size)
{
    int k;
    #pragma frequency_hint <init>
    if (!queue_malloc_initialized) {
        queue_malloc_initialized = 1;
        queue_malloc_init();
    }
    for (k = 0; k < MAX_MALLOC_QUEUE_NUM; ++k)
        if (malloc_buffer[k] <= ptr && ptr < malloc_buffer[k] + queue_len[k] * (queue_bits[k] / sizeof(long))) {
            if (size <= queue_bits[k])
                return ptr;
            queue_free(k, ptr);
            void *newptr = __real_malloc(size);
            memcpy(newptr, ptr, queue_bits[k]);
            return newptr;
        }
    assert(!(malloc_buffer[0] <= ptr && ptr < malloc_buffer[MAX_MALLOC_QUEUE_NUM-1] + queue_len[MAX_MALLOC_QUEUE_NUM-1] * (queue_bits[MAX_MALLOC_QUEUE_NUM-1] / sizeof(long))));
    ptr = __real_realloc(ptr, size);
    assert(ptr);
    return ptr;
}

void __attribute__((noinline)) __wrap_free(void *ptr)
{
    int k;
    #pragma frequency_hint <init>
    if (!queue_malloc_initialized) {
        queue_malloc_initialized = 1;
        queue_malloc_init();
    }
    for (k = 0; k < MAX_MALLOC_QUEUE_NUM; ++k) {
        if (malloc_buffer[k] <= ptr && ptr < malloc_buffer[k] + queue_len[k] * (queue_bits[k] / sizeof(long))) {
            queue_free(k, ptr);
            return;
        }
    }
    assert(!(malloc_buffer[0] <= ptr && ptr < malloc_buffer[MAX_MALLOC_QUEUE_NUM-1] + queue_len[MAX_MALLOC_QUEUE_NUM-1] * (queue_bits[MAX_MALLOC_QUEUE_NUM-1] / sizeof(long))));
    __real_free(ptr);
}
// */
