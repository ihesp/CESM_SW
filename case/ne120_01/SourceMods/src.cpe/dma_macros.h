// Author: Xiaohui Duan
// A set of Macros with DMA and athread_get/put implementation.
// #define DMA_FAST to use DMA implementation, by default, athread
// implementation is used.
// PLEASE DO NOT USE DMA_FAST at begining.
// DMA_FAST can be tried when athread_get/put implementation can
// produce excepted results.

#ifndef DMA_MACROS_H
#define DMA_MACROS_H
#ifdef DMA_FAST

#include <slave.h>
#include <dma.h>
// A magic macro that pretends that reply depends on dma instruction.
// Thus compiler will do less evil.
#define dma_rpl(desc, mem, ldm, reply)                  \
  asm("dma %0, %1, %2\n\t"                              \
      : : "r"(desc), "r"(mem), "r"(ldm), "r"(&reply)    \
      : "memory")

// Pretending to be initializing, actually, it is defining descriptors,
// reply and count variable.
#define dma_init()                                      \
  volatile int reply_shadow = 0;                        \
  int count_shadow = 0;                                 \
  dma_desc pe_get_desc = 0, pe_put_desc = 0;            \
  dma_desc row_get_desc = 0, row_put_desc = 0;          \
  dma_desc rank_get_desc = 0, rank_put_desc = 0;        \
  dma_desc bcast_get_desc = 0, brow_get_desc = 0;       \
  dma_set_op(&pe_get_desc, DMA_GET);                    \
  dma_set_op(&pe_put_desc, DMA_PUT);                    \
  dma_set_op(&row_get_desc, DMA_GET);                   \
  dma_set_op(&row_put_desc, DMA_PUT);                   \
  dma_set_op(&rank_get_desc, DMA_GET);                  \
  dma_set_op(&rank_put_desc, DMA_PUT);                  \
  dma_set_op(&bcast_get_desc, DMA_GET);                 \
  dma_set_op(&brow_get_desc, DMA_GET);                  \
                                                        \
  dma_set_mode(&pe_get_desc, PE_MODE);                  \
  dma_set_mode(&pe_put_desc, PE_MODE);                  \
  dma_set_mode(&row_get_desc, ROW_MODE);                \
  dma_set_mode(&row_put_desc, ROW_MODE);                \
  dma_set_mode(&rank_get_desc, RANK_MODE);              \
  dma_set_mode(&rank_put_desc, RANK_MODE);              \
  dma_set_mode(&bcast_get_desc, BCAST_MODE);            \
  dma_set_mode(&brow_get_desc, BROW_MODE);              \
                                                        \
  dma_set_reply(&pe_get_desc, &reply_shadow);           \
  dma_set_reply(&pe_put_desc, &reply_shadow);           \
  dma_set_reply(&row_get_desc, &reply_shadow);          \
  dma_set_reply(&row_put_desc, &reply_shadow);          \
  dma_set_reply(&rank_get_desc, &reply_shadow);         \
  dma_set_reply(&rank_put_desc, &reply_shadow);         \
  dma_set_reply(&bcast_get_desc, &reply_shadow);        \
  dma_set_reply(&brow_get_desc, &reply_shadow);         \
                                                        \
  dma_set_mask(&bcast_get_desc, 0xff);                  \
  dma_set_mask(&brow_get_desc, 0xff);                   \


// Set stride and bsize for strided get and rank mode.
#define dma_set_stride(mode, stride, bsize) {           \
    dma_set_stepsize(&(mode ## _desc), (long)(stride)); \
    dma_set_bsize(&(mode ## _desc), (long)(bsize));     \
  }                                                     \

// Set mask for brow and bcast mode
#define bcast_set_mask(mode, mask){                     \
    dma_set_mask(&(mode ## _desc), (long)(mask));       \
  }

// PE_MODE, DMA_GET, stride and bsize can be set by dma_set_stride
#define pe_get(mem, ldm, size) {                        \
    dma_set_size(&pe_get_desc, (long)(size));           \
    dma_rpl(pe_get_desc, (mem), (ldm), reply_shadow);   \
    count_shadow ++;                                    \
  }

// PE_MODE, DMA_PUT, stride and bsize can be set by dma_set_stride
#define pe_put(mem, ldm, size) {                        \
    dma_set_size(&pe_put_desc, (long)(size));           \
    dma_rpl(pe_put_desc, (mem), (ldm), reply_shadow);   \
    count_shadow ++;                                    \
  }

// ROW_MODE, DMA_GET, stride and bsize can be set by dma_set_stride
#define row_get(mem, ldm, size) {                       \
    dma_set_size(&row_get_desc, (long)(size));          \
    dma_rpl(row_get_desc, (mem), (ldm), reply_shadow);  \
    count_shadow ++;                                    \
  }

// ROW_MODE, DMA_PUT, stride and bsize can be set by dma_set_stride
#define row_put(mem, ldm, size) {                       \
    dma_set_size(&row_put_desc, (long)(size));          \
    dma_rpl(row_put_desc, (mem), (ldm), reply_shadow);  \
    count_shadow ++;                                    \
  }

// RANK_MODE, DMA_GET, bsize can be set by dma_set_stride
#define rank_get(mem, ldm, size) {                      \
    dma_set_size(&rank_get_desc, (long)(size));         \
    dma_rpl(rank_get_desc, (mem), (ldm), reply_shadow); \
    count_shadow ++;                                    \
  }

// RANK_MODE, DMA_PUT, stride and bsize can be set by dma_set_stride
#define rank_put(mem, ldm, size) {                      \
    dma_set_size(&rank_put_desc, (long)(size));         \
    dma_rpl(rank_put_desc, (mem), (ldm), reply_shadow); \
    count_shadow ++;                                    \
  }

// BCAST_MODE, DMA_GET, stride and bsize can be set by dma_set_stride,
// broadcasting mask can be set by bcast_set_mask
#define bcast_get(mem, ldm, size) {                             \
    dma_set_size(&bcast_get_desc, (long)(size));                \
    dma_rpl(bcast_get_desc, (mem), (ldm), reply_shadow);        \
    count_shadow ++;                                            \
  }

// BROW_MODE, DMA_GET, stride and bsize can be set by dma_set_stride,
// broadcasting mask can be set by bcast_set_mask
#define brow_get(mem, ldm, size) {                      \
    dma_set_size(&brow_get_desc, (long)(size));         \
    dma_rpl(brow_get_desc, (mem), (ldm), reply_shadow); \
    count_shadow ++;                                    \
  }

#else
#include <slave.h>
// Pretending to be initializing, actually, it is defining reply, count.
// No descriptor required, setting stride, bsize and mask are done by
// setting variables.
// REMINDER: please clear masks, strides and bsizes after using strided or
// rank mode, and, currently it is not so compatible with DMA_FAST since
// all modes share a same set of masks.
#define dma_init()                              \
  volatile int reply_shadow = 0;                \
  int count_shadow = 0;                         \
  long pe_get_stride_shadow = 0;                \
  long pe_put_stride_shadow = 0;                \
  long row_get_stride_shadow = 0;               \
  long row_put_stride_shadow = 0;               \
  long rank_get_stride_shadow = 0;              \
  long rank_put_stride_shadow = 0;              \
  long brow_get_stride_shadow = 0;              \
  long bcast_get_stride_shadow = 0;             \
                                                \
  long pe_get_bsize_shadow = 0;                 \
  long pe_put_bsize_shadow = 0;                 \
  long row_get_bsize_shadow = 0;                \
  long row_put_bsize_shadow = 0;                \
  long rank_get_bsize_shadow = 0;               \
  long rank_put_bsize_shadow = 0;               \
  long brow_get_bsize_shadow = 0;               \
  long bcast_get_bsize_shadow = 0;              \
                                                \
  int bcast_get_mask_shadow = 0xff;             \
  int brow_get_mask_shadow = 0xff;              \

// The comments for the following macros can be found in DMA_FAST mode.
#define dma_set_stride(mode, stride, bsize) {   \
    mode ## _stride_shadow = (long)(stride);    \
    mode ## _bsize_shadow = (long)(bsize);      \
  }

#define bcast_set_mask(mode, mask) {            \
    mode ## _mask_shadow = (int)(mask);          \
  }

#define pe_get(mem, ldm, size) {                                        \
    athread_get(PE_MODE, (mem), (ldm), (size), (void*)&reply_shadow     \
        , 0, pe_get_stride_shadow, pe_get_bsize_shadow);                \
    count_shadow ++;                                                    \
  }

#define row_get(mem, ldm, size) {                                       \
    athread_get(ROW_MODE, (mem), (ldm), (size), (void *)&reply_shadow   \
        , 0, row_get_stride_shadow, row_get_bsize_shadow);              \
    count_shadow ++;                                                    \
  }

#define rank_get(mem, ldm, size) {                                      \
    athread_get(RANK_MODE, (mem), (ldm), (size), (void *)&reply_shadow  \
        , 0, rank_get_stride_shadow, rank_get_bsize_shadow);            \
    count_shadow ++;                                                    \
  }

#define bcast_get(mem, ldm, size) {                                     \
    athread_get(BCAST_MODE, (mem), (ldm), (size), (void*)&reply_shadow  \
        , bcast_get_mask_shadow, bcast_get_stride_shadow                \
        , bcast_get_bsize_shadow);                                      \
    count_shadow ++;                                                    \
  }

#define brow_get(mem, ldm, size) {                                      \
    athread_get(BROW_MODE, (mem), (ldm), (size), (void *)&reply_shadow  \
        , brow_get_mask_shadow, brow_get_stride_shadow                  \
        , brow_get_bsize_shadow);                                       \
    count_shadow ++;                                                    \
  }

#define pe_put(mem, ldm, size) {                                        \
    athread_put(PE_MODE, (ldm), (mem), (size), (void*)&reply_shadow     \
        , pe_put_stride_shadow, pe_put_bsize_shadow);                   \
    count_shadow ++;                                                    \
  }

#define row_put(mem, ldm, size) {                                       \
    athread_put(ROW_MODE, (ldm), (mem), (size), (void *)&reply_shadow   \
        , row_put_stride_shadow, row_put_bsize_shadow);          \
    count_shadow ++;                                                    \
  }


#define rank_put(mem, ldm, size) {                                      \
    athread_put(RANK_MODE, (ldm), (mem), (size), (void *)&reply_shadow  \
        , rank_put_stride_shadow, rank_put_bsize_shadow);               \
    count_shadow ++;                                                    \
  }

#endif

#define bcast_reset_mask(mode) bcast_set_mask(mode, 0xff)
#define dma_reset_stride(mode) dma_set_stride(mode, 0, 0)

// Sync the reply and count. "memb" is always a good choice to avoid
// god's punishment.
#define dma_syn() {                                     \
    while (reply_shadow != count_shadow) {              \
    };                                                  \
    asm volatile("memb\n\t");                           \
  }

#endif
