// Author: Xiaohui Duan
// A set of Macros for dynamically sized LDM allocation.
// A safe implementation comes with error checking feature.
// The FAST implementation disabled error checking.
// #define LDM_FAST to use FAST implementation, while by default,
// the safe version is used.
// PLEASE DO NOT USE LDM_FAST at begining.
// LDM_FAST can be tried when the safe implementation can produce
// excepted results.

#ifndef LDM_ALLOC_H
#define LDM_ALLOC_H

// The default starting address for allocating LDM.
// 512B should be a safe value for programs without user defined
// __thread_local variables.
// If you use __thread_local variables, use ldmreport to check
// LDM size used by local variables, and define this value before
// including this header.
#ifndef LDM_ALLOC_BASE_ADDR
#define LDM_ALLOC_BASE_ADDR 2048
#endif

// The default alignment for allocating LDM.
// You can define this value before including this header.
// If you know what you are doing.
#ifndef LDM_ALLOC_ALIGNMENT
#define LDM_ALLOC_ALIGNMENT 32
#endif

#ifndef LDM_FAST
#include "cpe_print.h"
//#include <stdio.h>
//#include <stdlib.h>
#include <assert.h>
// Check is current allocated memory excceeded the sp, which should
// be warned.
#define check_sp(ldm_sp, name) {                                \
    void *sp;                                                   \
    asm volatile("ldi %0, 0($sp)\n\t" : "=r"(sp));              \
    if (ldm_sp > sp){                                           \
      cpe_printf("Too much LDM allocated for " #name "! "       \
          "LDM allocated to %x while sp is %x\n", ldm_sp, sp);  \
      assert(ldm_sp <= sp);                                     \
    }                                                           \
  }                                                             \

// Check is the dellocated variable is the last allocated pointer.
// If not, it should be warned.
#define check_last(ldm_sp, name) {                                      \
    if (name ## _end_shadow != ldm_sp){                                 \
      cpe_printf(#name " is not the last allocated pointer!\n");        \
      assert(name ## _end_shadow == ldm_sp);                            \
    }                                                                   \
  }

#define save_tail(ldm_sp, name) void *name = ldm_sp;
#else

// In fast mode, these two check is disabled.
#define check_sp(ldm_sp, name) {}
#define check_last(ldm_sp, name) {}
#define save_tail(ldm_sp, name){}
#endif

// Small tool to align pointer.
#define LDM_ALLOC_ALIGNMENT_MASK (LDM_ALLOC_ALIGNMENT - 1)
#define align_pointer(ptr, align_mask) \
  ((void*)((long)ptr + align_mask & ~align_mask))

// Initialize a pointer for LDM allocation
#define ldm_alloc_init()                                \
  void *ldm_sp_shadow = (void*)LDM_ALLOC_BASE_ADDR;     \

// Actually, we should define macros for passing allocation state between
// functions.
// Thus, defining void some_func(int arg1, int arg2, LDM_ALLOC_INHERIT)
// and calling some_func(arg1, arg2, arg3, LDM_ALLOC_PASS), will pass the
// LDM allocation pointer to callee.
// Of course, callee should not deallocate data from caller.
#define LDM_ALLOC_PASS ldm_sp_shadow
#define LDM_ALLOC_INHERIT void *ldm_sp_shadow

// Allocate size bytes for ptr
#define ldm_alloc_pointer(ptr, size)                                    \
  ptr = ldm_sp_shadow;                                                  \
  {                                                                     \
    ldm_sp_shadow += size;                                              \
    ldm_sp_shadow = align_pointer(ldm_sp_shadow, LDM_ALLOC_ALIGNMENT_MASK); \
    check_sp(ldm_sp_shadow, name);                                      \
  }                                                                     \
  save_tail(ldm_sp_shadow, ptr ## _end_shadow);

// Recycle an allocated pointer with safety checking.
// If it is not the last allocated pointer, error will be thrown.
#define ldm_dealloc(name) {                     \
    check_last(ldm_sp_shadow, name);            \
    ldm_sp_shadow = name;                       \
  }

// Recycle LDM after specific addr, no safty checking provided.
#define ldm_dealloc_after(addr) {               \
    ldm_sp_shadow = addr;                       \
  }
#endif

static inline void *ldm_alloc_pop(long size, void **ldm_sp) {
  long ret = *ldm_sp;
  *ldm_sp += size;
  *ldm_sp = align_pointer((*ldm_sp), LDM_ALLOC_ALIGNMENT_MASK);
#ifndef LDM_FAST
  check_sp((*ldm_sp), "pop2");
#endif
  return ret;
}

#ifdef POP2_LDM_MALLOC_COMPAT
#define ldm_malloc(size) ldm_alloc_pop(size, &ldm_sp_shadow)
#define ldm_free(ptr, size) {}
#endif
