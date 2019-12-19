#define FCAT_INNER(x, y) x ## _ ## y ## _
#define FCAT(x, y) FCAT_INNER(x, y)
#define K(x) FCAT(LWPF_UNIT, x)(){lwpf_start(x);}
#define U(x) void lwpf_start_ ## x
LWPF_KERNELS
#undef U
#undef K

#define K(x) FCAT(LWPF_UNIT, x)(){lwpf_stop(x);}
#define U(x) void lwpf_stop_ ## x
LWPF_KERNELS
#undef U
#undef K


#define U(x) lwpf_enter_ ## x ## _ () {lwpf_sync_counters_m2c(lwpf_global_counter_ ## x[_MYID], lwpf_kernel_count_ ## x);}
LWPF_UNIT
#undef U
#define U(x) lwpf_exit_ ## x ## _() {lwpf_sync_counters_c2m(lwpf_global_counter_ ## x[_MYID], lwpf_kernel_count_ ## x);}
LWPF_UNIT
#undef U
