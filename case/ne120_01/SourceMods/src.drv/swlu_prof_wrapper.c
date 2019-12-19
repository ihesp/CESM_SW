#include <stdlib.h>
#include <mpi.h>

static inline int env_swlu_prof_enabled(){
  char *swlu_enable = getenv("SWLU_PROF");
  if (swlu_enable != NULL && !strcmp(swlu_enable, "TRUE"))
    return 1;
  return 0;
}

static int swlu_prof_enabled = 0;
extern int swlu_prof_index;
extern int swlu_prof_set_comm;
extern MPI_Comm swlu_prof_comm;
extern int swlu_prof_print_every;

void cesm_swlu_prof_set_print_every_(const int *x)
{
    swlu_prof_print_every = *x;
}

void cesm_swlu_prof_init_withcomm_(const int *index, const MPI_Comm *comm)
{
  if (swlu_prof_enabled) {
    swlu_prof_index = *index;
    swlu_prof_set_comm = 1;
    swlu_prof_comm = *comm;
    swlu_prof_init_();
    swlu_prof_index = 0;
    swlu_prof_set_comm = 0;
  }
}

void cesm_swlu_prof_init_(const int *index)
{
  if (env_swlu_prof_enabled()) {
    swlu_prof_enabled = 1;
    swlu_prof_index = *index;
    swlu_prof_init_();
    swlu_prof_index = 0;
  }
}

void cesm_swlu_debug_init_()
{
  swlu_debug_init_();
}

void cesm_swlu_prof_print_(const int *index)
{
  if (swlu_prof_enabled) {
    swlu_prof_index = *index;
    swlu_prof_print_();
    swlu_prof_index = 0;
  }
}

void cesm_swlu_prof_start_(const int *index)
{
  if (swlu_prof_enabled) {
    swlu_prof_index = *index;
    swlu_prof_start_();
    swlu_prof_index = 0;
  }
}

void cesm_swlu_prof_stop_(const int *index)
{
  if (swlu_prof_enabled) {
    swlu_prof_index = *index;
    swlu_prof_stop_();
    swlu_prof_index = 0;
  }
}
