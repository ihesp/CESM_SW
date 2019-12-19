#ifndef EDGE_SLAVE_H
#define EDGE_SLAVE_H
typedef double FTYPE;
typedef int ITYPE;
typedef struct {
  double *buf;
  double *receive;
  double *addrs[3];
  int *putmap_addr;
  int *getmap_addr;
  int *reverse_addr;
  int *iwesn_addr;
  int np, nlev, nets, nete,max_corner_elem, iam;
} edge_param_t;

typedef struct {
  double *buf;
  double *receive;
  double *addrs[3];
  int *putmap_addr;
  int *getmap_addr;
  int *reverse_addr;
  int *iwesn_addr;
  int np, nlev, nets, nete;
  int max_corner_elem, offset, rsplit, iam;
} edge3p1_param_t;

typedef struct {
  double *buf;
  double *receive;
  double *addrs[3];
  int *putmap_addr;
  int *getmap_addr;
  int *reverse_addr;
  int *iwesn_addr;
  int np, nlev, nets, nete, qsize;
  int max_corner_elem, state_size, derived_size, elem_size;
  int iam;
} edgeAdvp1_param_t;

typedef struct {
  double *buf;
  double *receive;
  double *addrs[6];
  int *putmap_addr;
  int *getmap_addr;
  int *reverse_addr;
  int *iwesn_addr;
  int np, nlev, nets, nete, qsize, max_corner_elem, strid_size;
  int var_coef, hypervis_power, hypervis_scaling, rrearth;
  int iam;
} edgeq_param_t;
#endif
