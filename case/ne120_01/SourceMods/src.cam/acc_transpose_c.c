/*************************************************************************
	> File Name: ../src.cam/acc_transpose_c.c
	> Author: Xu Kai
	> Created Time: 2018年12月17日 星期一 14时50分37秒
 ************************************************************************/
	
	
#include<stdio.h>
#include<stdlib.h>

void acc_transpose_c_(double *src, double *dst, int *dim1_, int *dim2_)
{
    int shape_ori[2] = { *dim1_, *dim2_ };
    int swap_vec[2] = { 2, 1 };
    acc_master_transpose(src, dst, 8, 2, shape_ori, swap_vec);
}
void acc_transpose_c_v2_(double *src, double *dst, int *dim1_, int *dim2_, int *dim3_)
{
    int shape_ori[3] = { *dim1_, *dim2_, *dim3_ };
    int swap_vec[3] = { 2, 1, 3 };
    acc_master_transpose(src, dst, 8, 3, shape_ori, swap_vec);
}
