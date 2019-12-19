#!/bin/csh

#modify the source code directly
#csm_share/include/dynamic_vector_typedef.inc
#csm_share/shr/CMakeLists.txt

set SrcDIR = /home/export/online1/cesm01/wangl/cesm1_3_beta17_sehires02_asdphys/models/csm_share/shr

foreach file (shr_spfn_mod.F90 shr_sys_mod.F90 shr_infnan_mod.F90 shr_isnan.h  shr_assert_mod.F90 shr_strconvert_mod.F90 shr_log_mod.F90 shr_kind_mod.F90 shr_stream_mod.F90 shr_cal_mod.F90 )
cp $SrcDIR/$file .
chmod 750 $file
end



#  cesm1.3文件修改list
#  csm_share/shr/shr_spfn_mod.F90                     简单语法修改                  
#  csm_share/shr/shr_sys_mod.F90                      从1.2.1拷贝
#  csm_share/shr/CMakeLists.txt                       去掉里面关于shr_infan_mod.F90.in
#                                                     改用直接使用shr_infan_mod.F90 （1.2.1拷贝）
#  												   新文件缺少几个定义，定义需要参照1.2.1的
#  												   shr_test_infnan_mod.F90给出的定义
#  csm_share/shr/shr_infnan_mod.F90.in                不再使用
#  csm_share/shr/shr_infnan_mod.F90                   追加：1.2.1拷贝
#  csm_share/shr/shr_isnan.h                          追加：1.2.1拷贝
#  csm_share/shr/shr_assert_mod.F90                   追加：1.2.1拷贝               
#  csm_share/shr/shr_strconvert_mod.F90               f2003语法修改
#  csm_share/shr/shr_log_mod.F90                      f2003语法修改
#  csm_share/shr/shr_kind_mod.F90                     把1.2.1里的shr_test_infnan_mod.F90追加在这个文件后面
#  csm_share/shr/shr_stream_mod.F90                   使用的dynamic_vector_typedef.inc，包含很多f2003写法，
#                                                     暂时使用1.2.1的实现方式，若include文件修改完成后可以考虑替换回去
#  												   后面大气部分也遇到同样问题，且1.2.1中无同名文件
#  csm_share/shr/shr_cal_mod.F90                      简单语法修改，去掉了一个变量的F2003可定义属性
#  
#  csm_share/include/dynamic_vector_typedef.inc       语法修改（？） 尚未完成
#  
#  ./atm/cam/src/physics/cam/micro_mg_data.F90        语法修改F2003（？）
#  ./atm/cam/src/physics/cam/micro_mg_utils.F90       语法修改F2003（？）
#  
#  models/utils/esmf_wrf_timemgr                       esmf相关文件手动编译，然后拷贝，
#                                                      猜测MCT/noesmf/a1l1r1i1o1g1w1/csm_share/Depends
#                                                      这个文件是由cmake生成，依赖关系里面用到.mod文件，而我们编译器下
#                                                      所有.mod文件全为大写，依赖关系不能正确反映
#                                                      修改这个文件应该可以不需要全手动实现													
