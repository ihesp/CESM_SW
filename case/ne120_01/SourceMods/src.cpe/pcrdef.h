//从核性能计数器0的计数事件选择。无复位值。具体定义如下：
PCR_DEF_BEGIN(PCR0)
PCR_OPT_DEF(PC0_CYCLE             , 0x00) //5'h00:   周期计数
PCR_OPT_DEF(PC0_CNT_COND_JMP      , 0x01) //5'h01:   执行站台执行的条件转移指令计数
PCR_OPT_DEF(PC0_CNT_JMP           , 0x02) //5'h02:   执行站台执行的转移指令计数
PCR_OPT_DEF(PC0_CNT_INST          , 0x03) //5'h03:   执行指令计数
PCR_OPT_DEF(PC0_CNT_PIPE0_INST    , 0x04) //5'h04:   流水线0发射指令计数
PCR_OPT_DEF(PC0_CYC_LD_ST_BLOCK   , 0x05) //5'h05:   因LD/ST缓冲满而导致流水线停顿的周期计数
PCR_OPT_DEF(PC0_RESERVED0         , 0x06) //保留
PCR_OPT_DEF(PC0_RESERVED1         , 0x07) //保留
PCR_OPT_DEF(PC0_CYC_GLQ_BLOCK     , 0x08) //5'h08:   因GLQ队列满而导致流水线停顿的周期计数
PCR_OPT_DEF(PC0_CYC_GSQ_BLOCK     , 0x09) //5'h09:   因GSQ队列满而导致流水线停顿的周期计数
PCR_OPT_DEF(PC0_CYC_SBMD_BLOCK    , 0x0a) //5'h0a:   因SBMD接收缓冲满而导致流水线停顿的周期计数
PCR_OPT_DEF(PC0_CYC_CHNL_BLOCK    , 0x0b) //5'h0b:   因通道缓冲满而导致流水线停顿的周期计数
PCR_OPT_DEF(PC0_CNT_LDM_READ      , 0x0c) //5'h0c:   执行站台发起的LDM读访问次数计数
PCR_OPT_DEF(PC0_CNT_LDM_WRITE     , 0x0d) //5'h0d:   执行站台发起的LDM写访问次数计数
PCR_OPT_DEF(PC0_RESERVED2         , 0x0e) //保留
PCR_OPT_DEF(PC0_CNT_ICACHE_ACCESS , 0x0f) //5'h0f:   IBOX对L1 ICache的访问次数计数
PCR_OPT_DEF(PC0_CYC_MBOX_LDM_BLOCK, 0x10) //5'h10:   MBOX访问LDM的请求被阻塞的周期计数
PCR_OPT_DEF(PC0_CNT_JMP_PRED_FAIL , 0x11) //5'h11:   执行站台判断出JMP指令预测失败的次数
PCR_OPT_DEF(PC0_CNT_2_LDM_SCHED   , 0x12) //5'h12:   LDM同时仲裁上两个（流水线和核外）访问LDM请求的次数计数
PCR_OPT_DEF(PC0_CNT_2DATA_1ACCESS , 0x13) //5'h13:   LDM两个数据体一次只有一个被访问的次数计数
PCR_OPT_DEF(PC0_CNT_ADD_SUB_MUL   , 0x14) //5'h14:   浮点加减乘、乘加类指令等效操作计数（标量加减乘+1，标量乘加+2，向量加减乘+4，向量乘加+8）
PCR_DEF_END(PCR0)

//从核性能计数器1的计数事件选择。无复位值。具体定义如下：
PCR_DEF_BEGIN(PCR1)
PCR_OPT_DEF(PC1_CYCLE             , 0x00) //5'h00:   周期计数
PCR_OPT_DEF(PC1_CNT_JMP_FAIL      , 0x01) //5'h01:   执行站台判断出条件转移指令失败的次数
PCR_OPT_DEF(PC1_CNT_COND_JMP      , 0x02) //5'h02:   执行站台执行的条件转移指令计数
PCR_OPT_DEF(PC1_CNT_NOCOND_JMP    , 0x03) //5'h03:   执行站台执行的无条件转移指令计数
PCR_OPT_DEF(PC1_CNT_JMP           , 0x04) //5'h04:   执行站台执行的JMP指令计数
PCR_OPT_DEF(PC1_CNT_PIPE1_INST    , 0x05) //5'h05:   流水线1发射指令计数
PCR_OPT_DEF(PC1_CYC_SYN           , 0x06) //5'h06:   同步指令导致流水线停顿的周期计数
PCR_OPT_DEF(PC1_CYC_ACCESS_BLOCK  , 0x07) //5'h07:   因访存指令未处理完而导致流水线停顿的周期计数
PCR_OPT_DEF(PC1_CYC_PUT_BLOCK     , 0x08) //5'h08:   因PUT缓冲满导致流水线停顿的周期计数
PCR_OPT_DEF(PC1_RESERVED0         , 0x09) //保留
PCR_OPT_DEF(PC1_CYC_GETR_BLOCK    , 0x0a) //5'h0a:   因GETR缓冲空导致流水线停顿的周期计数
PCR_OPT_DEF(PC1_CYC_GETC_BLOCK    , 0x0b) //5'h0b:   因GETC缓冲空导致流水线停顿的周期计数
PCR_OPT_DEF(PC1_CYC_1INST         , 0x0c) //5'h0c:   发射站不能发射两条指令的周期计数
PCR_OPT_DEF(PC1_CYC_DATA_REL_BLOCK, 0x0d) //5'h0d:   数据相关导致流水线停顿的周期计数
PCR_OPT_DEF(PC1_CYC_SYN_JMP_BLOCK , 0x0e) //5'h0e:   因同步指令和条件转移指令导致流水线停顿的周期计数
PCR_OPT_DEF(PC1_CNT_ICACHE_MISS   , 0x0f) //5'h0f:   L1 ICache脱靶次数计数
PCR_OPT_DEF(PC1_CNT_LDM_ACCESS    , 0x10) //5'h10:   执行站台发起的LDM访问次数
PCR_OPT_DEF(PC1_CNT_TBOX_BLK_MBOX , 0x11) //5'h11:   MBOX访问LDM的请求被TBOX普通请求阻塞的周期计数
PCR_OPT_DEF(PC1_CNT_ATOM_BLK_MBOX , 0x12) //5'h12:   MBOX访问LDM的请求被原子修改写请求阻塞的周期计数
PCR_OPT_DEF(PC1_CNT_ADDR_BLK_MBOX , 0x13) //5'h13:   MBOX访问LDM的请求因为地址相关阻塞的周期计数
PCR_OPT_DEF(PC1_CNT_2DATA_2ACCESS , 0x14) //5'h14:   LDM两个数据体一次被同时访问的次数计数
PCR_OPT_DEF(PC1_CNT_DIV_SQRT      , 0x15) //5'h15:   浮点除法、平方根类指令等效操作计数（标量指令+1，向量指令+4）
PCR_DEF_END(PCR1)

//TA性能计数器的性能计数事件选择
PCR_DEF_BEGIN(PCR2)
PCR_OPT_DEF(PC2_CNT_CTRL_REQ      , 0x00) //5’h00:   阵列控制网络请求总数
PCR_OPT_DEF(PC2_CNT_GLD           , 0x01) //5’h01:   gld请求计数
PCR_OPT_DEF(PC2_CNT_GST           , 0x02) //5’h02:   gst请求计数
PCR_OPT_DEF(PC2_CNT_GF_AND_A      , 0x03) //5’h03:   gf&a请求计数
PCR_OPT_DEF(PC2_CNT_GUPDATE       , 0x04) //5’h04:   gupdate请求计数
PCR_OPT_DEF(PC2_CNT_ICACHE_MISS   , 0x05) //5’h05:   L1 ICache自主脱靶次数（由于有SBMD模式，该统计值与从核核心的统计值可能有所不同）
PCR_OPT_DEF(PC2_CNT_ROW_SYN       , 0x06) //5’h06:   行同步请求计数
PCR_OPT_DEF(PC2_CNT_COL_SYN       , 0x07) //5’h07:   列同步请求计数
PCR_OPT_DEF(PC2_CNT_USER_INT      , 0x08) //5’h08:   用户中断请求计数
PCR_OPT_DEF(PC2_CNT_DMA_REQ       , 0x09) //5’h09:   发起DMA请求计数
PCR_OPT_DEF(PC2_CNT_SBMD_START    , 0x0a) //5’h0a:   sbmd_start命令计数
PCR_OPT_DEF(PC2_CNT_SBMD_END      , 0x0b) //5’h0b:   sbmd_end命令计数
PCR_OPT_DEF(PC2_CNT_SBMD_BREAK    , 0x0c) //5’h0c:   sbmd_break命令计数
PCR_OPT_DEF(PC2_CYC_PIPE_MISS     , 0x0d) //5’h0d:   从核流水线脱靶等待的周期计数
PCR_OPT_DEF(PC2_CYC_SBMD          , 0x0e) //5’h0e:   SBMD状态周期计数
PCR_OPT_DEF(PC2_CYC_SELF          , 0x0f) //5’h0f:   自主运行状态周期计数
PCR_OPT_DEF(PC2_CNT_REG_PUT       , 0x10) //5’h10:   寄存器通信“PUT”请求包总数
PCR_OPT_DEF(PC2_CNT_REG_PUTR      , 0x11) //5’h11:   寄存器通信行向“PUT”请求包数
PCR_OPT_DEF(PC2_CNT_REG_PUTC      , 0x12) //5’h12:   寄存器通信列向“PUT”请求包数
PCR_OPT_DEF(PC2_CNT_FORE_LDM_READ , 0x13) //5’h13:   外部对LDM读总数
PCR_OPT_DEF(PC2_CNT_FORE_LDM_WRITE, 0x14) //5’h14:   外部对LDM写总数
PCR_OPT_DEF(PC2_CNT_REPLY_INC     , 0x15) //5’h15:   回答字自增1总数
PCR_OPT_DEF(PC2_CNT_INST_LOAD     , 0x16) //5’h16:   指令装填总数（16B为单位）
PCR_OPT_DEF(PC2_CNT_INST_LOAD_MISS, 0x17) //5’h17:   指令自主脱靶装填总数（16B为单位）
PCR_OPT_DEF(PC2_CNT_INST_LOAD_SBMD, 0x18) //5’h18:   SBMD指令装填总数（16B为单位）
PCR_OPT_DEF(PC2_CNT_MEM_READ_RESP , 0x19) //5’h19:   主存读响应总数
PCR_OPT_DEF(PC2_CNT_MEM_WRITE_RESP, 0x1a) //5’h1a:   主存写响应总数
PCR_OPT_DEF(PC2_CNT_OUTER_READ    , 0x1b) //5’h1b:   外部IO读总数（不含LDM访问）
PCR_OPT_DEF(PC2_CNT_OUTER_WRITE   , 0x1c) //5’h1c:   外部IO写总数（不含LDM访问）
PCR_OPT_DEF(PC2_CNT_REG_LDRC      , 0x1d) //5’h1d:   寄存器通信“ld&发送”请求包总数
PCR_OPT_DEF(PC2_CNT_REG_LDR       , 0x1e) //5’h1e:   寄存器通信行向“ld&发送”请求包数
PCR_OPT_DEF(PC2_CNT_REG_LDC       , 0x1f) //5’h1f:   寄存器通信列向“ld&发送”请求包数
PCR_DEF_END(PCR2)
