#ifndef _REGISTER_H
#define _REGISTER_H
/*get col id of cpes*/
#define GET_COL(col) \
	asm volatile ("rcsr%0, 2" : "=r"(col))
/*get row id of cpes*/
#define GET_ROW(row) \
	asm volatile ("rcsr%0, 1" : "=r"(row))

/*row communication -- send value to "row" cpe*/
#define REG_PUTR(var, dst) \
	asm volatile ("putr %0,%1"::"r"(var),"r"(dst))

/*row communication -- reciver value */
#define REG_GETR(var) \
	asm volatile ("getr %0":"=r"(var))
/*col communication -- send value to "row" cpe*/
#define REG_PUTC(var, dst) \
	asm volatile ("putc %0,%1"::"r"(var),"r"(dst))
/*col communication -- reciver value */
#define REG_GETC(var) \
	asm volatile ("getc %0":"=r"(var))
/* REG_SYNC or REG_SYNR 
 * can be repsented by athread_syn(scope cp,int mask)
 * description of scope & mase
 * typedef enum {
 *	ROW_SCOPE,COL_SCOPE,ARRAY_SCOPE
 * } scope;
 * mask - 0xffff
 * if ROW_SCOPE //low 8 bits is valied for 8 rows
 * if COL_SCOPE //low 8 bits is valied for 8 cols
 * if ARRAY_SCOPE // all bits is valied for 64 cpesl
 * */

/*col synchronization*/
#define REG_SYNR(mask) \
	asm volatile ( "synr%0"::"r"(mask))
/*row synchronization*/
#define REG_SYNC(mask) \
	asm volatile ( "sync%0"::"r"(mask))
#endif 
