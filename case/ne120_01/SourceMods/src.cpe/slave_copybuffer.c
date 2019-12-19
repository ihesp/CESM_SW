#include<slave.h>
typedef struct{
	double*a,*b;
	int length,length1;
} str_par;
void slave_copybuffer_(str_par* p){
	/*get pointer from fortran*/
	double* ptr_a=p->a;
	double* ptr_b=p->b;
	/*get threads id*/
	volatile unsigned int id=athread_get_id(-1);
	volatile unsigned long get_reply,put_reply;
	/*ldm buffer*/
	int N=p->length,i,jump;
	jump=N;
	if(id==63) N=p->length+p->length1;
	double a[N],b[N];
	get_reply=0;
	athread_get(PE_MODE,ptr_a+id*jump,
			a,sizeof(double)*N,&get_reply,0,0,0);
	athread_get(PE_MODE,ptr_b+id*jump,
			b,sizeof(double)*N,&get_reply,0,0,0);
	while(get_reply!=2);
	for(i=0;i<N;i++){
		b[i]=a[i];
	}
	put_reply=0;
	athread_put(PE_MODE,a,ptr_a+id*jump,
			sizeof(double)*N,&put_reply,0,0); 
	athread_put(PE_MODE,b,ptr_b+id*jump,
			sizeof(double)*N,&put_reply,0,0); 
	while(put_reply!=2);
}
