#include<math.h>
#define rex(x,y,z) (x+(y/(x*x)-x)*z)


//#define rtc(_x) asm volatile("rcsr %0, 4" : "=r"(_x))

double crt(double x)
{
   double old;
   double zero=0.,new=1.,one=1.,eig=8.,half=0.5,deig=0.125,three=3.,dt3=one/three;
   int i;
   int count=0;
   
   if(x==0.) return 0.;

   if(x<0.) 
    { 
      x=-x;
      if(x<one)
        { 
          while(x<deig) {x*=eig;count++;}
          for(i=0;i<7;i++)  
            { 
              old=new;
              new=rex(old,x,dt3);
            }
          for(i=0;i<count;i++) new*=half;
        }
      else
        {
          while(x>eig) {x*=deig;count++;}
          for(i=0;i<7;i++)  
            { 
              old=new;
              new=rex(old,x,dt3);
            }
          for(i=0;i<count;i++) new+=new;
        }
      new=-new;
      return(new);
    }
   else 
    {
     if(x<one)
        {
          while(x<deig) {x*=eig;count++;}
          for(i=0;i<7;i++)
            {
              old=new;
              new=rex(old,x,dt3);
            }
          for(i=0;i<count;i++) new*=half;
        }
      else
        {
          while(x>eig) {x*=deig;count++;}
          for(i=0;i<7;i++)
            {
              old=new;
              new=rex(old,x,dt3);
            }
          for(i=0;i<count;i++) new+=new;
        }
      return(new);
    }

}

