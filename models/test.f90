program test


     integer i,j
     integer b(3,3)


     do 100 i=1,10
          if(i<5) then
              write(*,*), "chenyuhu"
              go to 100
          endif


         write(*,*) "zhengchen"

    100 continue
    do j=1,3
      do i=1,3
        b(i,j)=i+j+i
      enddo
   enddo
   write(*,*),b
   write(*,*),b(1,:)
   


end program test
