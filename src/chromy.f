!----------------------------------------------------------------------------------------------------------------------------------
! Subroutine FORTRAN 
! Implementing Chromy algorithm in Bethel.R script
! Version 26 June 2008
!----------------------------------------------------------------------------------------------------------------------------------
      subroutine chromy 
     -           (cost, epsil, nstrat, nvar, maxiter, a, x)
	     integer 
     -		 nstrat, 
     -		 nvar, 
     -       iter,	 
     -       maxiter, 
     -       i, 
     -       j
	 
         double precision 
     -		 a(nstrat,nvar), 
     -	     t1(nvar,nstrat), 
     -   	 t2(nstrat,nvar), 
!     -   	 a1(nstrat,nvar),	 
     -		 x(nstrat),
     -       cost(nstrat), 
     -   	 alfa(nvar), 	
     -	     alfanext(nvar), 
     -	     den1(nstrat), 
     -	     den2s(nstrat), 
     -       c(nvar), 
     -   	 diff(nvar), 
     -       den2, 
     -       alfatot, 
     -	 	 massimo, 
     -	     epsil

      iter = 0
      massimo = dble(0)
      den2 = dble(0)
      alfatot = dble(0)
      do i=1,nstrat
	     do j=1,nvar
            t1(j,i) = dble(0)
            t2(i,j) = dble(0)
            diff(j) = dble(o)
	     end do
      end do	  
      do j=1,nvar
        alfa(j) = 1/dble(nvar)
        alfanext(j) = dble(0)
      end do
      do i=1,nstrat
        den1(i) = dble(0)
        den2s(i) = dble(0)
      end do		  
10    continue
      iter = iter + 1
!----------------------------------------------------------------------------------------------------------------------------------
!              den1 =sqrt(rowSums(t( t(a)*c(alfa)) ))  
      t1 = transpose(a)
      do i=1,nstrat
	     do j=1,nvar
            t1(j,i) = t1(j,i) * alfa(j)
	     end do
      end do
      t2 = transpose (t1)
      do i=1,nstrat
		    den1(i)=0
      end do
      do i=1,nstrat
        do j=1,nvar
		    den1(i)=den1(i)+t2(i,j)
        end do
      end do
      do i=1,nstrat
		    den1(i)=den1(i)**0.5
      end do
!----------------------------------------------------------------------------------------------------------------------------------
!              den2=sum(sqrt(rowSums(t(t(a*cost)*c(alfa)))))	
      do j=1,nvar
	     do i=1,nstrat
            t2(i,j) = a(i,j) * cost(i)
	     end do
      end do
      t1=transpose(t2)
      do i=1,nstrat
	     do j=1,nvar
            t1(j,i) = t1(j,i) * alfa(j)
	     end do
      end do
      t2 = transpose (t1)
      do i=1,nstrat
		    den2s(i)=0
      end do		  
      do i=1,nstrat
	     do j=1,nvar
		    den2s(i)=den2s(i)+t2(i,j)
         end do
      end do
      den2 = 0
      do i=1,nstrat
	     den2 = den2 + sqrt(den2s(i))
      end do
!----------------------------------------------------------------------------------------------------------------------------------
!	     x<-sqrt(cost)/(den1*den2+ epsilon)
      do i=1,nstrat
         x(i) = cost(i)**0.5 / ((den1(i) * den2) + epsil)
      end do
!----------------------------------------------------------------------------------------------------------------------------------	  
!	     alfatot <- sum( c(alfa)*(t(a)%*%x)**2 )   
      t1 = transpose(a)
      do j=1,nvar
          c(j)=0
      end do
      do j=1,nvar
        do i=1,nstrat
          c(j)=c(j)+t1(j,i)*x(i)
        end do
      end do  
      alfatot = 0	  
      do j=1,nvar
        alfatot = alfatot + (alfa(j) * (c(j)**2))
      end do
!----------------------------------------------------------------------------------------------------------------------------------	  
!	     alfanext <- c(alfa)*(t(a)%*%x)**2/alfatot 
      do j=1,nvar
        alfanext(j) = (alfa(j) * c(j)**2) / alfatot
      end do
!----------------------------------------------------------------------------------------------------------------------------------	  
!	     diff <- max(abs(alfanext-alfa))
      massimo = 0
      do j=1,nvar
        diff(j) = abs(alfanext(j)-alfa(j))
		if (diff(j) .gt. massimo) massimo = diff(j) 
      end do
!----------------------------------------------------------------------------------------------------------------------------------	  
!	     alfa <- alfanext  
      do j=1,nvar
        alfa(j) = alfanext(j)
      end do	   
!----------------------------------------------------------------------------------------------------------------------------------	  
!               while ( diff > epsilon && iter<maxiter)
      if (massimo .gt. epsil .and. iter .lt. maxiter) go to 10
      return
      end subroutine	  