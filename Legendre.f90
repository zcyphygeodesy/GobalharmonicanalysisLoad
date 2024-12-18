      subroutine BelPnm(pnm,maxn,t)
      !计算规格化Pnm
      !规格化Pnm采用改进的Belikov递推算法
	implicit none
	integer::maxn,n,m,k,kk,k1,k2,k3,astat(6)
	real*8::pnm((maxn+2)**2)
	real*8::t,u,a,b,c,d,e
	real*8::anm,bnm,cnm,da
!---------------------------------------------------------------------------
      pnm=0.d0;u=dsqrt(1-t**2);pnm(1)=1.d0
      pnm(2)=dsqrt(3.d0)*t;pnm(3)=dsqrt(3.d0)*u
	do n=1,maxn
        a=dsqrt(2.d0*n+1.d0)/dsqrt(2.d0*n-1.d0)
        b=dsqrt(2.d0*(n-1.d0)*(2.d0*n+1.d0))/dsqrt((2.d0*n-1.d0)*n)
        kk=n*(n+1)/2+1;k1=n*(n-1)/2+1;k2=n*(n-1)/2+2
        pnm(kk)=a*t*pnm(k1)-b*u/2.d0*pnm(k2)
	  do m=1,n
          kk=n*(n+1)/2+m+1;k1=n*(n-1)/2+m+1
          k2=n*(n-1)/2+m+2;k3=n*(n-1)/2+m
          c=a/n*dsqrt(dble(n**2-m**2));d=a/n/2.d0*dsqrt(dble((n-m)*(n-m-1)))
          e=a/n/2.d0*dsqrt(dble((n+m)*(n+m-1)))
          if(m==1)e=e*dsqrt(2.d0)
          pnm(kk)=c*t*pnm(k1)-d*u*pnm(k2)+e*u*pnm(k3)
	  enddo
      enddo
904   return
      end
!
!***************************************************************************
!
      subroutine Pndpn_dt(p,dp,n,t)
	implicit none
	integer::n,m,k
	real*8 ::p(n+5),dp(n+5),t
!---------------------------------------------------------------------------
!	计算全部Legendre函数pn及其导数dpn.p(1)=p0,dp(1)=dp0
      p(1)=1.d0;p(2)=t;p(3)=1.5d0*t**2-0.5d0;p(4)=2.5d0*t**3-1.5d0*t
      dp(1)=0.d0;dp(2)=1.d0;dp(3)=3.d0*t;dp(4)=7.5d0*t**2-1.5d0
      do m=3,n+1
         k=m-1
	   p(m)=real(2*k-1)/real(k)*t*p(m-1)-real(k-1)/real(k)*p(m-2)
	   dp(m)=t*dp(k)+real(m)*p(k)
      enddo
      end
