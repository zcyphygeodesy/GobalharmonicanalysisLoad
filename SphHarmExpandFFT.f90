      subroutine SphHarmExpandFFT(ewh,hd,cilm,nx,maxn,nlat,nlon,GRS)
      !输入等效水高球面格网ewh
  	implicit none
	include 'fftw3.f'
	integer::nx,maxn,nlat,nlon
	real*8::cilm(nx,maxn+1,maxn+1),ewh(nlat,nlon)
	real*8::GRS(6),hd(6),pi,RAD,dlon,b,lat,ae,bm
	complex*16,allocatable::cpx(:),sumL(:,:)
      real*8,allocatable::ewh1(:),bita(:),legI(:,:)
	integer ::n,m,i,j,k
      integer*8::plan  
!---------------------------------------------------------------------
	pi=dacos(-1.d0);RAD=pi/180.d0;ae=GRS(2);dlon=hd(5)*RAD
      allocate(ewh1(nlon),bita(maxn+1),cpx(nlon),legI(maxn+3,maxn+3),sumL(maxn+1,maxn+1))
      call dfftw_plan_dft_r2c_1d(plan,nlon,ewh1,cpx, FFTW_MEASURE)
      sumL=0.d0;cilm=0.d0
      do k=1,nlat
        !对纬度圈等效水高进行FFT变换
        ewh1(1:nlon)=ewh(k,1:nlon)
        call dfftw_execute(plan)  !提取cpx(1:maxn+1)
        !计算平滑因子bita(1:nlon)
        lat=hd(3)+(dble(k)-0.5d0)*hd(6)
        bita=1.d0
        !call smoothfactor(bita,maxn,hd,lat)!不使用平滑因子
        !计算Legendre函数数值积分legI(maxn+1,maxn+1)
        call integralPnm(legI,maxn,hd,lat)
        do n=1,maxn+1
          do m=1,maxn+1
            sumL(n,m)=sumL(n,m)+legI(n,m)*cpx(m)/bita(n) 
          enddo
        enddo
      enddo
      call dfftw_destroy_plan(plan) 
      !计算m维系数Bm*exp(-m*dlon)/4.d0/pi
      do m=1,maxn+1
        if(m==1)b=dlon
        if(m>1)b=2.d0/dble(m-1)*dsin(dble(m-1)*dlon/2.d0)
        Bm=b*exp(-dble(m-1)*dlon)/4.d0/pi
        sumL(1:maxn+1,m)=sumL(1:maxn+1,m)*Bm
      enddo
      do n=1,maxn+1
        do m=1,maxn+1
          cilm(1,n,m)=dble(sumL(n,m))/ae
          cilm(2,n,m)=-dimag(sumL(n,m))/ae
        enddo
      enddo
      deallocate(ewh1,bita,cpx,legI,sumL)
      end
!
!**********************************************************************************
!
      subroutine integralPnm(legI,maxn,hd,lat)
      !计算Legendre函数数值积分legI(maxn,maxn)/(nlon,nlon)
      !计算下三角及对角线，上三角赋零
  	implicit none
	integer ::maxn,n,m,i,j,k,mm
	real*8::hd(6),lat,legI(maxn+3,maxn+3),RAD,dlon,dlat,lat1,lat2
	real*8::u1,u2,t1,t2,u,t
	real*8::anm,anm1,bb,pi,b(10000),xk,y2,y1,tmp(10000)
	real*8,allocatable::pnm1(:),pnm2(:)
!---------------------------------------------------------------------
 	allocate(pnm1((maxn+2)**2),pnm2((maxn+2)**2))
	pi=dacos(-1.d0);RAD=pi/180.d0;legI=0.d0
      dlon=hd(5)*RAD;dlat=hd(6)*RAD
      lat1=lat-hd(6)/2.d0;lat2=lat+hd(6)/2.d0
      u=dcos(lat*RAD);u1=dcos(lat1*RAD);u2=dcos(lat2*RAD)
      t=dsin(lat*RAD);t1=dsin(lat1*RAD);t2=dsin(lat2*RAD)
      call BelPnm(pnm1,maxn,t1);call BelPnm(pnm2,maxn,t2)
      legI(1,1)=t2-t1;legI(2,1)=(t2**2-t1**2)*dsqrt(3.d0)/2.d0
      legI(2,2)=(u2*t2-u1*t1-dacos(t2)+dacos(t1))*dsqrt(3.d0)/2.d0
      do n=2,maxn !向前递推IPnm
        do m=0,n-1
          k=n*(n-1)/2+m+1
          anm=dsqrt(dble(4*n**2-1)/dble(n**2-m**2))
          anm1=dsqrt(dble(4*(n-1)**2-1)/dble((n-1)**2-m**2))
          legI(n+1,m+1)=(dble(n-2)*anm/anm1*legI(n-1,m+1)-anm*(u2**2*pnm2(k)-u1**2*pnm1(k)))/dble(n+1)
        enddo
      enddo
      !向前递推IPnn
      b(1)=dsqrt(3.d0)
      do n=2,maxn+2
        b(n)=dsqrt(dble(2*n+1)/dble(2*n))
      enddo
      do n=2,maxn
        k=n*(n+1)/2+n+1
        legI(n+1,n+1)=(dble(n)*b(n)*b(n-1)*legI(n-1,n-1)+(t2*pnm2(k)-t1*pnm1(k)))/dble(n+1)
      enddo
      if(dabs(lat)>80.d0)then !向后递推IPnn低阶部分
        bb=b(1);mm=1+nint((1.d0+dlog(1.d-20))/dlog((dcos(lat*RAD))**2))
        do n=2,maxn+1
          bb=bb*b(n)
        enddo
        n=maxn+1
        y1=0.d0;y2=0.d0
        do i=1,mm
          if(i==1)xk=1.d0
          if(i>1)xk=xk*dble(2*i-3)/dble(2*(i-1))
          y1=y1+xk/dble(n+2*i)*u1**(2*(i-1))
          y2=y2+xk/dble(n+2*i)*u2**(2*(i-1))
        enddo
        tmp(n+1)=-bb*u**(n+2)*(y2-y1)
        n=maxn+2;bb=bb*b(maxn+2)
        y1=0.d0;y2=0.d0
        do i=1,mm
          if(i==1)xk=1.d0
          if(i>1)xk=xk*dble(2*i-3)/dble(2*(i-1))
          y1=y1+xk/dble(n+2*i)*u1**(2*(i-1))
          y2=y2+xk/dble(n+2*i)*u2**(2*(i-1))
        enddo
        tmp(n+1)=-bb*u**(n+2)*(y2-y1)
        do n=maxn,2,-1
          k=(n+2)*(n+3)/2+n+3
          tmp(n+1)=(dble(n+3)*tmp(n+3)-t2*pnm2(k)+t1*pnm1(k))/dble(n+2)/b(n+2)/b(n+1)
        enddo
        do n=2,maxn/2
          legI(n+1,n+1)=tmp(n+1)
        enddo
      endif
3001  continue
      deallocate(pnm1,pnm2)
      end
!
!**********************************************************************************
!
      subroutine smoothfactor(bita,maxn,hd,lat)
      !计算平滑因子bita(nlon)
  	implicit none
	integer::maxn,i
	real*8::hd(6),lat,bita(maxn+1),RAD,dlon,dlat,lat1,lat2
	real*8::p(maxn+5),dp(maxn+5),cosf0,pi
!---------------------------------------------------------------------
	pi=dacos(-1.d0);RAD=pi/180.d0
      dlon=hd(5)*RAD
      lat1=(lat-hd(6)/2.d0)*RAD;lat2=(lat+hd(6)/2.d0)*RAD
      cosf0=(dsin(lat1)-dsin(lat2))*dlon/2.d0/pi+1.d0!!!dabs
      call Pndpn_dt(p,dp,maxn+3,cosf0)
      do i=1,maxn+1
        bita(i)=(p(i)-p(i+2))/dble(2*i+1)/(1.d0-cosf0)
        if(i<maxn/3)bita(i)=bita(i)**2
      enddo
      end
