      subroutine SphHarmPointFFT(ewh,cilm,nx,nlat,nlon,maxn,hd,GRS)
      !由等效水高球谐系数按FFT算法计算格网等效水高
  	implicit none
	include "fftw3.f"
	integer::nx,maxn,nlat,nlon
	real*8::cilm(nx,maxn+1,maxn+1),ewh(nlat,nlon)
	real*8::GRS(6),hd(6),ae,pi,RAD,dlon,lat,A0
	complex*16::mdln,mdlp,sumn,sump
	complex*16,allocatable::AN(:),AP(:),xN(:),xP(:)
	integer ::n,m,i,j,k
      integer*8::fwd = 0, bwd = 0
	real*8,allocatable::pnm(:)
!---------------------------------------------------------------------
      allocate(pnm((maxn+2)**2))
	pi=dacos(-1.d0);RAD=pi/180.d0;ae=GRS(2);dlon=hd(5)*RAD
      allocate(AN(nlon),AP(nlon),xN(nlon),xP(nlon))
      call dfftw_plan_dft_1d(fwd,nlon,AN,xN,FFTW_FORWARD, FFTW_MEASURE)
      call dfftw_plan_dft_1d(bwd,nlon,AP,xP,FFTW_BACKWARD,FFTW_MEASURE)
      ewh=0.d0
      do i=1,nlat
        lat=hd(3)+(dble(i)-0.5d0)*hd(6)
        call BelPnm(pnm,maxn,dsin(lat*RAD))
        AN=dcmplx(0.d0,0.d0);AP=dcmplx(0.d0,0.d0)
        do m=1,maxn
          mdln=dcmplx(dcos(m*dlon/2.d0),-dsin(m*dlon/2.d0))
          mdlp=dcmplx(dcos(m*dlon/2.d0),dsin(m*dlon/2.d0))
          sumn=dcmplx(0.d0,0.d0);sump=dcmplx(0.d0,0.d0)
          do n=m,maxn
            k=n*(n+1)/2+m+1
            sumn=sumn+dcmplx(cilm(1,n+1,m+1),cilm(2,n+1,m+1))*pnm(k)/2.d0
            sump=sump+dcmplx(cilm(1,n+1,m+1),-cilm(2,n+1,m+1))*pnm(k)/2.d0
          enddo
          AN(m+1)=mdln*sumn;AP(m+1)=mdlp*sump
        enddo
        !对纬度圈等效水高进行FFT变换
        call dfftw_execute(fwd)  !提取cpx(1:nlon)
        call dfftw_execute(bwd)  !提取cpx(1:nlon)
        A0=0.D0
        do n=0,maxn
          k=n*(n+1)/2+1
          A0=A0+cilm(1,n+1,1)*pnm(k)
        enddo
        do j=1,nlon
          ewh(i,j)=dble(xN(j))+dble(xP(j))+A0
        enddo
      enddo
      ewh=ewh*ae
      call dfftw_destroy_plan(fwd) 
      call dfftw_destroy_plan(bwd) 
      deallocate(AN,AP,xN,xP,pnm)
      end
