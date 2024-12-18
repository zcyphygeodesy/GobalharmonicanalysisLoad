      real*8 function Gauss2D(lon,lat,dt,row,col,hd)
!利用高斯基函数方法内插(lon,lat)处的函数值
!2021年4月20日,章传银
!-------------------------------------------------------------
      implicit none
	integer row,col,ii,jj,ki,kj
	real*8::dt(row,col),hd(6),fun
	real*8::qi,dm,rr,lat,lon,lat0,lon0,RAD,ae
!-------------------------------------------------------------
      RAD=datan(1.d0)/45.d0;ae=6378000
	ii=nint((lat-hd(3))/hd(6)+0.5d0)
	jj=nint((lon-hd(1))/hd(5)+0.5d0)
	qi=0.d0;dm=0.d0;Gauss2D=9999.d0
	do ki=ii-4,ii+4
	  lat0=hd(3)+(dble(ki)-0.5d0)*hd(6)
	  do kj=jj-4,jj+4
	    lon0=hd(1)+(dble(kj)-0.5d0)*hd(5)
	    rr=((lon-lon0)*dcos((lat+lat0)/2.d0*RAD))**2+(lat-lat0)**2
	    if(ki>0 .and. kj>0 .and. ki<=row .and. kj<=col)then
            if(dt(ki,kj)<9000.d0)then
             fun=dexp(-rr/hd(5)**2)
             dm=dm+dt(ki,kj)*fun;qi=qi+fun
            endif
	    endif
	  enddo
	enddo
	if(qi>1.d-18)Gauss2D=dm/qi
9100  continue
      end
!
!******************************************************************
!
      real*8 function CGrdPntD(lon,lat,dt,row,col,hd)
!利用距离反比方法内插(lon,lat)处的函数值
!2005年4月20日,章传银
!-------------------------------------------------------------
      implicit none
	integer row,col,ii,jj,ki,kj
	real*8::dt(row,col),hd(6),mdr
	real*8::qi,dm,rr,lat,lon,lat0,lon0,RAD,ae
!-------------------------------------------------------------
      RAD=datan(1.d0)/45.d0;ae=6378000
	ii=nint((lat-hd(3))/hd(6)+0.5d0)
	jj=nint((lon-hd(1))/hd(5)+0.5d0)
	qi=0.d0;dm=0.d0;CGrdPntD=9999.d0
      mdr=hd(5)*RAD*ae*1.d-3
	do ki=ii-4,ii+4
	  lat0=hd(3)+(dble(ki)-0.5d0)*hd(6)
	  do kj=jj-4,jj+4
	    lon0=hd(1)+(dble(kj)-0.5d0)*hd(5)
	    rr=((lon-lon0)*dcos((lat+lat0)/2.d0*RAD))**2+(lat-lat0)**2
	    rr=dsqrt(rr)*ae+mdr
	    if(ki>0 .and. kj>0 .and. ki<=row .and. kj<=col)then
            if(dt(ki,kj)<9000.d0)then
	        qi=qi+1.d0/rr;dm=dm+dt(ki,kj)/rr
            endif
	    endif
	  enddo
	enddo
	if(qi>1.d-16)CGrdPntD=dm/qi
9100  continue
      end
!
!******************************************************************
!
      real*8 function CGrdPntD2(lon,lat,dt,row,col,hd)
!利用距离平方反比方法内插(lon,lat)处的函数值
!2005年4月20日,章传银
!-------------------------------------------------------------
      implicit none
	integer row,col,ii,jj,ki,kj,kn
	real*8::dt(row,col),hd(6),mdr
	real*8::qi,dm,rr,lat,lon,lat0,lon0,RAD,ae
!-------------------------------------------------------------
      RAD=datan(1.d0)/45.d0;ae=6378000
	ii=nint((lat-hd(3))/hd(6)+0.5d0)
	jj=nint((lon-hd(1))/hd(5)+0.5d0)
	qi=0.d0;dm=0.d0;CGrdPntD2=0.d0
      mdr=hd(5)*RAD*ae*1.d-3
	do ki=ii-4,ii+4
	  lat0=hd(3)+(dble(ki)-0.5d0)*hd(6)
	  do kj=jj-4,jj+4
	    lon0=hd(1)+(dble(kj)-0.5d0)*hd(5)
	    rr=((lon-lon0)*dcos(lat*RAD))**2+(lat-lat0)**2
	    rr=rr*ae**2+mdr**2
	    if(ki>0 .and. kj>0 .and. ki<=row .and. kj<=col)then
            if(dt(ki,kj)<9000.d0)then
	        qi=qi+1.d0/rr;dm=dm+dt(ki,kj)/rr
            endif
          endif
	  enddo
	enddo
	if(qi>1.d-16)CGrdPntD2=dm/qi
9100  continue
      end
!
!******************************************************************
!
      real*8 function CShepard(lon,lat,dt,row,col,hd)
!利用Shepard方法内插(lon,lat)处的函数值
!2005年4月20日,章传银
!-------------------------------------------------------------
      implicit none
	integer row,col,ii,jj,ki,kj
	real*8::dt(row,col),hd(6),mdr,R0,pr
	real*8::qi,dm,rr,lat,lon,lat0,lon0,RAD,ae
!-------------------------------------------------------------
      RAD=datan(1.d0)/45.d0;ae=6378000
	ii=nint((lat-hd(3))/hd(6)+0.5d0)
	jj=nint((lon-hd(1))/hd(5)+0.5d0)
	qi=0.d0;dm=0.d0;CShepard=9999.d0;
      mdr=hd(5)*RAD*ae*1.d-3
	R0=4.d0*hd(5)*RAD*ae
	do ki=ii-4,ii+4
	  lat0=hd(3)+(dble(ki)-0.5d0)*hd(6)
	  do kj=jj-4,jj+4
	    lon0=hd(1)+(dble(kj)-0.5d0)*hd(5)
	    rr=(((lon-lon0)*dcos(lat*RAD))**2+(lat-lat0)**2)*RAD**2
	    rr=dsqrt(rr)*ae+mdr
	    if(ki>0 .and. kj>0 .and. ki<=row .and. kj<=col)then
	      if(rr<=R0/3.d0.and.dt(ki,kj)<9000.d0)then
	        qi=qi+(1.d0/rr)**2;dm=dm+dt(ki,kj)/rr**2
	      endif
	      if(rr>R0/3.d0.and.rr<=R0.and.dt(ki,kj)<9000.d0)then
	        pr=27.d0/4.d0/R0*(rr/R0-1.d0)**2
	        qi=qi+pr**2;dm=dm+dt(ki,kj)*pr**2
	      endif
	    endif
	  enddo
	enddo
	if(qi>1.d-12)CShepard=dm/qi
9100  continue
      end
!
!******************************************************************
!
      real*8 function GrdPnt2Spl(lon,lat,dt,row,col,hd,num)
!利用样条曲面内插(lon,lat)处的函数值
!2005年4月20日,章传银
!-------------------------------------------------------------
      implicit none
	integer row,col,ii,jj,iw,il,ir,jl,jr,num
	real*8::dt(row,col),hd(6),Gauss2D
	real*8::lat,lon,lat0,lon0,dm,val
	real*8,allocatable::H(:,:)
!-------------------------------------------------------------
	  GrdPnt2Spl=Gauss2D(lon,lat,dt,row,col,hd);goto 9001 !算法有bug，用高斯插值代替。2023.5
	ii=nint((lat-hd(3))/hd(6)+0.5d0)
	jj=nint((lon-hd(1))/hd(5)+0.5d0)
	if(ii<num+1.or.ii>row-num-1.or.jj<num+1.or.jj>col-num-1)then
	  GrdPnt2Spl=Gauss2D(lon,lat,dt,row,col,hd);goto 9001
	endif
	allocate(H(2*num+1,2*num+1))
	il=ii-num;ir=ii+num;jl=jj-num;jr=jj+num
	iw=2*num+1;dm=0.d0
	H=dt(il:ir,jl:jr)
	lat0=hd(3)+(dble(il)-.5d0)*hd(6)
	lon0=hd(1)+(dble(jl)-.5d0)*hd(5)
	call InterpSpline(iw,dm,H,lat0,lon0,hd(6),hd(5)
     >	,2*num+1,2*num+1,2*num+1,2*num+1,lat,lon,val)
	GrdPnt2Spl=val
	deallocate(H)
      if(GrdPnt2Spl>9.d3)GrdPnt2Spl=Gauss2D(lon,lat,dt,row,col,hd)
9001  if(GrdPnt2Spl>9.d3)GrdPnt2Spl=9999.d0
      end
