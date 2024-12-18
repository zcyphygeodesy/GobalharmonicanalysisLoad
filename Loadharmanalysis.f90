      subroutine Loadharmanalysis(dtmfl,loadfl,knd,dk,itd)
      implicit none
	character*800::infl,outfl,chfl,prfl,gcfl
	character*800::line,dtmfl,loadfl
	integer i,j,nlon,nlat,nlon0,nlat0,sn,kln,kk,nk,astat(8),maxn,nd,m,n,knd,pm,IY,IM,ID,ih,imn
	real*8 rec(800),hd(6),hd0(6),tmp,tm,tdate,yr0,mjd,mjd0,fk,dx(3),sec
	real*8::GRS(6),ae,rst0(4),rst(4),dk,itd,cm
      real*8::CGrdPntD2,BLH(3),err(40),bq,dp(40)!dp(4)迭代标准差
	real*8,allocatable::ewh(:,:),ewh0(:,:),ewh1(:,:),dtm(:,:),zero(:,:)
	real*8,allocatable::cilm(:,:,:),cilm0(:,:,:)
	integer::status=0
!---------------------------------------------------------------------
      GRS(1)= 3.986004415d14; GRS(2)=6378137.d0; GRS(3)=1.0826359d-3
      GRS(4) = 7.292115d-5; GRS(5)=1.d0/298.25641153d0;ae=GRS(2)
      fk=dsqrt(3.d0)*1.025d0/5.517d0*ae
      !输入迭代终止条件：残差标准差为原格值标准的dk倍
      !input the iteration condition parameters
      !Iteration termination condition: The standard deviation of the residual grid value is less than dk of the standard deviation
      !of the original grid value, or the difference of the residual standard deviation of the previous step iteration relative to 
      !the current step iteration is less than itd of the standard deviation of the original grid values.
      dk=1.d-3; itd=1.d-4
      !knd=0 landwater, =1 sea level variation =-1 surface atmosphere variation
      knd=1!knd=0陆地水，=1海平面变化，=-1地面大气压
	BLH(3)=0.d0;nk=0 !nk为程序批量构造的负荷球谐系数模型文件数
      !nk - file number of load spherical harmonic coefficient model time series
      cm=1.d-2;if(knd==1)cm=1.03d-2
      if(knd==-1)cm=1.0204d-2!大气压cm=hpa
      !打开陆海地形球坐标格网模型，用于分离陆海区域knd=0 or =1
      !The spatial resolution of the land-sea terrain grid should not be lower than that of the surface load grid.
      open(unit=8,file=dtmfl,status="old",iostat=status)
      if(status/=0)goto 904
      read(8,'(a)') line
      call PickRecord(line,kln,rec,sn)
      if(sn<6)then
          close(8);goto 904
      endif
	hd0(1:6)=rec(1:6)
	nlat0=nint((hd0(4)-hd0(3))/hd0(6))
	nlon0=nint((hd0(2)-hd0(1))/hd0(5))
	hd0(5)=(hd0(2)-hd0(1))/dble(nlon0)
	hd0(6)=(hd0(4)-hd0(3))/dble(nlat0)
	allocate(dtm(nlat0,nlon0), stat=astat(1))
	if (sum(astat(1:1)) /= 0) then
         close(8);goto 904
	endif
	do i=1,nlat0
	  read(8,*,end=705)(dtm(i,j),j=1,nlon0)
      enddo
705   close(8)
      open(unit=40,file='geocentricvariation.txt',status="replace")
      write(40,'(a12,3a10,a14)')"Year","Xcm","Ycm","Zcm","date"
      !打开球坐标格网时间序列文件名文件,球坐标格网头文件第七个数为长整型ETideLoad系统格式时间
      !The seventh numerical value of spherical coordinate load ewh grid file header is the long integer time agreed by ETideLoad
      open(unit=58,file=loadfl,status="old",iostat=status)
      if(status/=0) goto 902 
	do while(.not.eof(58))
        read(58,'(a)') infl
        open(unit=8,file=infl,status="old",iostat=status)
        if(status/=0)goto 901
        read(8,'(a)') line
        call PickRecord(line,kln,rec,sn)
        !pick the long integer time
        call PickRecord(line,kln,rec,sn)
        if(sn<7)then
          close(8);goto 901
        endif
	  tm=rec(7)
        call tmcnt(tm,IY,IM,ID,ih,imn,sec)
        call CAL2JD (IY,IM,ID,mjd,j)
        mjd=mjd+dble(ih)/24.d0+dble(imn)/1440.d0+dble(sec)/864.d2 !GPS_MJD
        if(nk==0)then
          call CAL2JD (IY,1,1,mjd0,j)
        endif
        tdate=dble(IY)+(mjd-mjd0)/365.25d0
        if(sn<6)goto 901
	  hd(1:6)=rec(1:6)
	  nlat=nint((hd(4)-hd(3))/hd(6))
	  nlon=nint((hd(2)-hd(1))/hd(5))
	  hd(5)=(hd(2)-hd(1))/dble(nlon)
	  hd(6)=(hd(4)-hd(3))/dble(nlat)
        maxn=nlat
        pm=nint(10.d0/hd(6))
	  if (maxn>2160) then
           close(8);goto 901
	  endif
        if(nk>0)goto 1111
	  allocate(ewh(nlat,nlon), stat=astat(1))
	  allocate(ewh0(nlat,nlon), stat=astat(2))
	  allocate(ewh1(nlat,nlon), stat=astat(3))
	  allocate(cilm(2,maxn+1,maxn+1), stat=astat(4))
	  allocate(cilm0(2,maxn+1,maxn+1), stat=astat(5))
	  allocate(zero(nlat,nlon), stat=astat(6))
	  if (sum(astat(1:6)) /= 0) then
           close(8);goto 901
	  endif
1111    ewh0=0.d0;ewh1=0.d0;zero=1.d0
	  do i=1,nlat
	    read(8,*,end=905)(ewh(i,j),j=1,nlon)
          BLH(1)=hd(3)+(dble(i)-0.5d0)*hd(6)
          do j=1,nlon
	      BLH(2)=hd(1)+(dble(j)-0.5d0)*hd(5)
            if(ewh(i,j)>9900.0)ewh(i,j)=0.d0
            if(knd>-0.5)then 
              tmp=CGrdPntD2(BLH(2),BLH(1),dtm,nlat0,nlon0,hd0)
              if(knd==0.and.tmp<0.d0)then!zero in sea for landwater variation
                zero(i,j)=0.d0; ewh(i,j)=0.d0
              endif
              if(knd==1.and.tmp>0.d0)then!zero on land for sea level variation
                zero(i,j)=0.d0; ewh(i,j)=0.d0
              endif
            endif
            ewh(i,j)=ewh(i,j)*cm!cm/hPa->m
          enddo
	  enddo
905	  close(8)
        call StatlsGrd(ewh,zero,nlat,nlon,pm,rst0)
        bq=rst0(2);ewh0=ewh;cilm=0;nd=0
        dp(1)=rst0(2);kk=1
	  write(prfl,*)infl(1:len_trim(infl)-4),"pro.ini" 
	  write(outfl,*)infl(1:len_trim(infl)-4),"cs.dat"
	  write(chfl,*)infl(1:len_trim(infl)-4),"rnt.dat"
        open(unit=14,file=prfl,status="replace")
        write(14,102)nd,(rst0(j)/cm,j=1,4)
        !采用FFT方法计算地面等效水高球谐系数
        !calculate surface load spherical harmonic coefficients by FFT method.
908     call SphHarmExpandFFT(ewh0,hd,cilm0,2,maxn,nlat,nlon,GRS)
        call SphHarmPointFFT(ewh1,cilm0,2,nlat,nlon,maxn,hd,GRS)
        ewh1=ewh0-ewh1 !计算残差等效水高格网 residual EWH
        call StatlsGrd(ewh1,zero,nlat,nlon,pm,rst)
        write(14,102)nd+1,(rst(j)/cm,j=1,4)
        err(nd+1)=rst(2);bq=dp(kk);kk=kk+1;dp(kk)=rst(2)
        if(rst(2)>rst0(2)*dk.and.bq-rst(2)>rst0(2)*itd.and.dp(kk)<dp(kk-1).and.nd<20)then  !!!!!迭代条件mnt
           cilm=cilm+cilm0;ewh0=ewh1;nd=nd+1;goto 908
        endif
        if(dp(kk)<dp(kk-1))then
          cilm=cilm+cilm0;ewh0=ewh1
        endif
        call StatlsGrd(ewh0,zero,nlat,nlon,pm,rst)
        close(14)
        rst(1)=cilm(1,1,1)*ae/cm!zero-degree term 零阶项
        open(unit=10,file=outfl,status="replace")
        write(10,'(F13.9,F13.2,F12.4,F8.3)')GRS(1)*1.d-14,GRS(2),rst(1),rst(2)/rst0(2)*1.d2
 	  do n=1,maxn
	    do m=0,n
	      write(10,'(2I6,2ES24.16)')n,m,cilm(1,n+1,m+1),cilm(2,n+1,m+1)
	    enddo
	  enddo
        close(10)
        open(unit=12,file=chfl,status="replace")
        write(12,101)(hd(i),i=1,6)
	  do i=1,nlat
	    write(12,'(15F12.4)')(ewh0(i,j)/cm,j=1,nlon)!m->cm/hPa
	  enddo
        close(12)
        !kk-1迭代次数
        dx(1)=cilm(1,2,2); dx(2)=cilm(2,2,2);dx(3)=cilm(1,2,1)
        dx=dx*fk*1.d3
        write(40,'(F12.4,3F10.4,I14)')tdate,dx(1),dx(2),dx(3),nint(tm)
        write(*, 104)nint(tm),rst(1),rst(2)/rst0(2)*1.d2
        nk=nk+1
901     continue
      enddo
	if(nk>0)deallocate(ewh,ewh0,ewh1,cilm,cilm0,zero)
902   continue
      close(40);close(58)
704   deallocate(dtm)
904   continue
101   format(4F6.1,2F13.8,F14.1)
102   format(I3,4F12.4)
104   format('     Epoch time',I12,':  zero-degree term',f10.4,',  relative error',f8.3,'%')
      end
!
!******************************************************************
      subroutine StatlsGrd(grid,zero,row,col,pm,rst)
!-------------------------------------------------------------
      implicit none
	integer::row,col,i,j,kk,pm
	real*8::grid(row,col),zero(row,col),rst(4),maxt,mint,pv,err,fr
!-----------------------------------------------------------------------
	pv=0.d0;err=0.d0;maxt=-9.d28;mint=9.d28
      fr=dble(row*col);kk=0
      do i=pm+1,row-pm
        do j=1,col
          if(grid(i,j)>9000.d0.or.zero(i,j)<0.5d0)goto 1001
          kk=kk+1;pv=pv+grid(i,j)/fr
	    if(maxt<grid(i,j))maxt=grid(i,j)
          if(mint>grid(i,j))mint=grid(i,j)
1001      continue
        enddo
      enddo
      pv=pv/dble(kk)*fr
      do i=1,row
        do j=1,col
          if(grid(i,j)>9000.d0.or.zero(i,j)<0.5d0)goto 1002
          err=err+(grid(i,j)-pv)**2/dble(kk)
1002      continue
        enddo
      enddo
      err=dsqrt(err)
	rst(1)=pv; rst(2)=err; rst(3)=mint;rst(4)=maxt
	return
      end
