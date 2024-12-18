      subroutine TimeJDday(tim,sgn)
!2006年10月1日,章传银
!时间相互转换
!-------------------------------------------------------------
      implicit none
	integer tim(4),year,month,day,DJ,sgn,j
	real*8 DJ0,DJ1,fd
!-----------------------------------------------------------------------
      if(sgn==0)then
	  year=tim(1);month=tim(2);day=tim(3)
	  call CAL2JD(year,month,day,DJ1,j)
	  tim(4)=nint(DJ1)
	  goto 9001
	endif
	DJ1=dble(tim(4))
      call JD2CAL(2400000.5D0,DJ1,year,month,day,FD,j)
	tim(1)=year;tim(2)=month;tim(3)=day
9001	return
      end
!
!***********************************************************************
!
	SUBROUTINE CAL2JD (IY0,IM0,ID0,DJM,J)
*  - - - - - - - - - - -
*   i a u _ C A L 2 J D
*  - - - - - - - - - - -
*  Gregorian Calendar to Julian Date.
*
*  Given:
*     IY,IM,ID    i     year, month, day in Gregorian calendar (Note 1)
*
*  Returned:
*     DJM0        d     MJD zero-point: always 2400000.5
*     DJM         d     Modified Julian Date for 0 hrs
*     J           i     status:
*                           0 = OK
*                          -1 = bad year   (Note 3: JD not computed)
*                          -2 = bad month  (JD not computed)
*                          -3 = bad day    (JD computed)
*-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER IY0, IM0, ID0,IY, IM, ID
      real*8 DJM0, DJM
      INTEGER J, MY, IYPMY

*  Earliest year allowed (4800BC)
      INTEGER IYMIN
      PARAMETER ( IYMIN = -4799 )

*  Month lengths in days
      INTEGER MTAB(12)
      DATA MTAB / 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IY=IY0;IM=IM0;ID=ID0
*  Preset status.
      J = 0

*  Validate year.
      IF ( IY.LT.IYMIN ) THEN
         J = -1
      ELSE

*     Validate month.
         IF ( IM.GE.1 .AND. IM.LE.12 ) THEN

*        Allow for leap year.
            IF ( MOD(IY,4) .EQ. 0 ) THEN
               MTAB(2) = 29
            ELSE
               MTAB(2) = 28
            END IF
            IF ( MOD(IY,100).EQ.0 .AND. MOD(IY,400).NE.0 ) MTAB(2) = 28

*        Validate day.
            IF ( ID.LT.1 .OR. ID.GT.MTAB(IM) ) J = -3

*        Result.
            MY = ( IM - 14 ) / 12
            IYPMY = IY + MY
            DJM0 = 2400000.5D0
            DJM = DBLE( ( 1461 * ( IYPMY + 4800 ) ) / 4
     :                + (  367 * ( IM-2 - 12*MY ) ) / 12
     :                - (    3 * ( ( IYPMY + 4900 ) / 100 ) ) / 4
     :                + ID - 2432076)
            DJM=DBLE(nint(DJM*1.d3))*1.d-3
*        Bad month
         ELSE
            J = -2
         END IF
      END IF

      END
!
!***************************************************************************
!
      SUBROUTINE JD2CAL ( DJ1, DJ2, IY, IM, ID, FD, J )
*  - - - - - - - - - - -
*   i a u _ J D 2 C A L
*  - - - - - - - - - - -
*  Julian Date to Gregorian year, month, day, and fraction of a day.
*
*  Given:
*     DJ1,DJ2     d     Julian Date (Notes 1, 2)
*
*  Returned:
*     IY          i     year
*     IM          i     month
*     ID          i     day
*     FD          d     fraction of day
*     J           i     status:
*                           0 = OK
*                          -1 = unacceptable date (Note 3)
*
*  Notes:
*
*  1) The earliest valid date is -68569.5 (-4900 March 1).  The
*     largest value accepted is 10^9.
*
*  2) The Julian Date is apportioned in any convenient way between
*     the arguments DJ1 and DJ2.  For example, JD=2450123.7 could
*     be expressed in any of these ways, among others:
*
*             DJ1            DJ2
*
*         2450123.7D0        0D0        (JD method)
*          2451545D0      -1421.3D0     (J2000 method)
*         2400000.5D0     50123.2D0     (MJD method)
*         2450123.5D0       0.2D0       (date & time method)
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      real*8 DJ1, DJ2
      INTEGER IY, IM, ID
      real*8 FD
      INTEGER J

*  Minimum and maximum allowed JD
      real*8 DJMIN, DJMAX
      PARAMETER ( DJMIN = -68569.5D0, DJMAX = 1D9 )

      INTEGER JD, L, N, I
      real*8 DJ, D1, D2, F1, F2, F, D

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Check if date is acceptable.
      DJ = DJ1 + DJ2
      IF ( DJ.LT.DJMIN .OR. DJ.GT.DJMAX ) THEN
         J = -1
      ELSE
         J = 0

*     Copy the date, big then small, and re-align to midnight.
         IF ( DJ1 .GE. DJ2 ) THEN
            D1 = DJ1
            D2 = DJ2
         ELSE
            D1 = DJ2
            D2 = DJ1
         END IF
         D2 = D2 - 0.5D0

*     Separate day and fraction.
         F1 = MOD(D1,1D0)
         F2 = MOD(D2,1D0)
         F = MOD(F1+F2,1D0)
         IF ( F .LT. 0D0 ) F = F+1D0
         D = ANINT(D1-F1) + ANINT(D2-F2) + ANINT(F1+F2-F)
         JD = NINT(D) + 1

*     Express day in Gregorian calendar.
         L = JD + 68569
         N = ( 4*L ) / 146097
         L = L - ( 146097*N + 3 ) / 4
         I = ( 4000 * (L+1) ) / 1461001
         L = L - ( 1461*I ) / 4 + 31
         J = ( 80*L ) / 2447
         ID = L - ( 2447*J ) / 80
         L = J / 11
         IM = J + 2 - 12*L
         IY = 100 * ( N-49 ) + I + L

         FD = F
         J = 0
      END IF
      END
!
!******************************************************
!
      real*8 FUNCTION EPJ ( DJ1, DJ2 )
*  - - - - - - - -
*   i a u _ E P J
*  - - - - - - - -
*  Julian Date to Julian Epoch.
*
*  Given:
*     DJ1,DJ2   d         Julian Date (see note)
*
*  The result is the Julian Epoch.
*
*  Note:
*
*     The Julian Date is supplied in two pieces, in the usual SOFA
*     manner, which is designed to preserve time resolution.  The
*     Julian Date is available as a single number by adding DJ1 and
*     DJ2.  The maximum resolution is achieved if DJ1 is 2451545D0
*     (J2000).
*-----------------------------------------------------------------------

      IMPLICIT NONE

      real*8 DJ1, DJ2

*  Reference epoch (J2000), JD
      real*8 DJ0
      PARAMETER ( DJ0 = 2451545D0 )

*  Days per Julian year
      real*8 DJY
      PARAMETER ( DJY = 365.25D0 )

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      EPJ = 2000D0 + ( ( DJ1-DJ0 ) + DJ2 ) / DJY
      END
!
!
!
      SUBROUTINE EPJ2JD ( EPJ, DJM0, DJM )
*  - - - - - - - - - - -
*   i a u _ E P J 2 J D
*  - - - - - - - - - - -
*  Julian Epoch to Julian Date.
*
*  Given:
*     EPJ         d     Julian Epoch (e.g. 1996.8D0)
*
*  Returned:
*     DJM0        d     MJD zero-point: always 2400000.5
*     DJM         d     Modified Julian Date
*-----------------------------------------------------------------------

      IMPLICIT NONE

      real*8 EPJ, DJM0, DJM

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      DJM0 = 2400000.5D0
      DJM  =   51544.5D0 + ( EPJ-2000D0 ) * 365.25D0

      END
!***************************************************************************
      double precision function gpsleap(tsec)

*** return total leap seconds since GPS epoch 1980jan06

*** note: does **NOT** return the full TAI-UTC delta
*** input time is GPS seconds -- initialized by setjd0()
*** Y2K -- only functional between 1980jan06-00:00:00  (GPS time start)
***                            and hard-coded date

      implicit double precision(a-h,o-z)

***** "Julian Date Converter"
***** http://aa.usno.navy.mil/data/docs/JulianDate.php
***** "Bulletin C"
***** http://hpiers.obspm.fr/eoppc/bul/bulc/bulletinc.dat
***** parameter(mjdhard=55196)            !*** cut-off date 2009dec31
***** parameter(mjdhard=55377)            !*** cut-off date 2010jun30
***** parameter(mjdhard=55561)            !*** cut-off date 2010dec31
***** parameter(mjdhard=55742)            !*** cut-off date 2011jun30
***** parameter(mjdhard=55926)            !*** cut-off date 2011dec31
***** parameter(mjdhard=56108)            !*** cut-off date 2012jun30
***** parameter(mjdhard=56292)            !*** cut-off date 2012dec31
***** parameter(mjdhard=56473)            !*** cut-off date 2013jun30
***** parameter(mjdhard=56657)            !*** cut-off date 2013dec31
***** parameter(mjdhard=56838)            !*** cut-off date 2014jun30
***** parameter(mjdhard=57022)            !*** cut-off date 2014dec31
***** parameter(mjdhard=57203)            !*** cut-off date 2015jun30
      parameter(mjdhard=57387)            !*** cut-off date 2015dec31

      save  /mjdoff/
      common/mjdoff/mjd0

*** clone for tests (and do any rollover)

      ttsec=tsec
      mjd0t=mjd0

    1 if(ttsec.ge.86400.d0) then
        ttsec=ttsec-86400.d0
        mjd0t=mjd0t+1
        go to 1
      endif

    2 if(ttsec.lt.0.d0) then
        ttsec=ttsec+86400.d0
        mjd0t=mjd0t-1
        go to 2
      endif

*** test date limits

      if(mjd0t.lt.44244) then             !*** 1980jan06
      !FATAL ERROR -cut-off date underflow in gpsleap()
        tai_utc = 0; stop 66767
      endif

*** http://maia.usno.navy.mil/ser7/tai-utc.dat
*** 1980 JAN  1 =JD 2444239.5  TAI-UTC=  19.0s
*** 1981 JUL  1 =JD 2444786.5  TAI-UTC=  20.0s
*** 1982 JUL  1 =JD 2445151.5  TAI-UTC=  21.0s
*** 1983 JUL  1 =JD 2445516.5  TAI-UTC=  22.0s
*** 1985 JUL  1 =JD 2446247.5  TAI-UTC=  23.0s
*** 1988 JAN  1 =JD 2447161.5  TAI-UTC=  24.0s
*** 1990 JAN  1 =JD 2447892.5  TAI-UTC=  25.0s
*** 1991 JAN  1 =JD 2448257.5  TAI-UTC=  26.0s
*** 1992 JUL  1 =JD 2448804.5  TAI-UTC=  27.0s
*** 1993 JUL  1 =JD 2449169.5  TAI-UTC=  28.0s
*** 1994 JUL  1 =JD 2449534.5  TAI-UTC=  29.0s
*** 1996 JAN  1 =JD 2450083.5  TAI-UTC=  30.0s
*** 1997 JUL  1 =JD 2450630.5  TAI-UTC=  31.0s
*** 1999 JAN  1 =JD 2451179.5  TAI-UTC=  32.0s
*** 2006 JAN  1 =JD 2453736.5  TAI-UTC=  33.0s
*** 2009 JAN  1 =JD 2454832.5  TAI-UTC=  34.0s
*** 2012 JUL  1 =JD 2456109.5  TAI-UTC=  35.0s
*** 2015 JUL  1 =JD 2457204.5  TAI-UTC=  36.0s

*** test against newest leaps first

      if    (mjd0t.ge.57204) then       !*** 2015 JUL 1 = 57204
        tai_utc = 36.d0
      elseif(mjd0t.ge.56109) then       !*** 2012 JUL 1 = 56109
        tai_utc = 35.d0
      elseif(mjd0t.ge.54832) then       !*** 2009 JAN 1 = 54832
        tai_utc = 34.d0
      elseif(mjd0t.ge.53736) then       !*** 2006 JAN 1 = 53736
        tai_utc = 33.d0
      elseif(mjd0t.ge.51179) then       !*** 1999 JAN 1 = 51179
        tai_utc = 32.d0
      elseif(mjd0t.ge.50630) then       !*** 1997 JUL 1 = 50630
        tai_utc = 31.d0
      elseif(mjd0t.ge.50083) then       !*** 1996 JAN 1 = 50083
        tai_utc = 30.d0
      elseif(mjd0t.ge.49534) then       !*** 1994 JUL 1 = 49534
        tai_utc = 29.d0
      elseif(mjd0t.ge.49169) then       !*** 1993 JUL 1 = 49169
        tai_utc = 28.d0
      elseif(mjd0t.ge.48804) then       !*** 1992 JUL 1 = 48804
        tai_utc = 27.d0
      elseif(mjd0t.ge.48257) then       !*** 1991 JAN 1 = 48257
        tai_utc = 26.d0
      elseif(mjd0t.ge.47892) then       !*** 1990 JAN 1 = 47892
        tai_utc = 25.d0
      elseif(mjd0t.ge.47161) then       !*** 1988 JAN 1 = 47161
        tai_utc = 24.d0
      elseif(mjd0t.ge.46247) then       !*** 1985 JUL 1 = 46247
        tai_utc = 23.d0
      elseif(mjd0t.ge.45516) then       !*** 1983 JUL 1 = 45516
        tai_utc = 22.d0
      elseif(mjd0t.ge.45151) then       !*** 1982 JUL 1 = 45151
        tai_utc = 21.d0
      elseif(mjd0t.ge.44786) then       !*** 1981 JUL 1 = 44786
        tai_utc = 20.d0
      elseif(mjd0t.ge.44239) then       !*** 1980 JAN 1 = 44239
        tai_utc = 19.d0

*** should never get here

      else
      !FATAL ERROR --fell thru tests in gpsleap()
        tai_utc = 0; stop 66768
      endif

*** convert TAI-UTC into GPS leap seconds

      gpsleap = tai_utc - 19.d0

      return
      end
!
!************************************************************************
!
      subroutine tmcnt(tm,iyr,imo,idy,ihr,imn,sec)
      implicit none
	integer::iyr,imo,idy,ihr,imn,kln
	real*8::tm,sec,tmp,dj0,fd
	character*40::tmstr,astr,aymd,ahms
!-------------------------------------------------------------------------
      write(astr,*)tm
      astr=trim(adjustl(astr));kln=len(astr)
      read(astr(1:4),*)iyr
      read(astr(5:6),*)imo
      read(astr(7:8),*)idy
      ihr=0;imn=0;sec=0.d0
      if(kln>9)read(astr(9:10),*)ihr!write(ihr,*)astr(9:10)
      if(kln>11)read(astr(11:12),*)imn!write(imn,*)astr(11:12)
      if(kln>13)read(astr(13:kln),*)sec!write(sec,*)astr(13:kln)
      continue
      end
!
!************************************************************************
!
      subroutine tmtostr(tm,tmstr)
      implicit none
	character*800::tmstr,astr
	real*8::tm
!-------------------------------------------------------------------------
      write(astr,'(F16.0)')tm
      astr=trim(adjustl(astr))
      write(tmstr,*)astr(1:len_trim(astr)-1)
      return
      end
!
!************************************************************************
!
      subroutine mjdtotm(mjd0,ltm)
      implicit none
	integer::IY,IM,ID,ihr,imn,k
	real*8::ltm,mjd0,mjd,sec,tmp
      mjd=mjd0
      call JD2CAL(2400000.5d0,mjd,IY,IM,ID,tmp,k)
      mjd=dble(nint(mjd*1.d6))*1.d-6;tmp=tmp+1.d-8
      ihr=floor(tmp*24.d0);tmp=tmp*24.d0-dble(ihr)
      imn=floor(tmp*6.d1);tmp=tmp*6.d1-dble(imn)
      sec=tmp*6.d1;
      if(nint(sec)-59.5>0)then
        sec=0.d0;imn=imn+1
      endif
      if(imn-59.5>0)then
        imn=0;ihr=ihr+1
      endif
      if(sec>5.d-2)ltm=IY*1.d10+IM*1.d8+ID*1.d6+ihr*1.d4+imn*1.d2+sec
      if(sec<5.d-2.and.imn>0.1)ltm=IY*1.d8+IM*1.d6+ID*1.d4+ihr*1.d2+imn
      if(sec<5.d-2.and.imn<0.1.and.ihr>0.1)ltm=IY*1.d6+IM*1.d4+ID*1.d2+ihr
      if(sec<5.d-2.and.imn<0.1.and.ihr<0.1)ltm=IY*1.d4+IM*1.d2+ID
      if(ltm<1901010131.d0)ltm=ltm*1.d2
      end
!
!************************************************************************
!
      subroutine mjdmnth(mjd1,mjd2,mjdn,nm)
      implicit none
	integer::IY,IM,ID,i,j,k,nm
	real*8::JD0,FD,mjd1,mjd2,mjdn(80000),mjd,mjd0
      JD0=2400000.5d0;ID=1
      call JD2CAL(JD0,mjd1,IY,IM,ID,FD,j)
      call CAL2JD (IY,IM,1,mjd0,j)
      if(IM==12)then
        IY=IY+1;IM=1
      else
        IY=IY;IM=IM+1
      endif
      call CAL2JD (IY,IM,ID,mjd,j)
      mjdn(1)=mjd;nm=1
      do while(mjd<mjd2)
        if(IM==12)then
          IY=IY+1;IM=1
        else
          IY=IY;IM=IM+1
        endif
        call CAL2JD (IY,IM,ID,mjd,j)
        nm=nm+1;mjdn(nm)=mjd
      enddo
      mjdn(nm+1)=mjd0
      end
!
!************************************************************************
!
      subroutine ltmtomathgls(stm)
!MS$ATTRIBUTES DLLEXPORT::LTMTOMATHGLS
      implicit none
	real*8::stm(2),ltm,mjd,mjd0,sec
	integer::IY,IM,ID,ih,imn,k
      ltm=stm(1)
      call tmcnt(ltm,IY,IM,ID,ih,imn,sec)
      call CAL2JD (IY,IM,ID,mjd,k)
      mjd=mjd+dble(ih)/24.d0+dble(imn)/1440.d0+sec/864.d2 !GPS_MJD
      call CAL2JD (1970,1,1,mjd0,k)
      stm(2)=(mjd-mjd0)*864.d2
      end
