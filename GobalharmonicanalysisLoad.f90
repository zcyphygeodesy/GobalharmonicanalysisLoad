!  GobalharmonicanalysisLoad.f90 
!
!  FUNCTIONS:
!  GobalharmonicanalysisLoad - Entry point of console application.
!
!****************************************************************************
      program GobalharmonicanalysisLoad
      implicit none
	character*800::dtmfl,loadfl
	integer knd
	real*8::dk,itd
!---------------------------------------------------------------------
      !����ر�������knd
      !knd=0 landwater EWH variation, =1 sea level variation =-1 surface atmosphere variation
      knd=1!knd=0½��ˮ�仯��=1��ƽ��仯��=-1�������ѹ�仯
      !����½����������������ļ��������ڷ���½������knd=0 or =1
      !The spatial resolution of the land-sea terrain grid should not be lower than that of the surface load grid.
      write(dtmfl,*)'sphetopo30m.dat'
      !���븺�����������ʱ�������ļ����ļ�,���������ͷ�ļ����߸���Ϊ������ETideLoadϵͳ��ʽʱ��
      !The seventh numerical value of spherical coordinate load ewh grid file header is the long integer time agreed by ETideLoad
      write(loadfl,*)'sphgridtmsquefls.txt'
      !���������ֹ�������в��׼��Ϊԭ��ֵ��׼��dk��
      !input the iteration condition parameters
      !Iteration termination condition: The standard deviation of the residual grid value is less than dk of the standard deviation
      !of the original grid value, or the difference of the residual standard deviation of the previous step iteration relative to 
      !the current step iteration is less than itd of the standard deviation of the original grid values.
      dk=1.d-3; itd=1.d-4
      write(*, *)"    Begin compulation......"
      call Loadharmanalysis(dtmfl,loadfl,knd,dk,itd)
      pause
      end
