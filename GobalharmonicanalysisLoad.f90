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
      !输入地表负荷类型knd
      !knd=0 landwater EWH variation, =1 sea level variation =-1 surface atmosphere variation
      knd=1!knd=0陆地水变化，=1海平面变化，=-1地面大气压变化
      !输入陆海地形球坐标格网文件名，用于分离陆海区域knd=0 or =1
      !The spatial resolution of the land-sea terrain grid should not be lower than that of the surface load grid.
      write(dtmfl,*)'sphetopo30m.dat'
      !输入负荷球坐标格网时间序列文件名文件,球坐标格网头文件第七个数为长整型ETideLoad系统格式时间
      !The seventh numerical value of spherical coordinate load ewh grid file header is the long integer time agreed by ETideLoad
      write(loadfl,*)'sphgridtmsquefls.txt'
      !输入迭代终止条件：残差标准差为原格值标准的dk倍
      !input the iteration condition parameters
      !Iteration termination condition: The standard deviation of the residual grid value is less than dk of the standard deviation
      !of the original grid value, or the difference of the residual standard deviation of the previous step iteration relative to 
      !the current step iteration is less than itd of the standard deviation of the original grid values.
      dk=1.d-3; itd=1.d-4
      write(*, *)"    Begin compulation......"
      call Loadharmanalysis(dtmfl,loadfl,knd,dk,itd)
      pause
      end
