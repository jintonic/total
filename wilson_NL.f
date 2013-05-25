c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      Subroutine wilson_nl(time_input, e_nu_input, 
     &     dNL_nue_out, dNL_anue_out, dNL_mutau_out)
      Implicit None
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     This gives the differential number luminosity [/MeV/sec] as a 
c     function of time after bounce [sec] and neutrino energy [MeV]
c     
c     Note: 1.1 msec < t < 18 sec
c-----------------------------------------------------------------------
c inputs:
      real*8 time_input, e_nu_input ! [sec], [MeV]
c output:
      real*8 dNL_nue_out, dNL_anue_out, dNL_mutau_out ! [/sec/MeV]

c constants:
      character*16 filename1, filename2
      parameter (filename1 = 
     &  'wilson-early.dat')
      parameter (filename2 = 
     &  'wilson-late.dat')
      integer Lunin
      parameter (Lunin = 7)
      integer time_num, spec_num
      parameter (time_num = 63)
      parameter (spec_num = 20)
      real*8 MeVerg
      parameter (MeVerg = 1.6D-6) ! MeV --> erg
      real*8 t_start
      parameter (t_start = 2d-3) !  Defining time_input = 0 as 2d-3

c local:
      integer i, j, j_dum
      real*8 time(0:time_num)  ! [sec] after the bounce
      save time
      real*8 E_nue(0:time_num), E_anue(0:time_num), E_mutau(0:time_num),
     &       E_tot(0:time_num) ! Energy emitted up to the present time [erg]
      real*8 N_nue(0:time_num), N_anue(0:time_num), N_mutau(0:time_num),
     &       N_tot(0:time_num) ! Number emitted up to the present time

      real*8 dum1, dum2, dum3, dum4, dum5, dum6
      real*8 neusph_nue(0:time_num), neusph_mutau(0:time_num) ! [cm]
      real*8 ave_nue(0:time_num), ave_anue(0:time_num), 
     &       ave_mutau(0:time_num) ! [MeV]

      real*8 group_e(0:spec_num) ! [MeV]
      real*8 log_group_e(1:spec_num)
      save group_e, log_group_e
      
      real*8 spec_nue(0:time_num, 0:spec_num-1), 
     &     spec_anue(0:time_num, 0:spec_num-1), 
     &     spec_mutau(0:time_num, 0:spec_num-1)

      real*8 egroup_nue(0:time_num, 0:spec_num), 
     &       egroup_anue(0:time_num, 0:spec_num), 
     &       egroup_mutau(0:time_num, 0:spec_num)

      real*8 E_integ_nue, E_integ_anue, E_integ_mutau 

      real*8 NL_nue(time_num-1), NL_anue(time_num-1), ! Number Luminosity
     &       NL_mutau(time_num-1), NL_tot(time_num-1) !           [/sec]
c      real*8 NL1_nue, NL2_nue
c      real*8 NL1_anue, NL2_anue
c      real*8 NL1_mutau, NL2_mutau
c      real*8 NL1_tot, NL2_tot

      real*8 dNLde_nue(time_num-1, 0:spec_num-1),  ! Number Luminosity Spectrum
     &       dNLde_anue(time_num-1, 0:spec_num-1), !                [/sec/MeV]
     &       dNLde_mutau(time_num-1, 0:spec_num-1)
      save dNLde_nue, dNLde_anue, dNLde_mutau

      real*8 log_dNL(1:spec_num) ! for cubic spline 
      real*8 dNL(1:spec_num) ! for cubic spline 
!!!!!!! range 1-spec_num for using nrecipe routine

      logical first  /.true./
      save first

c begin:
      character*1 slash
      parameter (slash = '/')
      character*128 DATA_DIR
      call get_environment_variable("TOTAL_DATA_DIR", DATA_DIR)

      time_input = time_input + t_start
c     (time_input = 0d0 ----> time_input = initial time of data)

      Do i = 1, time_num
         spec_nue(i,0) = 0d0
         spec_anue(i,0) = 0d0
         spec_mutau(i,0) = 0d0
         egroup_nue(i, 0) = 0d0
         egroup_anue(i, 0) = 0d0
         egroup_mutau(i, 0) = 0d0
      End do

      if (first) then 
         Goto 1 
      else
         Goto 1500
      end if

c----------------------------------------------------------------------
c Data Reading:
c----------------------------------------------------------------------
 1    group_e(0) = 0d0

      i = 1
      Open (Lunin, file=TRIM(DATA_DIR)//slash//filename1)
      Do while (.true.)
 10      read (Lunin, *, ERR = 10, END = 555) time(i)
 20      read (Lunin, *, ERR = 20, END = 555) 
     &          E_nue(i), E_anue(i), E_mutau(i), E_tot(i)
         E_mutau(i) = E_mutau(i) / 4d0
 30      read (Lunin, *, ERR = 30, END = 555) 
     &          N_nue(i), N_anue(i), N_mutau(i), N_tot(i)
         N_mutau(i) = N_mutau(i) / 4d0
 40      read (Lunin, *, ERR = 40, END = 555) 
     &          dum1, dum2, dum3
 50      read (Lunin, *, ERR = 50, END = 555) 
     &          neusph_nue(i), neusph_mutau(i), dum1
 60      read (Lunin, *, ERR = 60, END = 555) 
     &          dum1, dum2
 70      read (Lunin, *, ERR = 70, END = 555) 
     &          ave_nue(i), ave_anue(i), ave_mutau(i)
         ave_nue(i) = ave_nue(i) / 1d3
         ave_anue(i) = ave_anue(i) / 1d3 
         ave_mutau(i) = ave_mutau(i) / 1d3

         Do j = 1, spec_num
 80         read (Lunin, *, ERR = 80, END = 555) j_dum, group_e(j),
     &       dum1, dum2, dum3, dum4, dum5, dum6,
     &       egroup_nue(i,j), egroup_anue(i,j), egroup_mutau(i,j)
            group_e(j) = group_e(j) / 1d3 ! keV -> MeV
         End do

         i = i + 1
      End do
 555  close (Lunin)


      Open (Lunin, file=TRIM(DATA_DIR)//slash//filename2)
      Do while (.true.)
 610      read (Lunin, *, ERR = 610, END = 999) time(i)
 620      read (Lunin, *, ERR = 620, END = 999) 
     &          E_nue(i), E_anue(i), E_mutau(i), E_tot(i)
          E_mutau(i) = E_mutau(i) / 4d0
 630      read (Lunin, *, ERR = 630, END = 999) 
     &          N_nue(i), N_anue(i), N_mutau(i), N_tot(i)
          N_mutau(i) = N_mutau(i) / 4d0
 640      read (Lunin, *, ERR = 640, END = 999) 
     &          dum1, dum2, dum3
 650      read (Lunin, *, ERR = 650, END = 999) 
     &          neusph_nue(i), neusph_mutau(i), dum1
 660      read (Lunin, *, ERR = 660, END = 999) 
     &          dum1, dum2
 670      read (Lunin, *, ERR = 670, END = 999) 
     &          ave_nue(i), ave_anue(i), ave_mutau(i)
         ave_nue(i) = ave_nue(i) / 1d3
         ave_anue(i) = ave_anue(i) / 1d3 
         ave_mutau(i) = ave_mutau(i) / 1d3

         Do j = 1, spec_num
 680         read (Lunin, *, ERR = 680, END = 999) j_dum, group_e(j),
     &       egroup_nue(i,j), egroup_anue(i,j), egroup_mutau(i,j)
            group_e(j) = group_e(j) / 1d3 ! keV -> MeV
         End do

         read (Lunin, *)   ! reading dummy line `21'
c         write (*,*) i, time(i)
         i = i + 1
      End do
 999  close (Lunin)

c-----------------------------------------------------------------------
c Spectrum Generating:
c group distribution -> per unit energy
      Do i = 1, time_num 
         Do j = 1, spec_num - 1
            spec_nue(i, j) = egroup_nue(i,j) 
     &           / (group_e(j+1)-group_e(j-1))
            spec_anue(i, j) = egroup_anue(i,j) 
     &           / (group_e(j+1)-group_e(j-1))
            spec_mutau(i, j) = egroup_mutau(i,j) 
     &           / (group_e(j+1)-group_e(j-1))
         End do

c normalization:            
         E_integ_nue = 0d0
         E_integ_anue = 0d0
         E_integ_mutau = 0d0
         Do j = 1, spec_num - 1
            E_integ_nue = E_integ_nue 
     &           + (spec_nue(i,j-1) + spec_nue(i,j))
     &           * (group_e(j) - group_e(j-1)) * 0.5d0
            E_integ_anue = E_integ_anue
     &                  + (spec_anue(i,j-1) + spec_anue(i,j))
     &                  * (group_e(j) - group_e(j-1)) * 0.5d0
            E_integ_mutau = E_integ_mutau 
     &                  + (spec_mutau(i,j-1) + spec_mutau(i,j))
     &                  * (group_e(j) - group_e(j-1)) * 0.5d0
         End do
         Do j = 0, spec_num - 1
            spec_nue(i,j) = spec_nue(i,j) / E_integ_nue
            spec_anue(i,j) = spec_anue(i,j) / E_integ_anue
            spec_mutau(i,j) = spec_mutau(i,j) / E_integ_mutau
         End do
      End do

c Number Luminosity Calculation
c      1bin
      Do i = 2, time_num - 1
         NL_nue(i) = (N_nue(i) - N_nue(i-1)) 
     &        / (time(i) - time(i-1))
         NL_anue(i) = (N_anue(i) - N_anue(i-1)) 
     &        / (time(i) - time(i-1))
         NL_mutau(i) = (N_mutau(i) - N_mutau(i-1)) 
     &        / (time(i) - time(i-1))
         NL_tot(i) = (N_tot(i) - N_tot(i-1)) 
     &        / (time(i) - time(i-1))
c         write (*,*) i, time(i), NL_nue(i)
      End do
      

c Differential Number Luminosity Calculation:
      Do i = 2, time_num-1
         Do j = 0, spec_num - 1
            dNLde_nue(i,j) = NL_nue(i) * spec_nue(i,j)
            dNLde_anue(i,j) = NL_anue(i) * spec_anue(i,j)
            dNLde_mutau(i,j) = NL_mutau(i) * spec_mutau(i,j)
         End do
      End do


c Log-Data cal. for Natural Cubic Spline:
      Do j = 1, spec_num - 1
         log_group_e(j) = dlog10( group_e(j) )
      End do

      first = .false.


c----------------------------------------------------------------------
c     Output Calculation:
c----------------------------------------------------------------------
 1500 i = 2
      Do while (time_input .gt. time(i))
         i = i + 1
      End do
      j = 0
      Do while (e_nu_input .gt. group_e(j))
         j = j + 1
      End do

c for nue:
      Do j = 1, spec_num - 1
         Call linear_int(time(i-1), time(i), 
     &        dNLde_nue(i-1, j), dNLde_nue(i, j), 
     &        time_input, dNL(j))
         log_dNL(j) = dlog10(dNL(j))
      End do
      call nat_cub_sp(log_group_e, log_dNL, spec_num-1, 
     &     dlog10(e_nu_input), dNL_nue_out)
      dNL_nue_out = 10d0**dNL_nue_out

c for anue:
      Do j = 1, spec_num - 1
         Call linear_int(time(i-1), time(i), 
     &        dNLde_anue(i-1, j), dNLde_anue(i, j), 
     &        time_input, dNL(j))
         log_dNL(j) = dlog10(dNL(j))
c      write (*,*) j, log_dNL(j), 10d0**log_group_e(j)
      End do
      call nat_cub_sp(log_group_e, log_dNL, spec_num-1, 
     &     dlog10(e_nu_input), dNL_anue_out)
c      write (*,*) dNL_anue_out
      dNL_anue_out = 10d0**dNL_anue_out

c for mutau:
      Do j = 1, spec_num - 1
         Call linear_int(time(i-1), time(i), 
     &        dNLde_mutau(i-1, j), dNLde_mutau(i, j), 
     &        time_input, dNL(j))
         log_dNL(j) = dlog10(dNL(j))
      End do
      call nat_cub_sp(log_group_e, log_dNL, spec_num-1, 
     &     dlog10(e_nu_input), dNL_mutau_out)
      dNL_mutau_out = 10d0**dNL_mutau_out


c------------------------------------------------------------------------

      time_input = time_input - t_start

c--------------------------- Neutronization Burst --------------------
      If ((time_input .ge. 4d-2) .and. (time_input .le. 5d-2)) then
         call wilson_nb_NL(time_input, e_nu_input, 
     &     dNL_nue_out, dum1, dum2)

c         write (*,*) time_input, dNL_nue_out

         dNL_nue_out = dNL_nue_out * (1d0 * (5d-2 - time_input)/1d-2
     &        + 8.59d0/13.82d0 * (time_input - 4d-2) / 1d-2)

c         write (*,*) time_input, dNL_nue_out
      End if      

c---------------------------------------------------------------------

      return
      end
