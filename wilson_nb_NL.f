c      program test
c      real*8 dNL_nue, dNL_anue, dNL_mutau
c      Call wilson_nb_NL(0.045d0, 10d0, dNL_nue, dNL_anue, dNL_mutau)
c      write (*,*) dNL_nue, dNL_anue, dNL_mutau
c      end

c      include '../../../tools/interporlate/linear_int.f'
c      include '../../../tools/interporlate/nat_cub_sp.f'
c      include '../../../tools/nrecipe/spline.f'
c      include '../../../tools/nrecipe/splint.f'

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      Subroutine wilson_nb_NL(time_input, e_nu_input, 
     &     dNL_nue_out, dNL_anue_out, dNL_mutau_out)
      Implicit None
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     This gives the differential number luminosity [/MeV/sec] as a 
c     function of time after bounce [sec] and neutrino energy [MeV]
c     during NEUTRONIZATION BURST
c     Note: 40 msec < t < 50 msec
c-----------------------------------------------------------------------
c global:
c      include '../constants/neutrinoburst_const.f'
c inputs:
      real*8 time_input, e_nu_input ! [sec], [MeV]
c output:
      real*8 dNL_nue_out, dNL_anue_out, dNL_mutau_out ! [/sec/MeV]

c constants:
      character*60 filename
      parameter (filename = 
     &  'wilson-nb.dat')
      integer Lunin
      parameter (Lunin = 7)
      integer time_num, spec_num
      parameter (time_num = 100)
      parameter (spec_num = 21)
      real*8 MeVerg
      parameter (MeVerg = 1.6D-6) ! MeV --> erg
      real*8 t_start
      parameter (t_start = 2d-3) !  Defining time_input = 0 as 2d-3
      real*8 time_nb_start
      parameter (time_nb_start = 5.075d-1) ! sec

      real*8 erg_MeV
      parameter (erg_MeV = 1.602d-6)  ! erg/MeV

c local:
      integer i, j, j_dum
      real*8 time(0:time_num)  ! [sec] 
      real*8 tab(0:time_num)   ! [sec] after the bounce
      save time
      real*8 time_loc

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

      real*8 dNLde_nue(time_num-1, 0:spec_num-1),  ! Number Luminosity Spectrum
     &       dNLde_anue(time_num-1, 0:spec_num-1), !                [/sec/MeV]
     &       dNLde_mutau(time_num-1, 0:spec_num-1)
      save dNLde_nue, dNLde_anue, dNLde_mutau

      real*8 L_nue(1:time_num), L_anue(1:time_num), L_mutau(1:time_num)
      
      character*6 cdum1
      character*7 cdum2
      character*8 cdum3
      character*17 cdum4
      
      real*8 log_dNL(1:spec_num) ! for cubic spline 
      real*8 dNL(1:spec_num) ! for cubic spline 
!!!!!!! range 1-spec_num for using nrecipe routine

      logical first  /.true./
      save first

c begin:
c      write (*,*) 
c     &     'No.  Time, TaB [sec], L_nue, anue, mutau [10^50erg/s]'

      time_loc = time_input + time_nb_start - 0.04d0
c                                             ^^^^^^ 0.04--0.05 in origin

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
      Open (Lunin, file=filename)
      Do while (.true.)
 10      read (Lunin, 9001, ERR = 10, END = 555) 
     &        cdum1, time(i), cdum4, tab(i)
 20      read (Lunin, 9002, ERR = 20, END = 555) cdum1, L_nue(i),
     &        cdum2, L_anue(i), cdum3, L_mutau(i)
c         write (*, 5000) i, time(i), tab(i), 
c     &        L_nue(i)/1d50, L_anue(i)/1d50, L_mutau(i)/1d50

         Do j = 1, spec_num
 80         read (Lunin, *, ERR = 80, END = 555) j_dum, group_e(j),
     &       egroup_nue(i,j), egroup_anue(i,j), egroup_mutau(i,j)
            group_e(j) = group_e(j) / 1d3 ! keV -> MeV
c            write (*,*) j, group_e(j)
         End do

         Do j = 1, 2
            read (Lunin, *)
         End do

         i = i + 1
      End do
 555  close (Lunin)


c-------------------------------------------------------------------------
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
         E_integ_nue = 0d0    ! luminosity integration (not Ln)
         E_integ_anue = 0d0
         E_integ_mutau = 0d0
         Do j = 1, spec_num - 1
            E_integ_nue = E_integ_nue 
     & + (spec_nue(i,j-1)*group_e(j-1) + spec_nue(i,j)*group_e(j))
     &           * (group_e(j) - group_e(j-1)) * 0.5d0
            E_integ_anue = E_integ_anue
     & + (spec_anue(i,j-1)*group_e(j-1) + spec_anue(i,j)*group_e(j))
     &                  * (group_e(j) - group_e(j-1)) * 0.5d0
            E_integ_mutau = E_integ_mutau 
     & + (spec_mutau(i,j-1)*group_e(j-1) + spec_mutau(i,j)*group_e(j))
     &                  * (group_e(j) - group_e(j-1)) * 0.5d0
         End do
         Do j = 0, spec_num - 1
            spec_nue(i,j) = spec_nue(i,j) 
     &           / E_integ_nue * L_nue(i) / erg_MeV
            spec_anue(i,j) = spec_anue(i,j) 
     &           / E_integ_anue * L_anue(i)/ erg_MeV
            spec_mutau(i,j) = spec_mutau(i,j) 
     &           / E_integ_mutau * L_mutau(i)/ erg_MeV
         End do

      End do

c Differential Number Luminosity Calculation:
      Do i = 2, time_num-1
         Do j = 0, spec_num - 1
            dNLde_nue(i,j) = spec_nue(i,j)
            dNLde_anue(i,j) = spec_anue(i,j)
            dNLde_mutau(i,j) = spec_mutau(i,j)
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
      Do while (time_loc .gt. time(i))
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
     &        time_loc, dNL(j))
         log_dNL(j) = dlog10(dNL(j))
      End do
      call nat_cub_sp(log_group_e, log_dNL, spec_num-1, 
     &     dlog10(e_nu_input), dNL_nue_out)
      dNL_nue_out = 10d0**dNL_nue_out

c for anue:
      Do j = 1, spec_num - 1
         Call linear_int(time(i-1), time(i), 
     &        dNLde_anue(i-1, j), dNLde_anue(i, j), 
     &        time_loc, dNL(j))
         log_dNL(j) = dlog10(dNL(j))
      End do
      call nat_cub_sp(log_group_e, log_dNL, spec_num-1, 
     &     dlog10(e_nu_input), dNL_anue_out)
c      write (*,*) dNL_anue_out
      dNL_anue_out = 10d0**dNL_anue_out

c for mutau:
      Do j = 1, spec_num - 1
         Call linear_int(time(i-1), time(i), 
     &        dNLde_mutau(i-1, j), dNLde_mutau(i, j), 
     &        time_loc, dNL(j))
         log_dNL(j) = dlog10(dNL(j))
      End do
      call nat_cub_sp(log_group_e, log_dNL, spec_num-1, 
     &     dlog10(e_nu_input), dNL_mutau_out)
      dNL_mutau_out = 10d0**dNL_mutau_out


c------------------------------------------------------------------------
      return

 5000 format (1I4, 5D12.4)
 9001 format(1A6,1D12.4,2x,1A18,1D12.4)
 9002 format(5x,1A6,1D12.4,2x,1A7,1D12.4,2x,1A8,
     &        1D12.4)
      end
