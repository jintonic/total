      program test
      implicit none

      real*8 time               ! sec
      real*8 nu_energy          ! MeV
      real*8 dNL_nue, dNL_anue, dNL_mutau ! /sec/MeV
c       differential number luminosity for nu_e, nu_e_bar, nu_mutau
c        (nu_mutau = one of nu_mu, nu_tau, and their antiparticles)

c begin:
      time = 4.165d-2
      nu_energy = 10d0

      Call wilson_NL(time, nu_energy, dNL_nue, dNL_anue, dNL_mutau)
      write (*,*) dNL_nue, dNL_anue, dNL_mutau

      end

c subroutines required:
      include 'wilson_NL.f'
      include 'wilson_nb_NL.f'
      include 'linear_int.f'
      include 'nat_cub_sp.f'
      include 'spline.f'
      include 'splint.f'

