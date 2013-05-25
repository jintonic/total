c----------------------------------------------------------------------------
c----------------------------------------------------------------------------
      Subroutine nat_cub_sp(xa, ya, n, x, y)
      Implicit None
c----------------------------------------------------------------------------
c     Natural Cubic Spline
c----------------------------------------------------------------------------
c input:
      integer n
      real*8 xa(1:n), ya(1:n), x
c output:
      real*8 y

c local:
      integer n_max
      parameter (n_max = 10000)
      real*4 y2a_4(1:n_max) ! second derivative
      real*4 y4
      real*4 xa4(1:n_max), ya4(1:n_max)

      integer i

c begin:
      If (n .gt. n_max) then
         pause 'n_max > n  in nat_cub_sp.f'
      End if

      Do i = 1, n
         xa4(i) = real(xa(i))
         ya4(i) = real(ya(i))
      End do

      call spline(xa4, ya4, n, 1e30, 1e30, y2a_4)

      call splint(xa4, ya4, y2a_4, n, real(x), y4)
      y = dble(y4)

      return
      
      end
