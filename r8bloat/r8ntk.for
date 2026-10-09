      subroutine r8ntk(intk)
c
c  return the R,Z sequence grid size
c
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      include 'r8tkbparm.inc'
c
      intk=ntk
      return
      end
