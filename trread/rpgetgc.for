      subroutine rpgetgc(zchar)
C
C  get RPLOT calculator "guard" character for commands
C  i.e. leading character which means to parser, "this is a command"
C  (instead of an expression).
C
      use rpcalc_mod

      character*1 zchar
C
C----------------------
C
      zchar=gdchar
C
      return
      end
