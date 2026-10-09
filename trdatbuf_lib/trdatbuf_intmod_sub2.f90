
SUBROUTINE trdatbuf_intmod_sub2(s,ivarsize,idimvs,intbuf,idimi)
  ! ** gensis_f90.py generated code, do NOT edit **

  use trdatbuf_obj
  IMPLICIT NONE

  type(trdatbuf) :: s
  integer :: idimvs,idimi
  integer :: ivarsize(idimvs)
  integer :: intbuf(idimi)

  integer :: i,i1,i2

  !----------------------

  i=0
  i1=0
  i2=0

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LADR_ZA

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LADR_ZAF

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LAEPSP

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LAIMPS

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LAIMPX

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATALP

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATBDI

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATCUR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATDFL

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATDTG

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATDTS

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATDTX

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATEDI

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATEHP

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATELI

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATFIE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATFMN

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATFMX

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATGAS

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATGF3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATGF4

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATGF6

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATGFD

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATGFH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATGFT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATGIE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATL2B

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATLAD

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATLID

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATNTX

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATOGF

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATORC

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATPF0

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATPFH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATPFL

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATPLF

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATPOS

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATRBZ

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATRC3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATRC4

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATRC6

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATRCD

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATRCH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATRCT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATRCY

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATRMN

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATRTP

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATSBT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATSF3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATSFA

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATSFL

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATSFP

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATSFT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATSP3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATTET

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATTGF

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATTIT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATTPI

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATTQT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATTRC

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATTRF

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATTXI

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATVPH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATVSB

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATVSF

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATXKF

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATZEF

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATZFA

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATZIM

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDATZPL

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDMDF

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDMFR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDMMX

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDMSE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDPFC

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDPFC_PRE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDPFN

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDRBDY

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LDZBDY

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LEFD0

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LEFDB

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LEFDP

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LEFDQ

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LEFDR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LEFDS

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFBOL

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFBPA

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFBPB

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFD2F

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFDE2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFDF3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFDF4

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFDF6

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFDFD

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFDFH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFDFT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFECA

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFECB

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFECF

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFECQ

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFECX

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFFD0

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFFDB

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFFDP

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFFDQ

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFFDR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFFDS

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFGRB

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFLF3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFLF4

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFLF6

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFLFD

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFLFH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFLFT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFMSE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFNER

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFNI3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFNI4

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFNI6

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFNID

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFNIH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFNIM

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFNIT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFNMR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFOMG

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFPRS

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFPSI

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFQPR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFREE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFRQRFF

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSBI

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSIM

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYBOL

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYBPA

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYBPB

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYD2F

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYDE2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYDF3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYDF4

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYDF6

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYDFD

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYDFH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYDFT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYECF

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYGRB

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYLF3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYLF4

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYLF6

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYLFD

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYLFH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYLFT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYMSE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYNER

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYNI3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYNI4

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYNI6

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYNID

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYNIH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYNIM

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYNIT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYNMR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYOMG

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYPRS

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYQPR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYSBI

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYSIM

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYTER

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYTI2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYTI3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYTQI

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYV2F

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYVB2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYVC3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYVC4

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYVC6

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYVCD

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYVCH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYVCT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYVEE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYVIE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYVMO

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYVP2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYVPO

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYVPR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYVTR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFSYZF2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFTER

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFTI2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFTI3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFTQI

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFULNB

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFV2F

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFVB2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFVC3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFVC4

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFVC6

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFVCD

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFVCH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFVCT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFVEE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFVIE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFVMO

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFVP2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFVPO

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFVPR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFVTR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LFZF2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LHLFNB

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LLIM

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LMOMD

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LPELDA

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LPHSLH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LPMDF

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LPMSE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LPWREC

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LPWRLH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LPWRNB

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LPWRRF

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LRBQSP

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LRFS

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LRMFR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LRMP_BDY1

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LRMP_BDY2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LRMSE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LRP2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LRPL

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LRPSI

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSEVENT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYBOL

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYBPA

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYBPB

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYD2F

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYDE2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYDF3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYDF4

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYDF6

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYDFD

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYDFH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYDFT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYECF

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYGRB

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYLF3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYLF4

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYLF6

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYLFD

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYLFH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYLFT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYMSE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYNER

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYNI3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYNI4

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYNI6

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYNID

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYNIH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYNIM

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYNIT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYNMR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYOMG

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYPRS

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYQPR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYSBI

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYSIM

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYTER

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYTI2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYTI3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYTQI

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYV2F

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYVB2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYVC3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYVC4

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYVC6

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYVCD

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYVCH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYVCT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYVEE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYVIE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYVMO

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYVP2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYVPO

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYVPR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYVTR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LSSYZF2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LTABORT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LTFD0

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LTFDB

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LTFDP

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LTFDQ

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LTFDR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LTFDS

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LTHFS

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LTIME1

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LTIME2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LTIMEC

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LTIMECA

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LTIMECB

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LTIMLH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LTIMNB

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LTIMRF

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LTIMRFF

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LTPSI

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LTSAW

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LVLTNB

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXBOL

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXBPA

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXBPB

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXD2F

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXDE2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXDF3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXDF4

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXDF6

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXDFD

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXDFH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXDFT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXECF

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXFD0

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXFDB

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXFDP

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXFDQ

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXFDR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXFDS

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXFS

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXGRB

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXLF3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXLF4

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXLF6

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXLFD

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXLFH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXLFT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXMMX

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXMSE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXNER

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXNI3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXNI4

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXNI6

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXNID

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXNIH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXNIM

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXNIT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXNMR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXOMG

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXPFC

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXPRS

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXQPR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXRP2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXRPL

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSBI

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSIM

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYBOL

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYBPA

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYBPB

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYD2F

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYDE2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYDF3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYDF4

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYDF6

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYDFD

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYDFH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYDFT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYECF

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYGRB

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYLF3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYLF4

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYLF6

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYLFD

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYLFH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYLFT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYMSE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYNER

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYNI3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYNI4

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYNI6

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYNID

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYNIH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYNIM

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYNIT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYNMR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYOMG

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYPRS

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYQPR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYSBI

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYSIM

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYTER

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYTI2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYTI3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYTQI

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYV2F

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYVB2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYVC3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYVC4

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYVC6

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYVCD

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYVCH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYVCT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYVEE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYVIE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYVMO

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYVP2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYVPO

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYVPR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYVTR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXSYZF2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXTER

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXTI2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXTI3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXTQI

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXV2F

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXVB2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXVC3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXVC4

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXVC6

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXVCD

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXVCH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXVCT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXVEE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXVIE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXVMO

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXVP2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXVPO

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXVPR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXVTR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXZF2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXZIMPS

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXZIMPX

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LXZIMPXS

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LYRP2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LYRPL

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LYSIM

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LZFS

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%LZPSI

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%MMAX

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%MMX_MAXE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%MTRP2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%MTRPL

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NANTECH_D

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NANTICH_D

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NANTLH_D

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NBDATA

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NBDY

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NDEAEP

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NDFS

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NDMDF

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NDMFR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NDMMX

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NDMSE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NDPAEP

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NDPFC

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NDPFC_PRE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NEAEP

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NEFD0

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NEFDB

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NEFDP

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NEFDQ

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NEFDR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NEFDS

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NERBQ

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NFRBCODE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NFRBMODE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NGMAX_D

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NHECFT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NHFMX

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NIMDF

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NISSIM

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NISX

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NITSIM

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NITX

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NKAEP

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NKRBQ

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NLIM

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NMAPSIM

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NMDXKF

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NMIMP

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NMOMD

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NMUAEP

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NMURBQ

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NMZEFF

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NNSIM

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NPACK

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NPACK_RBQ

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NPAEP

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NPELDA

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NPMDF

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NPMSE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NPRBQ

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRHIS

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRHIX

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRHORBQ

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRIBOL

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRIBPA

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRIBPB

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRID2F

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRIDE2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRIDF3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRIDF4

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRIDF6

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRIDFD

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRIDFH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRIDFT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRIECF

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRIFD0

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRIFDB

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRIFDP

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRIFDQ

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRIFDR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRIFDS

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRIGRB

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRILF3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRILF4

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRILF6

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRILFD

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRILFH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRILFT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRIMSE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRINER

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRINI3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRINI4

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRINI6

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRINID

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRINIH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRINIM

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRINIT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRINMR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRIOMG

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRIPRS

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRIQPR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRISBI

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRISIM

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRITER

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRITI2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRITI3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRITQI

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRIV2F

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRIVB2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRIVC3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRIVC4

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRIVC6

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRIVCD

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRIVCH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRIVCT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRIVEE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRIVIE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRIVMO

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRIVP2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRIVPO

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRIVPR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRIVTR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRIZF2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRMFR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRMSE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRP2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRPL

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NRPSI

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSAWFLAG

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSC

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSC_FMT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSEVENT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSEVPER

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSEVTBL

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSEVWDS

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSFAST_D

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYBOL

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYBPA

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYBPB

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYD2F

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYDE2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYDF3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYDF4

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYDF6

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYDFD

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYDFH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYDFT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYECF

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYGRB

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYLF3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYLF4

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYLF6

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYLFD

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYLFH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYLFT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYMSE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYNER

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYNI3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYNI4

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYNI6

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYNID

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYNIH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYNIM

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYNIT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYNMR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYOMG

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYPRS

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYQPR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYSBI

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYSIM

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYTER

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYTI2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYTI3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYTQI

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYV2F

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYVB2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYVC3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYVC4

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYVC6

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYVCD

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYVCH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYVCT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYVEE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYVIE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYVMO

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYVP2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYVPO

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYVPR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYVTR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NSYZF2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NTABORT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NTAEP

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NTFD0

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NTFDB

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NTFDP

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NTFDQ

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NTFDR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NTFDS

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NTHFS

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NTIME1

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NTIME2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NTIMEC

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NTIMECA

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NTIMECB

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NTIMLH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NTIMNB

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NTIMRF

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NTIMRFF

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NTPSI

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NTRBQ

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NTSAW

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXBOL

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXBPA

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXBPB

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXD2F

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXDE2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXDF3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXDF4

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXDF6

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXDFD

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXDFH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXDFT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXECF

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXFD0

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXFDB

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXFDP

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXFDQ

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXFDR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXFDS

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXFS

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXGRB

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXLF3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXLF4

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXLF6

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXLFD

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXLFH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXLFT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXMMX

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXMSE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXNE0

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXNER

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXNI3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXNI4

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXNI6

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXNID

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXNIH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXNIM

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXNIT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXNM0

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXNMR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXOMG

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXPFC

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXPRS

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXQPR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXRP2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXRPL

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSBI

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSIM

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYBOL

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYBPA

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYBPB

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYD2F

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYDE2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYDF3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYDF4

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYDF6

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYDFD

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYDFH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYDFT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYECF

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYGRB

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYLF3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYLF4

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYLF6

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYLFD

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYLFH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYLFT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYMSE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYNER

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYNI3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYNI4

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYNI6

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYNID

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYNIH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYNIM

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYNIT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYNMR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYOMG

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYPRS

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYQPR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYSBI

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYSIM

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYTER

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYTI2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYTI3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYTQI

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYV2F

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYVB2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYVC3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYVC4

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYVC6

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYVCD

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYVCH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYVCT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYVEE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYVIE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYVMO

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYVP2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYVPO

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYVPR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYVTR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXSYZF2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXTE0

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXTER

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXTI2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXTI3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXTQI

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXV2F

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXVB2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXVC3

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXVC4

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXVC6

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXVCD

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXVCH

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXVCT

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXVEE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXVIE

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXVMO

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXVP2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXVPO

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXVPR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXVTR

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXZEF2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NXZF2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NYRP2

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NYRPL

      i=i+1
      i1=i2+1
      i2=i2+ivarsize(i)
      intbuf(i1:i2) = s%NZPSI

END SUBROUTINE trdatbuf_intmod_sub2
