#! /bin/csh -f
#
#  generate r8bloat (real*8 precision) from asymmetric single precision
#  TRANSP bloat routines
#
#    1.  create the working directory...
#
  if ( ! -d $TMPDIR_TR/r8bloat ) then
    mkdir $TMPDIR_TR/r8bloat
  endif
#
  echo " "
  echo 'r8bloat_lib.csh starting ("no match" message may follow...)'
#
  set wkdir = $TMPDIR_TR/r8bloat
  rm $wkdir/*.*
#
#    2.  fetch the sources into the working directory...
#
  set csdir = $CODESYSDIR/source
#
  cp $csdir/fftsub/*.for $wkdir
  cp $csdir/comput/bloata0.for $wkdir
  cp $csdir/comput/tk*.* $wkdir
  cp $csdir/comput/fpolar.for $wkdir
  cp $csdir/comput/lokjaca.for $wkdir
  cp $csdir/comput/sincos.for $wkdir
#
  echo " "
  echo "r8bloat_lib.csh:  work directory has r4 sources:"
  cd $wkdir
  ls
#
#    3.  get list of subroutine & function names; create substitution table
#
  foreach file ( *.for )
    fgtok $file >> r8bloat_names.tmp
  end
#
  set r8bloat_names = `cat r8bloat_names.tmp`
  rm r8bloat_names.tmp
#
  @ i = 0
  while ( $i < $#r8bloat_names ) 
    @ i++
    echo "$r8bloat_names[$i] r8$r8bloat_names[$i]" >> r8bloat_names.sub
  end
#
  echo " "
  echo "r8bloat_lib.csh:  substitution table succesfully generated."
#
# cat r8bloat_names.sub
#
  set ilist = ( *.inc )
#
  foreach file ( *.for ) 
    echo " converting:  $file ..."
    r8_convert $file r8$file r8bloat_names.sub
    @ j = 0
    while ( $j < $#ilist ) 
      @ j++
      supsub r8$file $ilist[$j] r8$ilist[$j]
    end
    rm $file
    if ( -f r8$file~ ) then
      rm r8$file~
    endif
  end
#
  foreach file ( $ilist ) 
    echo " converting include file:  $file ..."
    r8_convert_include $file r8$file r8bloat_names.sub
    rm $file
  end
#
  echo " "
  echo " conversion done; copying back to r8bloat source library directory."
#
  rm $csdir/r8bloat/*.for
  rm $csdir/r8bloat/*.inc
#
  cp $wkdir/*.* $csdir/r8bloat
#
  cd $csdir/r8bloat
  ls -l
