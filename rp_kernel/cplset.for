      subroutine cplset_exec
      use cplotr_mod
      use mfblok_mod
C
C  SET AXES DEFAULTS - ALL LINEAR
C
      NAXISC = 1
C
C  SET SCALE DEFAULTS - ALL AUTOSCALE
C
      NSCALC = 1
C
C  default disk/dir string lengths
C
      lfdir = 0
      lfdisk = 0
C
C  clear "extra" run pointers
C
      nrun_x = 0
      krun_x = 0
      lrun_x = 0
C
      nu_mks = 0
C
C  LUNs for extra open MF.PLN files
C    caution, RPLOT uses 77 & 88, 89, 90, 91!
C    UREAD uses 5,6, 10--20 (normally)
C    LUN mgmt is risky in the current state of the software...
C
      MFLUNI_X = (/78,79,80,81,82,83,84,85,86,87/)
C
C  MDSplus flag default:  .FALSE.
C
      NLMDS = .FALSE.
      NLMDS_X = .FALSE.
C
C  NetCDF flag default:  .FALSE.
C
      nlcdf = .false.
      nlcdf_x = .false.
C
      mds_cache = .false.
      lun_tf = 90
      lun_nf = 91
      lun_mf = 92
C
      end
