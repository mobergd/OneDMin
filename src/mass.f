      subroutine mass(x,y)

      implicit none

      character*2 x
      double precision y

      y=0.d0
      if (x.eq."H".or.x.eq."h") y= 1.007825d0
      if (x.eq."D".or.x.eq."d") y= 2.014102d0
      if (x.eq."C".or.x.eq."c") y=12.d0
      if (x.eq."N".or.x.eq."n") y=14.003074d0
      if (x.eq."O".or.x.eq."o") y=15.994915d0
      if (x.eq."P".or.x.eq."p") y=30.973762d0

      if (x.eq."He".or.x.eq."he".or.x.eq."HE") y=4.002603d0
      if (x.eq."Ne".or.x.eq."ne".or.x.eq."NE") y=19.992440d0
      if (x.eq."Ar".or.x.eq."ar".or.x.eq."AR") y=39.962383d0
      if (x.eq."Kr".or.x.eq."kr".or.x.eq."KR") y=83.911507d0

      if (x.eq."Si".or.x.eq."si".or.x.eq."SI") y=27.9769265d0
      if (x.eq."Ti".or.x.eq."ti".or.x.eq."TI") y=47.947947d0

      if (x.eq."F".or.x.eq."f")                y=18.998403d0
      if (x.eq."Cl".or.x.eq."cl".or.x.eq."CL") y=34.968853d0

c special
      if (x.eq."Aa".or.x.eq."aa".or.x.eq."AA") y=39.962383d0
      if (x.eq."Hx".or.x.eq."hx".or.x.eq."HX") y=4.002603d0
      if (x.eq."Ha".or.x.eq."he".or.x.eq."HA") y=4.002603d0
      if (x.eq."N2".or.x.eq."n2") y=14.003074d0
      if (x.eq."Na".or.x.eq."na") y=14.003074d0
      if (x.eq."Np".or.x.eq."np") y=14.003074d0
      if (x.eq."LJ".or.x.eq."LJ") y=14.003074d0
      if (x.eq."Ca".or.x.eq."ca") y=12.d0
      if (x.eq."Ho".or.x.eq."ho") y= 1.007825d0
      if (x.eq."Oh".or.x.eq."oh") y=15.994915d0

      if (y.eq.0.d0) then
        write(6,*)"Don't know the mass of ",x
        write(6,*)"Please add it to the subroutine 'mass'"
        stop
      endif

      return

      end

