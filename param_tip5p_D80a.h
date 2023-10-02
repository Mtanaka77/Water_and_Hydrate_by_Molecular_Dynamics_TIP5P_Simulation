!  param_tip5p_D80a.h
!
      integer(C_INT) np0,nq0,npq0,npq5,nq1,            &
                     ip0,mesh,mintpol,mx,my,mz,        &
                     mx1,my1,mz1,kstart,brillouin
      logical        if_xyz1,if_xyz2,if_obsv
!
! /home2, or /lv01 for photon
      character     praefixs*29,praefixc*24,praefixe*24, &
                    praefixi*24,suffix2*2,suffix1*2,suffix0*1
!     character*4   praefixs*31,praefixc*29,praefixe*29, &
!                   praefixi*29,suffix2*2,suffix1*2,suffix0*1
!
      real(C_DOUBLE) epsilon
!     real(C_DOUBLE) t_init,t_wipe_sta,t_wipe_end
!
      parameter  (epsilon=88.d0/88.d0)  ! 273 K
!     parameter  (t_init=1000.d0,t_wipe_sta=1700.d0,  &
!                 t_wipe_end=4700.d0)
!
      parameter  (kstart=0,suffix2='0a', & ! 120a is created
                           suffix1='0a', & ! TIP501__0 
                           suffix0='0')    ! 
!     parameter  (kstart=2,suffix2='0b', & ! 120b, kstart=2
!                          suffix1='0a', & ! TIP501__0
!                          suffix0='0')    ! 
!     parameter  (kstart=2,suffix2='0d', & ! 120d, kstart=2
!                          suffix1='0c', & ! TIP501__0
!                          suffix0='0')    ! 
!     parameter  (kstart=1,suffix2='0f', & ! 120f, kstart=1
!                          suffix1='0e', & ! TIP501__1
!                          suffix0='1')    ! 
!     parameter  (kstart=3,suffix2='0g', & ! 120g, kstart=3
!                          suffix1='0f', & ! TIP501__1, exc>0
!                          suffix0='1')    ! 
! /home2, /lv01                          +++++++ short 
      parameter (praefixs='/home/tanakam/MPI_wat5/TIP580',  & ! LXwat3
                 praefixi='/data/sht/tanakam/tip580',      &
                 praefixc='/data/sht/tanakam/tip580',      &
                 praefixe='/data/sht/tanakam/tip580')
      parameter (if_xyz1=.true., if_xyz2=.false.,          &
                 if_obsv=.false.)  !! if .false,, save to a file 
!
! moldyn  L > 32 Ang
!  Total number is nq0; for mh3: nq0=6210,np0=216
!                              +++++ null
      parameter  (nq0=6210,np0=216) ! 6210 atoms by 5-body CS1
!
      parameter  (npq5=nq0+np0)
      parameter  (npq0=nq0/5+np0)
      parameter  (nq1=nq0/5)
!
      parameter  (ip0=3,mesh=32,mintpol=4*50048)
      parameter  (brillouin=1)
      parameter  (mx=mesh,my=mesh,mz=mesh,mx1=mx+1,my1=my+1,mz1=mz+1)
!
!  forces.h
!
      real(C_DOUBLE) AA1,AA2,AA3,AA4,AA5,PPP,sqrtpi
      parameter ( AA1 = 0.254829592d0, AA2 = -0.284496736d0 )
      parameter ( AA3 = 1.421413741d0, AA4 = -1.453152027d0 )
      parameter ( AA5 = 1.061405429d0, PPP =  0.3275911d0   )
      parameter ( sqrtpi = sqrt(3.14159265358979323846d0)   ) 

