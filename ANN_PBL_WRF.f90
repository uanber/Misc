MODULE ANN_PBL_WRF

implicit none

CONTAINS

public :: relu, relu_prime

SUBROUTINE ANN_INPUTS (pbl_input_mean, pbl_input_scale, pbl_diagnostic_mean, pbl_diagnostic_scale, &
              state_mean, state_scale, sl_real_mean, rt_real_mean, sl_latent_to_real, &
              rt_latent_to_real, w_latent_to_real, rcld_latent_to_real, rrain_latent_to_real, &
              cld_latent_to_real, pbl_encoder_W, pbl_encoder_b, pbl_hidden_W, pbl_hidden_b, &
              pbl_tend_decoder_W, pbl_tend_decoder_b, pbl_diag_decoder_W, pbl_diag_decoder_b, &
              !Real WRF variables
              
              sl_real, th_phy, rt_real, sl_adv_real, rt_adv_real, w_real, sl_domain_top, &
              rt_domain_top, lh, shf, tsk, SWDNT, swdn_tod, PSFC, &
              cldmid, cldhigh, qr, qc, qv, &
              ph, phb, & ! geopotential hieght and its perturbation 
              !theta_T
              ims,ime,jms,jme,kms,kme, znu)
              
              
              
              real, parameter:: g= 9.81
              real, parameter:: Cpd= 1005.
              real, parameter:: Rd= 287.
              real, parameter:: Lv= 2.5e6
              
              integer , parameter:: cld=10
              integer , parameter:: rcld=10
              integer , parameter:: rrain=2
              integer , parameter:: rt=8
              integer , parameter:: rt_adv=8
              integer , parameter:: sl=9
              integer , parameter:: sl_adv=9
              integer , parameter:: sl_rad_clr=9
              integer , parameter:: w=5
              
              

!integer, INTENT(IN) ::  npz
!real, intent(in), dimension(is:ie,js:je,npz+1) :: phalf, zhalf
INTEGER, INTENT(IN) ::                                             &
     &           ids,ide,jds,jde,kds,kde                              &
     &          ,ims,ime,jms,jme,kms,kme                              &
     &          ,kts,kte,num_tiles

real, intent(in), dimension(kms:kme) :: znu ! hieght

real, intent(in), dimension(33) :: pbl_input_mean, pbl_input_scale
real, intent(in), dimension(34) :: pbl_diagnostic_mean, pbl_diagnostic_scale
real, intent(in), dimension(17) :: state_mean, state_scale
real, intent(in), dimension(20) :: sl_real_mean, rt_real_mean

real, intent(in), dimension(sl, 20) :: sl_latent_to_real
real, intent(in), dimension(rt, 20) :: rt_latent_to_real
real, intent(in), dimension(w, 20) :: w_latent_to_real

real, intent(in), dimension(rcld, 20) :: rcld_latent_to_real
real, intent(in), dimension(rrain, 20) :: rrain_latent_to_real
real, intent(in), dimension(cld, 20) :: cld_latent_to_real


real, intent(in), dimension(33, 128) :: pbl_encoder_W
real, intent(in), dimension(128) :: pbl_encoder_b
real, intent(in), dimension(128,128) :: pbl_hidden_W
real, intent(in), dimension(128) :: pbl_encoder_b
real, intent(in), dimension(128,17) :: pbl_tend_decoder_W
real, intent(in), dimension(17) :: pbl_tend_decoder_b
real, intent(in), dimension(128,34) :: pbl_diag_decoder_W
real, intent(in), dimension(34) :: pbl_diag_decoder_b


!integer, INTENT(IN) ::  npz
real, intent(INOUT), dimension(ims:ime,kms:kme,jms:jme) :: sl_real, rt_real, w_real
real, intent(in), dimension(ims:ime,kms:kme+1,jms:jme) ::  w ! model vertical wind (defined at model grid)
real, intent(in), dimension(ims:ime,kms:kme,jms:jme) :: th_phy ! potential temperature perturbation
!real, intent(INOUT), dimension(ims:ime,kms:kme,jms:jme) :: theta_T ! potential temperature 

real, intent(in), dimension(ims:ime,jms:jme) :: hfx, lh ! sensible and latent heat flux in w/m^2
real, intent(in), dimension(ims:ime,jms:jme) :: tsk ! sfc temp.

real, intent(in), dimension(ims:ime,jms:jme) :: SWDNT ! downwelling SW at TOA
real, intent(in), dimension(ims:ime,jms:jme) :: PSFC ! sfc pressure 


real, intent(in), dimension(ims:ime,kms:kme,jms:jme) :: qr ! 3D rain water mixing ratio. 
                                                           !For top of domain just specify some level at 3km

real, intent(in), dimension(ims:ime,kms:kme,jms:jme) :: qc ! 3D cloud water mixing ratio. 
                                                           
real, intent(in), dimension(ims:ime,kms:kme,jms:jme) :: qv ! 3D water vapor mixing ratio. 
                                                           

              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
! calculating liquid water static energy sl_real and rt_real and w_real!             
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do k=1, km
   do i=ims,ime
      do j=jms,jme
      
          sl_real(i,k,j) = Cpd* th_phy(i,k,j)+300 - Lv * qc(i,k,j) + g * (PHB(i,k,j) + PH(i,k,j))/9.81 ! last term is g*z
          rt_real(i,k,j) = qv(i,k,j) + qc(i,kj)
          w_real(i,j,k) = 0.5*(W(i,j,k) + W(i,j,k+1))
          
      enddo
    enddo
enddo





SUBROUTINE ANN_INPUTS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! defining linear activation function and its derivative !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

pure function relu(x) result(res)
    !! Rectified Linear Unit (RELU) activation function.
    real(rk), intent(in) :: x(:)
    real(rk) :: res(size(x))
    res = max(0., x)
  end function relu

  pure function relu_prime(x) result(res)
    ! First derivative of the REctified Linear Unit (RELU) activation function.
    real(rk), intent(in) :: x(:)
    real(rk) :: res(size(x))
    where (x > 0)
      res = 1
    elsewhere
      res = 0
    end where
  end function relu_prime
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

