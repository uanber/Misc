MODULE ANN_PBL_WRF

CONTAINS

SUBROUTINE ANN_INPUTS(pbl_input_mean, pbl_input_scale, pbl_diagnostic_mean, pbl_diagnostic_scale, &
              state_mean, state_scale, sl_real_mean, rt_real_mean, sl_latent_to_real, &
              rt_latent_to_real, w_latent_to_real, rcld_latent_to_real, rrain_latent_to_real, &
              cld_latent_to_real, pbl_encoder_W, pbl_encoder_b, pbl_hidden_W, pbl_hidden_b, &
              pbl_tend_decoder_W, pbl_tend_decoder_b, pbl_diag_decoder_W, pbl_diag_decoder_b, &
              !Real WRF variables
              
              sl_real, rt_real, sl_adv_real, rt_adv_real, w_real, sl_domain_top, &
              rt_domain_top, lhf, shf, sst, swdn_toa, swdn_tod, p_surface, &
              cldmid, cldhigh, rrain_domain_top)
              
              
              
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


integer, INTENT(IN) ::  npz
real, intent(in), dimension(is:ie,js:je,npz+1) :: sl_real, rt_real

              
              
              
              
