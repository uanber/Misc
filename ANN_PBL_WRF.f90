MODULE ANN_PBL_WRF

implicit none

CONTAINS

public :: relu, latent

SUBROUTINE ANN_INPUTS (pbl_input_mean, pbl_input_scale, pbl_diagnostic_mean, pbl_diagnostic_scale, &
              state_mean, state_scale, sl_real_mean, rt_real_mean, sl_latent_to_real, &
              rt_latent_to_real, w_latent_to_real, rcld_latent_to_real, rrain_latent_to_real, &
              cld_latent_to_real, pbl_encoder_W, pbl_encoder_b, pbl_hidden_W, pbl_hidden_b, &
              pbl_tend_decoder_W, pbl_tend_decoder_b, pbl_diag_decoder_W, pbl_diag_decoder_b, &
              !Real WRF variables
              
              th_phy, sl_adv_real, rt_adv_real, & ! what's this?
              rt_domain_top, lh, shf, tsk, SWDNT, swdn_tod, PSFC, &
              CLDFRA, qr, qc, qv, &
              ph, phb, & ! geopotential hieght and its perturbation 
              !theta_T
              sl_avg, rt_avg, w_avg, qr_avg, cld_avg, & ! domain mean profiles sl, rt, w, and qr (rain water mixing ratio needed for qr(d_top)
              sst_avg, shf_avg, lh_avg, wsdnt_avg, psfc_avg, & ! domain mean 2d variables
              ids,ide,jds,jde,kds,kde, &
              ims,ime,jms,jme,kms,kme, &
              kts,kte,num_tiles, znu, z, &
              d_top, d_low, d_mid, d_high )
              
              
              
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
              
!
USE module_dm, ONLY: wrf_dm_sum_real ! import the domain avg function
USE module_init_utilities, ONLY : interp_0 ! import interpolation function

!integer, INTENT(IN) ::  npz
!real, intent(in), dimension(is:ie,js:je,npz+1) :: phalf, zhalf
INTEGER, INTENT(IN) ::                                             &
     &           ids,ide,jds,jde,kds,kde                              &
     &          ,ims,ime,jms,jme,kms,kme                              &
     &          ,kts,kte,num_tiles

!!!!!!!
INTEGER, INTENT(IN) ::  d_top ! top of the BPL say 3 km to be specified in the namelist
INTEGER, INTENT(IN) ::  d_low, d_mid, d_high ! index of low, mid, and high clouds to be specified in the namelist

!!!!!!!


real, intent(in), dimension(kms:kme) :: znu ! hieght

REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN   ) :: z ! this is the height in x,y,z

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
!
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

real, intent(in), dimension(ims:ime,kms:kme,jms:jme) :: cldfra ! cloud fraction

                                                           
real, intent(INOUT), dimension(kms:kme) :: sl_avg, rt_avg, w_avg, qr_avg, cld_avg, ! domain avg sl, rt, and w and qr ! can also be local                                                      
real, intent(INOUT)                     :: sst_avg, shf_avg, lh_avg, swdnt_avg, psfc_avg ! domain mean

real ::  no_point
real ::  sl_sum, rt_sum, w_sum, qr_sum, cld_sum, z_sum, sst_sum, shf_sum, lh_sum, swdnt_sum, psfc_sum

! local
real, DIMENSION( ims:ime , jms:jme ) :: sl_2d, rt_2d, w_2d , qr_2d, cld_2d
real, DIMENSION( kms:kme ) :: z_avg 

! local  
real, dimension(ims:ime,kms:kme,jms:jme) :: sl_real, rt_real, w_real
real :: cldlow, cldmid, cldhigh, rrain_domain_top, swdn_tod, sl_domain_top, rt_domain_top, rrain_domain_top
              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
! calculating liquid water static energy sl_real and rt_real and w_real!             
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do k=1, km
   do i=ims,ime
      do j=jms,jme
      
          sl_real(i,k,j) = Cpd* th_phy(i,k,j)+300 - Lv * qc(i,k,j) + g * (PHB(i,k,j) + PH(i,k,j))/9.81 !last term is g*z (could also use z)
          rt_real(i,k,j) = qv(i,k,j) + qc(i,k,j)
          w_real(i,k,j) = 0.5*(W(i,k,j) + W(i,k+1,j))
          
      enddo
    enddo
enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
! calculating profile averages !             
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

no_points = float((ide-ids)*(jde-jds)) ! number of grid points
 
DO k=kts,kde-1
      !DO k=kts,min(kte+1,kde)
         sl_2d = 0.
         rt_2d = 0.
         w_2d = 0.
         qr_2d =0.
         cld_2d = 0.
         
         sl_sum = 0.0
         rt_sum=0.0
         w_sum = 0.0
         qr_sum = 0.0
         cld_sum = 0.0
         
         z_sum = 0.0

         DO ij = 1 , num_tiles;
            DO j=j_start(ij),j_end(ij); DO i=i_start(ij),i_end(ij)
            
               sl_2d(i,j) = sl_real(i,k,j)
               sl_sum = sl_sum + sl_2d(i,j)

               rt_2d(i,j) = rt_real(i,k,j)
               rt_sum = rt_sum + rt_2d(i,j)

               w_2d(i,j) = w_real(i,k,j)
               w_sum = w_sum + w_2d(i,j)
               
               qr_2d(i,j) = qr(i,k,j)
               qr_sum = qr_sum + qr_2d(i,j)
               
               cld_2d(i,j) = cldfra(i,k,j)
               cld_sum = cld_sum + cld_2d(i,j)
               
               z_sum = z_sum + z(i,k,j)

            ENDDO; ENDDO
         ENDDO
         

         !domain mean sl:
         sl_avg(k) = wrf_dm_sum_real ( sl_sum )
         sl_avg(k) = sl_avg(k) / no_points
         
         !domain mean rt:
         rt_avg(k) = wrf_dm_sum_real ( drt_sum )
         rt_avg(k) = rt_avg(k) / no_points
         
         !domain mean w:
         w_avg(k) = wrf_dm_sum_real ( w_sum )
         w_avg(k) = w_avg(k) / no_points
         
         !domain mean qr:
         qr_avg(k) = wrf_dm_sum_real ( qr_sum )
         qr_avg(k) = qr_avg(k) / no_points
         
         !domain mean cldfra:
         cld_avg(k) = wrf_dm_sum_real ( cld_sum )
         cld_avg(k) = cld_avg(k) / no_points

         !domain mean z (useful for interpolation)
         z_avg(k) = wrf_dm_sum_real ( z_sum )
         z_avg(k) = z_avg(k) / no_points
         
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
! calculating 2D variables averages !             
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

sst_sum = 0.
shf_sum = 0.
lh_sum = 0.
swdnt_sum = 0.0
psfc_sum = 0.0


DO ij = 1 , num_tiles;
   DO j=j_start(ij),j_end(ij); DO i=i_start(ij),i_end(ij)
   
      sst_sum = sst_sum + tsk(i,j)
               
      shf_sum = shf_sum + shf(i,j)

      lh_sum = lh_sum + lh(i,j)
               
      swdnt_sum = swdnt_sum + swdnt(i,j)
               
      psfc_sum = psfc_sum + psfc(i,j)
               
      ENDDO; ENDDO
ENDDO


!domain mean :
       
sst_avg = sst_sum / no_points

shf_avg = shf_sum / no_points
         
lh_avg = lh_sum / no_points

swdnt_avg = swdnt_sum / no_points

psfc_avg = psfc_sum / no_points


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! variables at low, mid, high, and domain top levels !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

cldlow = cld_avg(d_low)
cldmid = cld_avg(d_mid)
cldhigh = cld_avg(d_high)

swdn_tod = swdnt(d_top)

rrain_domain_top = qr_avg(d_top)

sl_domain_top = sl_avg(d_top)

rt_domain_top = rt_avg (d_top)

rrain_domain_top = ar_avg(d_top)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Concatenate variables in one input array  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

input_array = [sl_avg, rt_avg, w_avg, sl_domain_top, rt_domain_top, & 
               lh_avg, shf_avg, sst_avg, cldmid, cldhigh, &
               swdnt_avg, swdn_tod, psfc_avg, rrain_domain_top]





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Interpolate the weights, biases, and inputs into the models grid !! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

do k=kts,kte
   th_scm_target_modellevels(k) = interp_0(th_scm_target, z_force, zzz_avg(k), num_force_layers) 
enddo





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Transform variables from physical space to to latent space !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

sl_latent = matmul(sl_avg(:) - sl_real_mean(:) , transpose(sl_latent_to_real) )
rt_latent = matmul(rt_avg(:) - sl_real_mean(:) , transpose(rt_latent_to_real) )
w_latent = matmul(w_avg(:) - sl_real_mean(:) , transpose(w_latent_to_real) )









SUBROUTINE ANN_INPUTS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! defining linear activation function                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

pure function relu(x) result(res)
    !! Rectified Linear Unit (RELU) activation function.
    real , intent(in) :: x(:)
    real :: res(size(x))
    res = max(0., x)
  end function relu

  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Function to transform variables from physical space to latent space !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 pure function latent(x, x_c) result(x_latent)
    
    real, intent(in) :: x(:)
    real, intent(in) :: x_c(:)
    
    real, intent(out) :: x_latent(:)
    
    !sl_latent = tf.tensordot(sl_real - sl_real_mean, tf.transpose(sl_latent_to_real), [0, 0])
     x_latent = matmul(x , x_c)
    
  end function latent

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Function to transform variables from latent space to physical space !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 pure function physical(x, x_c) result(x_physical)
    
    real, intent(in) :: x(:)
    real, intent(in) :: x_c(:)
    
    real, intent(out) :: x_latent(:)
    
     x_physical = matmul(x_c , x)
    
  end function physical

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Function to concatenate arrays for ANN inputs  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


function Concat_Arrays
implicit none
 
  ! Note: in Fortran 90 you must use the old array delimiters (/ , /)
  real, dimension(:) :: a 
  real, dimension(:) :: b 
  real, dimension(:), allocatable :: c
 
  allocate(c(size(a)+size(b)))
  c(1 : size(a)) = a
  c(size(a)+1 : size(a)+size(b)) = b
  write(*,*) c
 
  ! alternative
  d = [a, b] ! (/a, b/)
  write(*,*) d
 
  deallocate(c)
  deallocate(d)
 
end function Concat_Arrays



! flatten 2D array into 1D : new_array = pack(old_array,.true.)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Forward Prop. Subroutine  !!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine feedforward (a, z, w, b)

real, intent(in), dimensions(:,:) :: w
real, intent(in), dimensions(:) :: b
real, intent(inout), dimensions(:) :: a, z


do n = 1, 2 !size(layers)
        z = matmul(a, transpose(w(:,n))) +  b
        a = relu(z)
enddo
      
end subroutine feedforward




