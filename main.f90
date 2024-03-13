
!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

program fd_1d
!----------------------------------------------------------------------------------------------------------!
  use constants, only: iow
  use eos, only: eos_df_drho, eos_type
  use write_helper, only: adjl
  use flags, only: F_bc_dirichlet_eq_0, F_bc_dirichlet_eq_1, F_sphere, F_incompressible
  use parser_vars, only: beta, k_gr, delta,      &
                        & mxa_kind, mxb_kind, ghi_kind, glo_kind, &
                        & bc_lo_mxa, bc_lo_mxb, bc_lo_grafted, bc_hi_mxa, bc_hi_mxb, bc_hi_grafted, &
                        & chainlen_mxa, chainlen_mxb, chainlen_glo, chainlen_ghi,  &
                        & ns_mxa, ns_mxb, ns_glo, ns_ghi, &
                        & exist_mxa, exist_mxb, exist_glo, exist_ghi,            &
                        & Rg2_per_mon_mxa, Rg2_per_mon_mxb, Rg2_per_mon_glo, Rg2_per_mon_ghi, &
                        & edwards_solver, linear_solver, geometry,    &
                        & check_stability_every, compute_every, field_every, frac, nx, &
                        & gnode_lo, gnode_hi, rho_seg_bulk,&
                        & gdens_lo, gdens_hi, max_iter, max_wa_error, square_gradient, thermo_every, &
                        & chi12
  use arrays, only: qmxa, qmxb, qglo, qghi, qglo_aux, qghi_aux, &
                        & qfinal_mxa, qfinal_mxb, qfinal_glo, qfinal_ghi, qfinal_glo_aux, qfinal_ghi_aux, &
                        & dir_nodes_id, dir_nodes_rdiag, n_dir_nodes, &
                        & phi_mxa, phi_mxb, phi_glo, phi_ghi, phi_kd1, phi_kd2, phi_tot, &
                        & dphi_dr, d2phi_dr2, &
                        & dx, coeff_nx, &
                        & coeff_ns_mxa, coeff_ns_mxb, coeff_ns_glo, coeff_ns_ghi, &
                        & ds_mxa, ds_mxb, ds_glo, ds_ghi, &
                        & Ufield, &
                        & wa_kd1, wa_bulk_kd1, wa_ifc_kd1, wa_ifc_new_kd1, wa_ifc_backup_kd1,     &
                        & wa_kd2, wa_bulk_kd2, wa_ifc_kd2, wa_ifc_new_kd2, wa_ifc_backup_kd2,     &
                        & wa_ifc_mxa, wa_ifc_mxb, wa_ifc_glo, wa_ifc_ghi, &
                        & surface_area, rr, irr, layer_area,pressure
!----------------------------------------------------------------------------------------------------------!
  implicit none
!----------------------------------------------------------------------------------------------------------!
  logical :: convergence = .false., restart = .false., restore = .false.

  integer :: iter, jj, ii, tt

 
  real(8) :: get_nchains
  real(8) :: wa_error_new = 1.d10, wa_error_old = 1.d10
  real(8) :: qinit_lo = 0.d0, qinit_hi = 0.d0
  real(8) :: free_energy = 0.d0
  real(8) :: nchglo = 0.d0, nchghi = 0.d0, nch_mxa = 0.d0, nch_mxb = 0.d0
!----------------------------------------------------------------------------------------------------------!
  open (unit=iow, file="o.log")

  call parser
  call init_scf_params
  call init_arrays
  call init_mesh
  call init_geom
  call init_solid
  call init_field
   

  write (iow, '(A85)') adjl("---------------------------------BEGIN THE SIMULATION--------------------------------", 85)
  write (*, '(A85)') adjl("---------------------------------BEGIN THE SIMULATION--------------------------------", 85)

  write (iow, '(2X,A15,10(2X,A15))') "Iteration", "energy (mN/m)", "error (k_B T)", "gdens_lo", 'gdens_hi', "nch_m", &
  &                               "nch_mxb", "nch_glo", "nch_ghi", "fraction"
  write (*, '(2X,A15,10(2X,A15))') "Iteration", "energy (mN/m)", "error (k_B T)", "gdens_lo", 'gdens_hi', "nch_m", &
  &                               "nch_mxb", "nch_glo", "nch_ghi", "fraction"
  
  do iter = 0, max_iter

    if (exist_mxa) then
      !set the dirichlet boundary conditions for matrix chains
      n_dir_nodes = 0
      !dirichlet lower bound
      if (bc_lo_mxa .eq. F_bc_dirichlet_eq_0) then
        dir_nodes_id(n_dir_nodes) = 0
        dir_nodes_rdiag(n_dir_nodes) = 0.d0
        n_dir_nodes = n_dir_nodes + 1
      else if (bc_lo_mxa .eq. F_bc_dirichlet_eq_1) then
        dir_nodes_id(n_dir_nodes) = 0
        dir_nodes_rdiag(n_dir_nodes) = 1.0d0
        n_dir_nodes = n_dir_nodes + 1
      end if
      !dirichlet upper bound
      if (bc_hi_mxa .eq. F_bc_dirichlet_eq_0) then
        dir_nodes_id(n_dir_nodes) = nx
        dir_nodes_rdiag(n_dir_nodes) = 0.d0
        n_dir_nodes = n_dir_nodes + 1
      else if (bc_hi_mxa .eq. F_bc_dirichlet_eq_1) then
        dir_nodes_id(n_dir_nodes) = nx
        dir_nodes_rdiag(n_dir_nodes) = 1.0d0
        n_dir_nodes = n_dir_nodes + 1
      end if

      if (geometry .eq. F_sphere) then
        do ii = 0, n_dir_nodes - 1
          dir_nodes_rdiag(ii) = dir_nodes_rdiag(ii)*rr(dir_nodes_id(ii))
        end do
      end if

      !matrix chains
      do ii = 0, nx
        qmxa(ii, 1) = 1.d0
        qfinal_mxa(ii, 0) = 1.d0
      end do

      if (geometry .eq. F_sphere) then
        do ii = 0, nx
          qmxa(ii, 1) = qmxa(ii, 1)*rr(ii)
          qfinal_mxa(ii, 0) = qfinal_mxa(ii, 0)*rr(ii)
        end do
      end if

      call solver_edwards(bc_lo_mxa, bc_hi_mxa, n_dir_nodes, dir_nodes_id, dir_nodes_rdiag, &
&                           Rg2_per_mon_mxa, nx, ns_mxa, dx, ds_mxa, edwards_solver,      &
&                           linear_solver, wa_ifc_mxa, qmxa, qfinal_mxa)

      if (geometry .eq. F_sphere) then
        do tt = 0, ns_mxa
          do ii = 0, nx
            qfinal_mxa(ii, tt) = qfinal_mxa(ii, tt)*irr(ii)
          end do
        end do
      end if

      call contour_convolution(chainlen_mxa, nx, ns_mxa, coeff_ns_mxa, &
&                                qfinal_mxa, qfinal_mxa, phi_mxa)
    end if

    if (exist_mxb) then
      !set the dirichlet boundary conditions for matrix chains
      n_dir_nodes = 0
      !dirichlet lower bound
      if (bc_lo_mxb .eq. F_bc_dirichlet_eq_0) then
        dir_nodes_id(n_dir_nodes) = 0
        dir_nodes_rdiag(n_dir_nodes) = 0.d0
        n_dir_nodes = n_dir_nodes + 1
      else if (bc_lo_mxb .eq. F_bc_dirichlet_eq_1) then
        dir_nodes_id(n_dir_nodes) = 0
        dir_nodes_rdiag(n_dir_nodes) = 1.0d0
        n_dir_nodes = n_dir_nodes + 1
      end if
      !dirichlet upper bound
      if (bc_hi_mxb .eq. F_bc_dirichlet_eq_0) then
        dir_nodes_id(n_dir_nodes) = nx
        dir_nodes_rdiag(n_dir_nodes) = 0.d0
        n_dir_nodes = n_dir_nodes + 1
      else if (bc_hi_mxb .eq. F_bc_dirichlet_eq_1) then
        dir_nodes_id(n_dir_nodes) = nx
        dir_nodes_rdiag(n_dir_nodes) = 1.0d0
        n_dir_nodes = n_dir_nodes + 1
      end if

      if (geometry .eq. F_sphere) then
        do ii = 0, n_dir_nodes - 1
          dir_nodes_rdiag(ii) = dir_nodes_rdiag(ii)*rr(dir_nodes_id(ii))
        end do
      end if

      !matrix chains
      do ii = 0, nx
        qmxb(ii, 1) = 1.d0
        qfinal_mxb(ii, 0) = 1.d0
      end do

      if (geometry .eq. F_sphere) then
        do ii = 0, nx
          qmxb(ii, 1) = qmxb(ii, 1)*rr(ii)
          qfinal_mxb(ii, 0) = qfinal_mxb(ii, 0)*rr(ii)
        end do
      end if

      call solver_edwards(bc_lo_mxb, bc_hi_mxb, n_dir_nodes, dir_nodes_id, dir_nodes_rdiag, &
&                           Rg2_per_mon_mxb, nx, ns_mxb, dx, ds_mxb, edwards_solver,      &
&                           linear_solver, wa_ifc_mxb, qmxb, qfinal_mxb)

      if (geometry .eq. F_sphere) then
        do tt = 0, ns_mxb
          do ii = 0, nx
            qfinal_mxb(ii, tt) = qfinal_mxb(ii, tt)*irr(ii)
          end do
        end do
      end if

      call contour_convolution(chainlen_mxb, nx, ns_mxb, coeff_ns_mxb, &
&                                qfinal_mxb, qfinal_mxb, phi_mxb)
    end if

    !edwards diffusion for grafted chains in the lower boundary
    if (exist_glo) then

      !set the dirichlet boundary conditions
      n_dir_nodes = 0
      !dirichlet lower bound
      if (bc_lo_grafted .eq. F_bc_dirichlet_eq_0) then
        dir_nodes_id(n_dir_nodes) = 0
        dir_nodes_rdiag(n_dir_nodes) = 0.d0
        n_dir_nodes = n_dir_nodes + 1
      else if (bc_lo_grafted .eq. F_bc_dirichlet_eq_1) then
        dir_nodes_id(n_dir_nodes) = 0
        dir_nodes_rdiag(n_dir_nodes) = 1.0d0
        n_dir_nodes = n_dir_nodes + 1
      end if
      !dirichlet upper bound
      if (bc_hi_grafted .eq. F_bc_dirichlet_eq_0) then
        dir_nodes_id(n_dir_nodes) = nx
        dir_nodes_rdiag(n_dir_nodes) = 0.d0
        n_dir_nodes = n_dir_nodes + 1
      else if (bc_hi_grafted .eq. F_bc_dirichlet_eq_1) then
        dir_nodes_id(n_dir_nodes) = nx
        dir_nodes_rdiag(n_dir_nodes) = 1.0d0
        n_dir_nodes = n_dir_nodes + 1
      end if

      if (geometry .eq. F_sphere) then
        do ii = 0, n_dir_nodes - 1
          dir_nodes_rdiag(ii) = dir_nodes_rdiag(ii)*rr(dir_nodes_id(ii))
        end do
      end if

      !glo_aux
      do ii = 0, nx
        qglo_aux(ii, 1) = 1.d0
        qfinal_glo_aux(ii, 0) = 1.d0
      end do

      if (geometry .eq. F_sphere) then
        do ii = 0, nx
          qglo_aux(ii, 1) = qglo_aux(ii, 1)*rr(ii)
          qfinal_glo_aux(ii, 0) = qfinal_glo_aux(ii, 0)*rr(ii)
        end do
      end if

      call solver_edwards(bc_lo_grafted, bc_hi_grafted, n_dir_nodes, dir_nodes_id, dir_nodes_rdiag, &
&                           Rg2_per_mon_glo, nx, ns_glo, dx, ds_glo, edwards_solver,        &
&                           linear_solver, wa_ifc_glo, qglo_aux, qfinal_glo_aux)

      if (geometry .eq. F_sphere) then
        do tt = 0, ns_glo
          do ii = 0, nx
            qfinal_glo_aux(ii, tt) = qfinal_glo_aux(ii, tt)*irr(ii)
          end do
        end do
      end if

      ! grafted chains
      do ii = 0, nx
        qglo(ii, 1) = 0.d0
        qfinal_glo(ii, 0) = 0.d0
      end do

      delta = (1.0/(dx(gnode_lo)*layer_area(gnode_lo)))*1e+30 !1/m3
      qinit_lo = chainlen_glo*delta*(gdens_lo*surface_area) &
&                     /(rho_seg_bulk*qfinal_glo_aux(gnode_lo, ns_glo))

      qglo(gnode_lo, 1) = qinit_lo
      qfinal_glo(gnode_lo, 0) = qinit_lo

      if (geometry .eq. F_sphere) then
        do ii = 0, nx
          qglo(ii, 1) = qglo(ii, 1)*rr(ii)
          qfinal_glo(ii, 0) = qfinal_glo(ii, 0)*rr(ii)
        end do
      end if

      call solver_edwards(bc_lo_grafted, bc_hi_grafted, n_dir_nodes, dir_nodes_id, dir_nodes_rdiag, &
&                           Rg2_per_mon_glo, nx, ns_glo, dx, ds_glo, edwards_solver,        &
&                           linear_solver, wa_ifc_glo, qglo, qfinal_glo)

      if (geometry .eq. F_sphere) then
        do tt = 0, ns_glo
          do ii = 0, nx
            qfinal_glo(ii, tt) = qfinal_glo(ii, tt)*irr(ii)
          end do
        end do
      end if

      call contour_convolution(chainlen_glo, nx, ns_glo, coeff_ns_glo,  &
&                                qfinal_glo_aux, qfinal_glo, phi_glo)
    end if

    !edwards diffusion for grafted chains in the lower boundary
    if (exist_ghi) then

      !set the dirichlet boundary conditions
      n_dir_nodes = 0
      !dirichlet lower bound
      if (bc_lo_grafted .eq. F_bc_dirichlet_eq_0) then
        dir_nodes_id(n_dir_nodes) = 0
        dir_nodes_rdiag(n_dir_nodes) = 0.d0
        n_dir_nodes = n_dir_nodes + 1
      else if (bc_lo_grafted .eq. F_bc_dirichlet_eq_1) then
        dir_nodes_id(n_dir_nodes) = 0
        dir_nodes_rdiag(n_dir_nodes) = 1.0d0
        n_dir_nodes = n_dir_nodes + 1
      end if
      !dirichlet upper bound
      if (bc_hi_grafted .eq. F_bc_dirichlet_eq_0) then
        dir_nodes_id(n_dir_nodes) = nx
        dir_nodes_rdiag(n_dir_nodes) = 0.d0
        n_dir_nodes = n_dir_nodes + 1
      else if (bc_hi_grafted .eq. F_bc_dirichlet_eq_1) then
        dir_nodes_id(n_dir_nodes) = nx
        dir_nodes_rdiag(n_dir_nodes) = 1.0d0
        n_dir_nodes = n_dir_nodes + 1
      end if

      if (geometry .eq. F_sphere) then
        do ii = 0, n_dir_nodes - 1
          dir_nodes_rdiag(ii) = dir_nodes_rdiag(ii)*rr(dir_nodes_id(ii))
        end do
      end if

      !ghi_aux
      do ii = 0, nx
        qghi_aux(ii, 1) = 1.d0
        qfinal_ghi_aux(ii, 0) = 1.d0
      end do

      if (geometry .eq. F_sphere) then
        do ii = 0, nx
          qghi_aux(ii, 1) = qghi_aux(ii, 1)*rr(ii)
          qfinal_ghi_aux(ii, 0) = qfinal_ghi_aux(ii, 0)*rr(ii)
        end do
      end if

      call solver_edwards(bc_lo_grafted, bc_hi_grafted, n_dir_nodes, dir_nodes_id, dir_nodes_rdiag, &
&                           Rg2_per_mon_ghi, nx, ns_ghi, dx, ds_ghi, edwards_solver,        &
&                           linear_solver, wa_ifc_ghi, qghi_aux, qfinal_ghi_aux)

      if (geometry .eq. F_sphere) then
        do tt = 0, ns_ghi
          do ii = 0, nx
            qfinal_ghi_aux(ii, tt) = qfinal_ghi_aux(ii, tt)*irr(ii)
          end do
        end do
      end if

      ! grafted chains
      do ii = 0, nx
        qghi(ii, 1) = 0.d0
        qfinal_ghi(ii, 0) = 0.d0
      end do

      delta = (1.0/(dx(gnode_hi)*layer_area(gnode_hi)))*1e+30 !1/m3
      qinit_hi = chainlen_ghi*delta*(gdens_hi*surface_area) &
&                     /(rho_seg_bulk*qfinal_ghi_aux(gnode_hi, ns_ghi))

      qghi(gnode_hi, 1) = qinit_hi
      qfinal_ghi(gnode_hi, 0) = qinit_hi

      if (geometry .eq. F_sphere) then
        do ii = 0, nx
          qghi(ii, 1) = qghi(ii, 1)*rr(ii)
          qfinal_ghi(ii, 0) = qfinal_ghi(ii, 0)*rr(ii)
        end do
      end if

      call solver_edwards(bc_lo_grafted, bc_hi_grafted, n_dir_nodes, dir_nodes_id, dir_nodes_rdiag, &
&                           Rg2_per_mon_ghi, nx, ns_ghi, dx, ds_ghi, edwards_solver,        &
&                           linear_solver, wa_ifc_ghi, qghi, qfinal_ghi)

      if (geometry .eq. F_sphere) then
        do tt = 0, ns_ghi
          do ii = 0, nx
            qfinal_ghi(ii, tt) = qfinal_ghi(ii, tt)*irr(ii)
          end do
        end do
      end if

      call contour_convolution(chainlen_ghi, nx, ns_ghi, coeff_ns_ghi,  &
&                                qfinal_ghi_aux, qfinal_ghi, phi_ghi)
    end if

    phi_tot = 0.d0
    if (exist_mxa) phi_tot = phi_tot + phi_mxa
    if (exist_mxb) phi_tot = phi_tot + phi_mxb
    if (exist_glo) phi_tot = phi_tot + phi_glo
    if (exist_ghi) phi_tot = phi_tot + phi_ghi

    phi_kd1 = 0.d0
    phi_kd2 = 0.d0
    if (mxa_kind==1) phi_kd1 = phi_kd1 + phi_mxa
    if (mxb_kind==1) phi_kd1 = phi_kd1 + phi_mxb
    if (glo_kind==1) phi_kd1 = phi_kd1 + phi_glo
    if (ghi_kind==1) phi_kd1 = phi_kd1 + phi_ghi
    if (mxa_kind==2) phi_kd2 = phi_kd2 + phi_mxa
    if (mxb_kind==2) phi_kd2 = phi_kd2 + phi_mxb
    if (glo_kind==2) phi_kd2 = phi_kd2 + phi_glo
    if (ghi_kind==2) phi_kd2 = phi_kd2 + phi_ghi
    

    !calculate new fields
    wa_kd1 = 0.d0
    wa_kd2 = 0.d0

    !solid field
    do jj = 0, nx
      wa_kd1(jj) = wa_kd1(jj) + Ufield(jj) ! TODO: different solid field for kd1/2
      wa_kd2(jj) = wa_kd2(jj) + Ufield(jj)
    end do

    !mixing term
    if (dabs(chi12) .gt. 1e-7) then
      do jj = 0, nx
        wa_kd1(jj) = wa_kd1(jj) + chi12 * phi_kd2(jj)
        wa_kd2(jj) = wa_kd2(jj) + chi12 * phi_kd1(jj)
      end do
    end if

    !eos term
    
      do jj = 0, nx
        wa_kd1(jj) = wa_kd1(jj) + eos_df_drho(phi_tot(jj))*beta
        wa_kd2(jj) = wa_kd2(jj) + eos_df_drho(phi_tot(jj))*beta
      end do
    


    ! gradient term
    if (square_gradient) then
      do jj = 1, nx - 1
        d2phi_dr2(jj) = (phi_tot(jj - 1) - 2.d0*phi_tot(jj) + phi_tot(jj + 1))/(dx(jj)*1.e-10)**2
        dphi_dr(jj) = (phi_tot(jj + 1) - phi_tot(jj))/(dx(jj)*1.e-10)
      end do
      d2phi_dr2(0) = (phi_tot(0) - 2.d0*phi_tot(0) + phi_tot(1))/(dx(1)*1.e-10)**2
      d2phi_dr2(nx) = (phi_tot(nx - 1) - 2.d0*phi_tot(nx) + phi_tot(nx))/(dx(nx)*1.e-10)**2

      dphi_dr(0) = (phi_tot(1) - phi_tot(0))/(dx(1)*1.e-10)
      dphi_dr(nx) = 0.d0

      do jj = 0, nx
        wa_kd1(jj) = wa_kd1(jj) - k_gr*(rho_seg_bulk*d2phi_dr2(jj))*beta
        wa_kd2(jj) = wa_kd2(jj) - k_gr*(rho_seg_bulk*d2phi_dr2(jj))*beta
      end do
    end if

    !subtract bulk contribution
    wa_bulk_kd1 = eos_df_drho(1.d0)*beta
    wa_bulk_kd2 = eos_df_drho(1.d0)*beta
    wa_ifc_new_kd1 = wa_kd1 - wa_bulk_kd1
    wa_ifc_new_kd2 = wa_kd2 - wa_bulk_kd2
    

    

    !apply field mixing rule and update field
   
    do jj = 0, nx
      wa_ifc_kd1(jj) = (1.d0 - frac)*wa_ifc_kd1(jj) + frac*wa_ifc_new_kd1(jj)
      wa_ifc_kd2(jj) = (1.d0 - frac)*wa_ifc_kd2(jj) + frac*wa_ifc_new_kd2(jj)
    end do
    if (eos_type.eq.F_incompressible) then
      do jj = 0, nx
        pressure(jj)= chi12*(phi_kd1(jj) + phi_kd2(jj) - 1)
        wa_ifc_kd1(jj) = wa_ifc_kd1(jj) + 0.5*pressure(jj)
        wa_ifc_kd2(jj) = wa_ifc_kd2(jj) + 0.5*pressure(jj)
      end do
    end if
   
   
   
    if (mxa_kind==1) wa_ifc_mxa = phi_glo*chi12+wa_ifc_kd1
    if (mxb_kind==1) wa_ifc_mxb = wa_ifc_kd1
    if (glo_kind==1) wa_ifc_glo = wa_ifc_kd1
    if (ghi_kind==1) wa_ifc_ghi = wa_ifc_kd1
    if (mxa_kind==2) wa_ifc_mxa = wa_ifc_kd2
    if (mxb_kind==2) wa_ifc_mxb = wa_ifc_kd2
    if (glo_kind==2) wa_ifc_glo = phi_mxa*chi12+wa_ifc_kd2
    if (ghi_kind==2) wa_ifc_ghi = wa_ifc_kd2
    
    
    
   
   
    do jj = 0, nx
      wa_error_new = max(wa_error_new,dabs(pressure(jj)))
    end do

    !The present section checks the behavior of the field.
    !In case it diverges or it converges to a bad solution, the field is
    !restored and the fraction is decreased
    if (check_stability_every .gt. 0) then
    if (mod(iter, check_stability_every) .eq. 0 .and. iter .gt. 0) then
      restart = .false.
      restore = .false.

      if (wa_error_new > wa_error_old) then
        write (iow, '(3X,"*wa_error_new > wa_error_old:",E16.7,">",E16.7)') wa_error_new, wa_error_old
        write (6, '(3X,"*wa_error_new > wa_error_old:",E16.7,">",E16.7)') wa_error_new, wa_error_old
        restore = .true.
      end if
      if (isnan(wa_error_new)) then
        write (iow, '(3X,"*current field is not a number!")')
        write (6, '(3X,"*current field is not a number!")')
        restart = .true.
      end if

      if (exist_mxa .or. exist_mxb) then
        do jj = 1, nx - 1
          if (phi_tot(jj) .lt. 1.e-10) then
            restart = .true.
            write (iow, '(3X,"*Convergence to unphysical solution!")')
            write (6, '(3X,"*Convergence to unphysical solution!")')
            exit
          end if
        end do
      end if

      if (restart) then
        wa_ifc_kd1 = 0.d0
        wa_ifc_backup_kd1 = 0.d0
        frac = frac*0.5d0
        wa_error_old = 1.d10
        write (iow, '(3X,"Restart simulation with a new fraction:",E16.7)') frac
        write (6, '(3X,"Restart simulation with a new fraction:",E16.7)') frac
      elseif (restore) then
        frac = frac*0.5d0
        wa_ifc_kd1 = wa_ifc_backup_kd1
        wa_error_old = 1.d10
        write (iow, '(3X,"Restore previous field with a new fraction:",E16.7)') frac
        write (6, '(3X,"Restore previous field with a new fraction:",E16.7)') frac
      else
        wa_error_old = wa_error_new
        wa_ifc_backup_kd1 = wa_ifc_kd1
      end if
    end if
    end if

    convergence = (wa_error_new .lt. 1e-10)

    if (mod(iter, field_every) .eq. 0 .or. convergence) call export_field_binary(wa_ifc_kd1, nx)
    if (compute_every .gt. 0) then
      if (mod(iter, compute_every) .eq. 0 .or. convergence) call export_computes(qinit_lo, qinit_hi)
    end if
    if (mod(iter, thermo_every) .eq. 0 .or. convergence) then
      if (exist_mxa) nch_mxa = get_nchains(coeff_nx, nx, layer_area, phi_mxa, &
&                                                      rho_seg_bulk, chainlen_mxa)
      if (exist_mxb) nch_mxb = get_nchains(coeff_nx, nx, layer_area, phi_mxb, &
&                                                      rho_seg_bulk, chainlen_mxb)
      if (exist_glo) nchglo = get_nchains(coeff_nx, nx, layer_area, phi_glo,  &
&                                                      rho_seg_bulk, chainlen_glo)
      if (exist_ghi) nchghi = get_nchains(coeff_nx, nx, layer_area, phi_ghi,  &
&                                                      rho_seg_bulk, chainlen_ghi)

      call compute_energies(free_energy)

      ! flush the log file
      close (iow)
      open (unit=iow, file="o.log", position='append')

      write (iow, '(2X,I15,10(2X,E15.7))') iter, free_energy, wa_error_new, nchglo/surface_area, &
&                                         nchghi/surface_area, nch_mxa, nch_mxb, nchglo, nchghi, frac
      write (*, '(2X,I15,10(2X,E15.7))') iter, free_energy, wa_error_new, nchglo/surface_area, &
&                                         nchghi/surface_area, nch_mxa, nch_mxb, nchglo, nchghi, frac
    end if
    if (convergence) exit
  end do

  write (iow, '(A85)')
  write (*, '(A85)')

  if (iter .eq. max_iter) write (iow, '("Convergence of max iterations",I16," iterations")') iter
  if (wa_error_new .lt. max_wa_error) write (iow, '("Convergence of max field error",E16.7," k_B T")') wa_error_new
  if (iter .eq. max_iter) write (*, '("Convergence of max iterations",I16," iterations")') iter
  if (wa_error_new .lt. max_wa_error) write (*, '("Convergence of max field error",E16.7," k_B T")') wa_error_new

  write (iow, '(A85)') adjl('-------------------------------------FINAL OUTPUT------------------------------------', 85)
  write (*, '(A85)') adjl('-------------------------------------FINAL OUTPUT------------------------------------', 85)
  write (iow, '(3X,A17,E16.7," mN/m")') adjl("Free energy:", 17), free_energy
  write (*, '(3X,A17,E16.7," mN/m")') adjl("Free energy:", 17), free_energy
  if (exist_mxa) then
    write (iow, '(A85)') adjl('------------------------------------mxa CHAINS-----------------------------------', 85)
    write (*, '(A85)') adjl('------------------------------------mxa CHAINS-----------------------------------', 85)
    write (iow, '(3X,A17,E16.7," chains")') adjl("matrix chains:", 17), nch_mxa
    write (*, '(3X,A17,E16.7," chains")') adjl("matrix chains:", 17), nch_mxa
  end if
  if (exist_mxb) then
    write (iow, '(A85)') adjl('------------------------------------mxb CHAINS-----------------------------------', 85)
    write (*, '(A85)') adjl('------------------------------------mxb CHAINS-----------------------------------', 85)
    write (iow, '(3X,A17,E16.7," chains")') adjl("matrix chains:", 17), nch_mxb
    write (*, '(3X,A17,E16.7," chains")') adjl("matrix chains:", 17), nch_mxb
  end if
  if (exist_glo) then
    write (iow, '(A85)') adjl('---------------------------------glo CHAINS-----------------------------------', 85)
    write (*, '(A85)') adjl('---------------------------------glo CHAINS-----------------------------------', 85)
    write (iow, '(3X,A17,E16.7," chains")') adjl("glo chains:", 17), nchglo
    write (*, '(3X,A17,E16.7," chains")') adjl("glo chains:", 17), nchglo

    write (iow, '(3X,A17,E16.7)') adjl("q_f(r_gi,0):", 17), qfinal_glo_aux(gnode_lo, ns_glo)
    write (*, '(3X,A17,E16.7)') adjl("q_f(r_gi,0):", 17), qfinal_glo_aux(gnode_lo, ns_glo)
    write (iow, '(3X,A17,E16.7)') adjl("q_g(r_gi,N_g):", 17), qfinal_glo(gnode_lo, 0)
    write (*, '(3X,A17,E16.7)') adjl("q_g(r_gi,N_g):", 17), qfinal_glo(gnode_lo, 0)
    write (iow, '(3X,A17,E16.7," chains/angstrom")') adjl("grafting density:", 17), nchglo/surface_area
    write (*, '(3X,A17,E16.7," chains/angstrom")') adjl("grafting density:", 17), nchglo/surface_area
    write (iow, '(3X,A17,E16.7)') adjl("error:", 17), (nchglo/surface_area - gdens_lo)/gdens_lo*100.0
    write (*, '(3X,A17,E16.7)') adjl("error:", 17), (nchglo/surface_area - gdens_lo)/gdens_lo*100.0
  end if
  if (exist_ghi) then
    write (iow, '(A85)') adjl("---------------------------------ghi CHAINS-----------------------------------", 85)
    write (*, '(A85)') adjl("---------------------------------ghi CHAINS-----------------------------------", 85)
    write (iow, '(3X,A17,E16.7," chains")') adjl("ghi chains:", 17), nchghi
    write (*, '(3X,A17,E16.7," chains")') adjl("ghi chains:", 17), nchghi

    write (iow, '(3X,A17,E16.7)') adjl("q_f(r_gi,0):", 17), qfinal_ghi_aux(gnode_hi, ns_ghi)
    write (*, '(3X,A17,E16.7)') adjl("q_f(r_gi,0):", 17), qfinal_ghi_aux(gnode_hi, ns_ghi)
    write (iow, '(3X,A17,E16.7)') adjl("q_g(r_gi,N_g):", 17), qfinal_ghi(gnode_hi, 0)
    write (*, '(3X,A17,E16.7)') adjl("q_g(r_gi,N_g):", 17), qfinal_ghi(gnode_hi, 0)
    write (iow, '(3X,A17,E16.7," chains/angstrom")') adjl("grafting density:", 17), nchghi/surface_area
    write (*, '(3X,A17,E16.7," chains/angstrom")') adjl("grafting density:", 17), nchghi/surface_area
    write (iow, '(3X,A17,E16.7)') adjl("error:", 17), (nchghi/surface_area - gdens_hi)/gdens_hi*100.0
    write (*, '(3X,A17,E16.7)') adjl("error:", 17), (nchghi/surface_area - gdens_hi)/gdens_hi*100.0
  end if
!----------------------------------------------------------------------------------------------------------!
end program fd_1d
# diplo
