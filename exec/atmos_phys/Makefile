# Makefile created by mkmf 19.3.0

CPPDEFS = -DGFDL_COMPATIBLE_MICROP -DINTERNAL_FILE_NML

OTHERFLAGS = -I$(BUILDROOT)fms

include $(MK_TEMPLATE)


.DEFAULT:
	-echo $@ does not exist.
all: libatmos_phys.a
CLUBB_3D_var.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_3D_var.F90 parameters_model.o constants_clubb.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_3D_var.F90
CLUBB_driver_SCM.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_driver_SCM.F90 aerosol_types.o aer_ccn_act.o aer_ccn_act_k.o ice_nucl.o alt_cloud.o clubb_core.o clubb_precision.o constants_clubb.o error_code.o grid_class.o parameter_indices.o parameters_model.o parameters_tunable.o stats_subs.o stats_variables.o T_in_K_module.o variables_prognostic_module.o polysvp.o aerosol_params.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_driver_SCM.F90
MG_microp_3D.o: $(SRCROOT)atmos_param_am3/clubb/MG_microp_3D.F90 cloud_rad.o cldwat2m_micro.o micro_mg.o strat_cloud_utilities.o strat_netcdf.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/MG_microp_3D.F90
MISR_simulator.o: $(SRCROOT)atmos_param_am3/cosp/MISR_simulator/MISR_simulator.F $(SRCROOT)atmos_param_am3/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/cosp	$(SRCROOT)atmos_param_am3/cosp/MISR_simulator/MISR_simulator.F
Skw_module.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/Skw_module.F90 constants_clubb.o clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/Skw_module.F90
T_in_K_module.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/T_in_K_module.F90 constants_clubb.o clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/T_in_K_module.F90
advance_helper_module.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/advance_helper_module.F90 clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/advance_helper_module.F90
advance_sclrm_Nd_module.o: $(SRCROOT)atmos_param_am3/clubb/advance_sclrm_Nd_module.F90 clubb_precision.o grid_class.o parameters_model.o lapack_wrap.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/advance_sclrm_Nd_module.F90
advance_windm_edsclrm_module.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/advance_windm_edsclrm_module.F90 grid_class.o parameters_model.o parameters_tunable.o model_flags.o clubb_precision.o stats_type.o stats_variables.o clip_explicit.o error_code.o constants_clubb.o sponge_layer_damping.o lapack_wrap.o diffusion.o mean_adv.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/advance_windm_edsclrm_module.F90
advance_wp2_wp3_module.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/advance_wp2_wp3_module.F90 grid_class.o parameters_tunable.o stats_type.o stats_variables.o constants_clubb.o model_flags.o clubb_precision.o error_code.o lapack_wrap.o fill_holes.o clip_explicit.o csr_matrix_class_3array.o gmres_wrap.o gmres_cache.o diffusion.o mean_adv.o hyper_diffusion_4th_ord.o advance_helper_module.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/advance_wp2_wp3_module.F90
advance_xm_wpxp_module.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/advance_xm_wpxp_module.F90 parameters_tunable.o constants_clubb.o parameters_model.o grid_class.o model_flags.o mono_flux_limiter.o pdf_parameter_module.o clubb_precision.o error_code.o stats_type.o stats_variables.o sponge_layer_damping.o diffusion.o mean_adv.o clip_semi_implicit.o advance_helper_module.o lapack_wrap.o pos_definite_module.o clip_explicit.o fill_holes.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/advance_xm_wpxp_module.F90
advance_xp2_xpyp_module.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/advance_xp2_xpyp_module.F90 constants_clubb.o model_flags.o parameters_tunable.o parameters_model.o grid_class.o clubb_precision.o clip_explicit.o stats_type.o error_code.o stats_variables.o array_index.o diffusion.o mean_adv.o advance_helper_module.o lapack_wrap.o interpolation.o fill_holes.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/advance_xp2_xpyp_module.F90
aer_ccn_act.o: $(SRCROOT)atmos_shared_am3/tracer_driver/aer_ccn_act/aer_ccn_act.F90 aer_ccn_act_k.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/aer_ccn_act/aer_ccn_act.F90
aer_ccn_act_k.o: $(SRCROOT)atmos_shared_am3/tracer_driver/aer_ccn_act/aer_ccn_act_k.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/aer_ccn_act/aer_ccn_act_k.F90
aer_in_act.o: $(SRCROOT)atmos_shared_am3/tracer_driver/aer_ccn_act/aer_in_act.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/aer_ccn_act/aer_in_act.F90
aerosol.o: $(SRCROOT)atmos_shared_am3/aerosol/aerosol.F90 aerosol_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/aerosol/aerosol.F90
aerosol_cloud.o: $(SRCROOT)atmos_param_am3/strat_cloud/aerosol_cloud.F90 aerosol_types.o polysvp.o aerosol_params.o aer_ccn_act.o ice_nucl.o strat_cloud_utilities.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/strat_cloud/aerosol_cloud.F90
aerosol_params.o: $(SRCROOT)atmos_shared_am3/tracer_driver/aer_ccn_act/aerosol_params.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/aer_ccn_act/aerosol_params.F90
aerosol_types.o: $(SRCROOT)atmos_shared_am3/aerosol/aerosol_types.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/aerosol/aerosol_types.F90
aerosolrad_driver.o: $(SRCROOT)atmos_param_am3/radiation/aerosols/aerosolrad_driver.F90 aerosol_types.o aerosol.o aerosolrad_types.o aerosolrad_package.o sea_esf_rad.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/radiation/aerosols/aerosolrad_driver.F90
aerosolrad_package.o: $(SRCROOT)atmos_param_am3/radiation/aerosols/aerosolrad_package.F90 aerosol_types.o aerosolrad_types.o sealw99.o esfsw_driver.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/radiation/aerosols/aerosolrad_package.F90
aerosolrad_types.o: $(SRCROOT)atmos_param_am3/radiation/aerosols/aerosolrad_types.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/radiation/aerosols/aerosolrad_types.F90
alt_cloud.o: $(SRCROOT)atmos_param_am3/clubb/alt_cloud.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/alt_cloud.F90
anl_erf.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/anl_erf.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/anl_erf.F90
array_index.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/array_index.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/array_index.F90
array_lib.o: $(SRCROOT)atmos_param_am3/cosp/quickbeam/array_lib.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_defs.H mrgrnk.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/cosp	$(SRCROOT)atmos_param_am3/cosp/quickbeam/array_lib.F90
atmos_age_tracer.o: $(SRCROOT)atmos_shared_am3/tracer_driver/atmos_age_tracer.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/atmos_age_tracer.F90
atmos_carbon_aerosol.o: $(SRCROOT)atmos_shared_am3/tracer_driver/atmos_carbon_aerosol.F90 atmos_tracer_utilities.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/atmos_carbon_aerosol.F90
atmos_ch3i.o: $(SRCROOT)atmos_shared_am3/tracer_driver/atmos_ch3i.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/atmos_ch3i.F90
atmos_cmip_diag.o: $(SRCROOT)atmos_shared_am3/atmos_cmip_diag/atmos_cmip_diag.F90 $(SRCROOT)atmos_shared_am3/atmos_cmip_diag/atmos_cmip_interp.inc
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_shared_am3/atmos_cmip_diag	$(SRCROOT)atmos_shared_am3/atmos_cmip_diag/atmos_cmip_diag.F90
atmos_co2.o: $(SRCROOT)atmos_shared_am3/tracer_driver/atmos_co2.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/atmos_co2.F90
atmos_convection_tracer.o: $(SRCROOT)atmos_shared_am3/tracer_driver/atmos_convection_tracer.F90 atmos_tracer_utilities.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/atmos_convection_tracer.F90
atmos_dust.o: $(SRCROOT)atmos_shared_am3/tracer_driver/atmos_dust.F90 atmos_tracer_utilities.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/atmos_dust.F90
atmos_lib.o: $(SRCROOT)atmos_param_am3/cosp/quickbeam/atmos_lib.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/cosp	$(SRCROOT)atmos_param_am3/cosp/quickbeam/atmos_lib.F90
atmos_nudge.o: $(SRCROOT)atmos_shared_am3/atmos_nudge/atmos_nudge.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/atmos_nudge/atmos_nudge.F90
atmos_radon.o: $(SRCROOT)atmos_shared_am3/tracer_driver/atmos_radon.F90 atmos_tracer_utilities.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/atmos_radon.F90
atmos_regional_tracer_driver.o: $(SRCROOT)atmos_shared_am3/tracer_driver/atmos_regional_tracer_driver.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/atmos_regional_tracer_driver.F90
atmos_sea_salt.o: $(SRCROOT)atmos_shared_am3/tracer_driver/atmos_sea_salt.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/atmos_sea_salt.F90
atmos_soa.o: $(SRCROOT)atmos_shared_am3/tracer_driver/atmos_soa.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/atmos_soa.F90
atmos_sulfate.o: $(SRCROOT)atmos_shared_am3/tracer_driver/atmos_sulfate.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/atmos_sulfate.F90
atmos_sulfur_hex.o: $(SRCROOT)atmos_shared_am3/tracer_driver/atmos_sulfur_hex.F90 atmos_tracer_utilities.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/atmos_sulfur_hex.F90
atmos_tracer_driver.o: $(SRCROOT)atmos_shared_am3/tracer_driver/atmos_tracer_driver.F90 atmos_tracer_utilities.o atmos_radon.o atmos_carbon_aerosol.o atmos_convection_tracer.o atmos_sulfur_hex.o atmos_ch3i.o atmos_sea_salt.o atmos_dust.o atmos_sulfate.o atmos_soa.o tropchem_driver.o atmos_regional_tracer_driver.o strat_chem_driver.o atmos_age_tracer.o atmos_co2.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/atmos_tracer_driver.F90
atmos_tracer_utilities.o: $(SRCROOT)atmos_shared_am3/tracer_driver/atmos_tracer_utilities.F90 monin_obukhov.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/atmos_tracer_utilities.F90
betaDistribution.o: $(SRCROOT)atmos_param_am3/cloud_generator/betaDistribution.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/cloud_generator/betaDistribution.F90
betts_miller.o: $(SRCROOT)atmos_param_am3/betts_miller/betts_miller.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/betts_miller/betts_miller.F90
bm_massflux.o: $(SRCROOT)atmos_param_am3/betts_miller/bm_massflux.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/betts_miller/bm_massflux.F90
bm_omp.o: $(SRCROOT)atmos_param_am3/betts_miller/bm_omp.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/betts_miller/bm_omp.F90
bulkphys_rad.o: $(SRCROOT)atmos_param_am3/radiation/clouds/bulkphys_rad.F90 cloudrad_types.o strat_clouds_W.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/radiation/clouds/bulkphys_rad.F90
calendar.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/calendar.F90 clubb_precision.o constants_clubb.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/calendar.F90
cg_drag.o: $(SRCROOT)atmos_param_am3/cg_drag/cg_drag.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/cg_drag/cg_drag.F90
check_nan.o: $(SRCROOT)atmos_param_am3/strat_cloud/check_nan.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/strat_cloud/check_nan.F90
cldwat2m_micro.o: $(SRCROOT)atmos_param_am3/strat_cloud/cldwat2m_micro.F90 gamma_mg.o strat_cloud_utilities.o simple_pdf.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/strat_cloud/cldwat2m_micro.F90
clip_explicit.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/clip_explicit.F90 grid_class.o parameters_model.o model_flags.o clubb_precision.o stats_type.o stats_variables.o constants_clubb.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/clip_explicit.F90
clip_semi_implicit.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/clip_semi_implicit.F90 clubb_precision.o constants_clubb.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/clip_semi_implicit.F90
cloud_generator.o: $(SRCROOT)atmos_param_am3/cloud_generator/cloud_generator.F90 betaDistribution.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/cloud_generator/cloud_generator.F90
cloud_obs.o: $(SRCROOT)atmos_param_am3/cloud_obs/cloud_obs.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/cloud_obs/cloud_obs.F90
cloud_rad.o: $(SRCROOT)atmos_param_am3/cloud_rad/cloud_rad.F90 gamma_mg.o mg_const.o aerosol_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/cloud_rad/cloud_rad.F90
cloud_spec.o: $(SRCROOT)atmos_param_am3/radiation/clouds/cloud_spec.F90 get_random_number_stream.o aerosol_types.o physics_radiation_exch.o sea_esf_rad.o cloudrad_types.o strat_clouds_W.o donner_deep_clouds_W.o uw_clouds_W.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/radiation/clouds/cloud_spec.F90
cloud_zonal.o: $(SRCROOT)atmos_param_am3/cloud_zonal/cloud_zonal.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/cloud_zonal/cloud_zonal.F90
cloudrad_diagnostics.o: $(SRCROOT)atmos_param_am3/radiation/clouds/cloudrad_diagnostics.F90 cloudrad_types.o microphys_rad.o isccp_clouds.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/radiation/clouds/cloudrad_diagnostics.F90
cloudrad_driver.o: $(SRCROOT)atmos_param_am3/radiation/clouds/cloudrad_driver.F90 aerosol_types.o physics_radiation_exch.o sea_esf_rad.o cloudrad_types.o cloud_spec.o cloudrad_package.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/radiation/clouds/cloudrad_driver.F90
cloudrad_package.o: $(SRCROOT)atmos_param_am3/radiation/clouds/cloudrad_package.F90 sea_esf_rad.o cloudrad_types.o cloudrad_diagnostics.o bulkphys_rad.o microphys_rad.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/radiation/clouds/cloudrad_package.F90
cloudrad_types.o: $(SRCROOT)atmos_param_am3/radiation/clouds/cloudrad_types.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/radiation/clouds/cloudrad_types.F90
clouds.o: $(SRCROOT)atmos_param_am3/clouds/clouds.F90 cloud_rad.o cloud_zonal.o cloud_obs.o rh_clouds.o strat_cloud.o diag_cloud.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clouds/clouds.F90
clubb_core.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/clubb_core.F90 constants_clubb.o parameters_tunable.o parameters_model.o model_flags.o grid_class.o numerical_check.o variables_diagnostic_module.o pdf_parameter_module.o advance_sclrm_Nd_module.o advance_xm_wpxp_module.o advance_xp2_xpyp_module.o surface_varnce_module.o pdf_closure_module.o mixing_length.o advance_windm_edsclrm_module.o saturation.o advance_wp2_wp3_module.o clubb_precision.o error_code.o Skw_module.o clip_explicit.o T_in_K_module.o stats_subs.o stats_type.o stats_variables.o fill_holes.o sigma_sqd_w_module.o parameter_indices.o variables_prognostic_module.o csr_matrix_class_3array.o gmres_wrap.o gmres_cache.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/clubb_core.F90
clubb_precision.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/clubb_precision.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/clubb_precision.F90
constants_clubb.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/constants_clubb.F90 clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/constants_clubb.F90
conv_closures.o: $(SRCROOT)atmos_param_am3/shallow_cu/conv_closures.F90 conv_utilities_k.o conv_plumes_k.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/shallow_cu/conv_closures.F90
conv_plumes.o: $(SRCROOT)atmos_param_am3/shallow_cu/conv_plumes.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/shallow_cu/conv_plumes.F90
conv_plumes_k.o: $(SRCROOT)atmos_param_am3/shallow_cu/conv_plumes_k.F90 $(SRCROOT)fms/include/fms_platform.h aer_ccn_act_k.o conv_utilities_k.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)atmos_param_am3/shallow_cu/conv_plumes_k.F90
conv_utilities.o: $(SRCROOT)atmos_param_am3/shallow_cu/conv_utilities.F90 conv_utilities_k.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/shallow_cu/conv_utilities.F90
conv_utilities_k.o: $(SRCROOT)atmos_param_am3/shallow_cu/conv_utilities_k.F90 $(SRCROOT)fms/include/fms_platform.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)atmos_param_am3/shallow_cu/conv_utilities_k.F90
corr_matrix_module.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/corr_matrix_module.F90 input_reader.o matrix_operations.o constants_clubb.o clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/corr_matrix_module.F90
cosp.o: $(SRCROOT)atmos_param_am3/cosp/cosp.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_defs.H cosp_types.o cosp_simulator.o cosp_modis_simulator.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/cosp	$(SRCROOT)atmos_param_am3/cosp/cosp.F90
cosp_constants.o: $(SRCROOT)atmos_param_am3/cosp/cosp_constants.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/cosp	$(SRCROOT)atmos_param_am3/cosp/cosp_constants.F90
cosp_diagnostics.o: $(SRCROOT)atmos_param_am3/cosp/cosp_diagnostics.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_defs.H cosp_types.o cosp_io.o cosp_constants.o lmd_ipsl_stats.o cosp_modis_simulator.o modis_simulator.o cosp_utils.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/cosp	$(SRCROOT)atmos_param_am3/cosp/cosp_diagnostics.F90
cosp_driver.o: $(SRCROOT)atmos_param_am3/cosp/cosp_driver.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_defs.H cosp_types.o cosp.o cosp_io.o cosp_constants.o cosp_modis_simulator.o cosp_diagnostics.o cosp_utils.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/cosp	$(SRCROOT)atmos_param_am3/cosp/cosp_driver.F90
cosp_io.o: $(SRCROOT)atmos_param_am3/cosp/cosp_io.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_defs.H cosp_constants.o cosp_types.o cosp_modis_simulator.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/cosp	$(SRCROOT)atmos_param_am3/cosp/cosp_io.F90
cosp_isccp_simulator.o: $(SRCROOT)atmos_param_am3/cosp/cosp_isccp_simulator.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_defs.H cosp_constants.o cosp_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/cosp	$(SRCROOT)atmos_param_am3/cosp/cosp_isccp_simulator.F90
cosp_lidar.o: $(SRCROOT)atmos_param_am3/cosp/cosp_lidar.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_defs.H cosp_constants.o cosp_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/cosp	$(SRCROOT)atmos_param_am3/cosp/cosp_lidar.F90
cosp_misr_simulator.o: $(SRCROOT)atmos_param_am3/cosp/cosp_misr_simulator.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_defs.H cosp_constants.o cosp_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/cosp	$(SRCROOT)atmos_param_am3/cosp/cosp_misr_simulator.F90
cosp_modis_simulator.o: $(SRCROOT)atmos_param_am3/cosp/cosp_modis_simulator.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_defs.H cosp_constants.o cosp_types.o modis_simulator.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/cosp	$(SRCROOT)atmos_param_am3/cosp/cosp_modis_simulator.F90
cosp_radar.o: $(SRCROOT)atmos_param_am3/cosp/llnl/cosp_radar.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_defs.H cosp_constants.o cosp_types.o cosp_utils.o radar_simulator_types.o array_lib.o atmos_lib.o format_input.o mrgrnk.o math_lib.o optics_lib.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/cosp	$(SRCROOT)atmos_param_am3/cosp/llnl/cosp_radar.F90
cosp_rttov_simulator.o: $(SRCROOT)atmos_param_am3/cosp/cosp_rttov_simulator.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_defs.H cosp_constants.o cosp_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/cosp	$(SRCROOT)atmos_param_am3/cosp/cosp_rttov_simulator.F90
cosp_simulator.o: $(SRCROOT)atmos_param_am3/cosp/cosp_simulator.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_defs.H cosp_types.o cosp_radar.o cosp_lidar.o cosp_isccp_simulator.o cosp_modis_simulator.o cosp_misr_simulator.o cosp_rttov_simulator.o cosp_stats.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/cosp	$(SRCROOT)atmos_param_am3/cosp/cosp_simulator.F90
cosp_stats.o: $(SRCROOT)atmos_param_am3/cosp/cosp_stats.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_defs.H cosp_constants.o cosp_types.o llnl_stats.o lmd_ipsl_stats.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/cosp	$(SRCROOT)atmos_param_am3/cosp/cosp_stats.F90
cosp_types.o: $(SRCROOT)atmos_param_am3/cosp/cosp_types.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_defs.H cosp_constants.o cosp_utils.o radar_simulator_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/cosp	$(SRCROOT)atmos_param_am3/cosp/cosp_types.F90
cosp_utils.o: $(SRCROOT)atmos_param_am3/cosp/cosp_utils.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_defs.H cosp_constants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/cosp	$(SRCROOT)atmos_param_am3/cosp/cosp_utils.F90
csr_matrix_class_3array.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/csr_matrix_class_3array.F90 constants_clubb.o grid_class.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/csr_matrix_class_3array.F90
cu_mo_trans.o: $(SRCROOT)atmos_param_am3/cu_mo_trans/cu_mo_trans.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/cu_mo_trans/cu_mo_trans.F90
cumulus_closure_k.o: $(SRCROOT)atmos_param_am3/donner_deep/cumulus_closure_k.F90 donner_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/donner_deep/cumulus_closure_k.F90
damping_driver.o: $(SRCROOT)atmos_param_am3/damping_driver/damping_driver.F90 mg_drag.o cg_drag.o topo_drag.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/damping_driver/damping_driver.F90
deep_conv.o: $(SRCROOT)atmos_param_am3/shallow_cu/deep_conv.F90 conv_utilities.o conv_utilities_k.o conv_plumes_k.o conv_closures.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/shallow_cu/deep_conv.F90
detr_ice_num.o: $(SRCROOT)atmos_param_am3/moist_processes/detr_ice_num.F90 mg_const.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/moist_processes/detr_ice_num.F90
diag_cloud.o: $(SRCROOT)atmos_param_am3/diag_cloud/diag_cloud.F90 cloud_zonal.o diag_cloud_rad.o shallow_conv.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/diag_cloud/diag_cloud.F90
diag_cloud_rad.o: $(SRCROOT)atmos_param_am3/diag_cloud_rad/diag_cloud_rad.F90 cloud_rad.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/diag_cloud_rad/diag_cloud_rad.F90
diag_integral.o: $(SRCROOT)atmos_param_am3/diag_integral/diag_integral.F90 $(SRCROOT)fms/include/fms_platform.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)atmos_param_am3/diag_integral/diag_integral.F90
diffusion.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/diffusion.F90 grid_class.o clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/diffusion.F90
diffusivity.o: $(SRCROOT)atmos_param_am3/diffusivity/diffusivity.F90 monin_obukhov.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/diffusivity/diffusivity.F90
donner_cape_k.o: $(SRCROOT)atmos_param_am3/donner_deep/donner_cape_k.F90 donner_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/donner_deep/donner_cape_k.F90
donner_cloud_model_k.o: $(SRCROOT)atmos_param_am3/donner_deep/donner_cloud_model_k.F90 donner_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/donner_deep/donner_cloud_model_k.F90
donner_deep.o: $(SRCROOT)atmos_param_am3/donner_deep/donner_deep.F90 donner_types.o conv_utilities_k.o conv_plumes_k.o fms_donner.o nonfms_donner.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/donner_deep/donner_deep.F90
donner_deep_clouds_W.o: $(SRCROOT)atmos_param_am3/radiation/clouds/donner_deep_clouds_W.F90 cloudrad_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/radiation/clouds/donner_deep_clouds_W.F90
donner_deep_k.o: $(SRCROOT)atmos_param_am3/donner_deep/donner_deep_k.F90 donner_types.o conv_utilities_k.o conv_plumes_k.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/donner_deep/donner_deep_k.F90
donner_deep_miz.o: $(SRCROOT)atmos_param_am3/donner_deep/donner_deep_miz.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/donner_deep/donner_deep_miz.F90
donner_lite_k.o: $(SRCROOT)atmos_param_am3/donner_deep/donner_lite_k.F90 donner_types.o conv_utilities_k.o conv_plumes_k.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/donner_deep/donner_lite_k.F90
donner_lscloud_k.o: $(SRCROOT)atmos_param_am3/donner_deep/donner_lscloud_k.F90 donner_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/donner_deep/donner_lscloud_k.F90
donner_meso_k.o: $(SRCROOT)atmos_param_am3/donner_deep/donner_meso_k.F90 donner_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/donner_deep/donner_meso_k.F90
donner_rad_k.o: $(SRCROOT)atmos_param_am3/donner_deep/donner_rad_k.F90 donner_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/donner_deep/donner_rad_k.F90
donner_types.o: $(SRCROOT)atmos_param_am3/donner_deep/donner_types.F90 $(SRCROOT)atmos_param_am3/donner_deep/donner_types.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/donner_deep	$(SRCROOT)atmos_param_am3/donner_deep/donner_types.F90
donner_utilities_k.o: $(SRCROOT)atmos_param_am3/donner_deep/donner_utilities_k.F90 donner_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/donner_deep/donner_utilities_k.F90
dry_adj.o: $(SRCROOT)atmos_param_am3/dry_adj/dry_adj.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/dry_adj/dry_adj.F90
dsd.o: $(SRCROOT)atmos_param_am3/cosp/quickbeam/dsd.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_defs.H array_lib.o math_lib.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/cosp	$(SRCROOT)atmos_param_am3/cosp/quickbeam/dsd.F90
edt.o: $(SRCROOT)atmos_param_am3/edt/edt.F90 monin_obukhov.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/edt/edt.F90
endian.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/endian.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/endian.F90
entrain.o: $(SRCROOT)atmos_param_am3/entrain/entrain.F90 monin_obukhov.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/entrain/entrain.F90
error_code.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/error_code.F90 constants_clubb.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/error_code.F90
esfsw_bands.o: $(SRCROOT)atmos_param_am3/radiation/shortwave/esfsw_bands.F90 esfsw_utilities.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/radiation/shortwave/esfsw_bands.F90
esfsw_driver.o: $(SRCROOT)atmos_param_am3/radiation/shortwave/esfsw_driver.F90 esfsw_parameters.o esfsw_bands.o esfsw_utilities.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/radiation/shortwave/esfsw_driver.F90
esfsw_parameters.o: $(SRCROOT)atmos_param_am3/radiation/shortwave/esfsw_parameters.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/radiation/shortwave/esfsw_parameters.F90
esfsw_utilities.o: $(SRCROOT)atmos_param_am3/radiation/shortwave/esfsw_utilities.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/radiation/shortwave/esfsw_utilities.F90
extrapolation.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/extrapolation.F90 clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/extrapolation.F90
file_functions.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/file_functions.F90 clubb_precision.o constants_clubb.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/file_functions.F90
fill_holes.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/fill_holes.F90 grid_class.o clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/fill_holes.F90
fms_donner.o: $(SRCROOT)atmos_param_am3/donner_deep/fms_donner.F90 $(SRCROOT)atmos_param_am3/donner_deep/donner_nml.h atmos_tracer_utilities.o donner_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/donner_deep	$(SRCROOT)atmos_param_am3/donner_deep/fms_donner.F90
format_input.o: $(SRCROOT)atmos_param_am3/cosp/quickbeam/format_input.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_defs.H array_lib.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/cosp	$(SRCROOT)atmos_param_am3/cosp/quickbeam/format_input.F90
gamma_mg.o: $(SRCROOT)atmos_param_am3/strat_cloud/gamma_mg.F90 mg_const.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/strat_cloud/gamma_mg.F90
gas_tf.o: $(SRCROOT)atmos_param_am3/radiation/longwave/gas_tf.F90 longwave_utilities.o longwave_params.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/radiation/longwave/gas_tf.F90
gases.o: $(SRCROOT)atmos_param_am3/cosp/quickbeam/gases.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/cosp	$(SRCROOT)atmos_param_am3/cosp/quickbeam/gases.F90
get_random_number_stream.o: $(SRCROOT)atmos_param_am3/radiation/clouds/get_random_number_stream.F90 cloudrad_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/radiation/clouds/get_random_number_stream.F90
gmres_cache.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/gmres_cache.F90 clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/gmres_cache.F90
gmres_wrap.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/gmres_wrap.F90 gmres_cache.o clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/gmres_wrap.F90
grey_radiation.o: $(SRCROOT)atmos_param_am3/grey_radiation/grey_radiation.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/grey_radiation/grey_radiation.F90
grid_class.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/grid_class.F90 clubb_precision.o constants_clubb.o error_code.o file_functions.o interpolation.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/grid_class.F90
hydrostatic_module.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/hydrostatic_module.F90 constants_clubb.o grid_class.o clubb_precision.o interpolation.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/hydrostatic_module.F90
hyper_diffusion_4th_ord.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/hyper_diffusion_4th_ord.F90 clubb_precision.o grid_class.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/hyper_diffusion_4th_ord.F90
icarus.o: $(SRCROOT)atmos_param_am3/cosp/icarus-scops-3.7/icarus.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/cosp	$(SRCROOT)atmos_param_am3/cosp/icarus-scops-3.7/icarus.F90
ice_nucl.o: $(SRCROOT)atmos_shared_am3/tracer_driver/aer_ccn_act/ice_nucl.F90 aer_ccn_act_k.o aerosol_params.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/aer_ccn_act/ice_nucl.F90
input_names.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/input_names.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/input_names.F90
input_reader.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/input_reader.F90 clubb_precision.o constants_clubb.o input_names.o interpolation.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/input_reader.F90
interpolation.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/interpolation.F90 clubb_precision.o constants_clubb.o model_flags.o error_code.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/interpolation.F90
isccp_cloud_types.o: $(SRCROOT)atmos_param_am3/cosp/icarus-scops-3.7/isccp_cloud_types.F $(SRCROOT)atmos_param_am3/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/cosp	$(SRCROOT)atmos_param_am3/cosp/icarus-scops-3.7/isccp_cloud_types.F
isccp_clouds.o: $(SRCROOT)atmos_param_am3/radiation/clouds/isccp_clouds.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/radiation/clouds/isccp_clouds.F90
lapack_wrap.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/lapack_wrap.F90 constants_clubb.o error_code.o clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/lapack_wrap.F90
lidar_simulator.o: $(SRCROOT)atmos_param_am3/cosp/actsim/lidar_simulator.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/cosp	$(SRCROOT)atmos_param_am3/cosp/actsim/lidar_simulator.F90
lin_cloud_microphys.o: $(SRCROOT)atmos_param_am3/lin_cloud_microphys/lin_cloud_microphys.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/lin_cloud_microphys/lin_cloud_microphys.F90
llnl_stats.o: $(SRCROOT)atmos_param_am3/cosp/llnl/llnl_stats.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_defs.H cosp_constants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/cosp	$(SRCROOT)atmos_param_am3/cosp/llnl/llnl_stats.F90
lmd_ipsl_stats.o: $(SRCROOT)atmos_param_am3/cosp/actsim/lmd_ipsl_stats.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_defs.H llnl_stats.o cosp_constants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/cosp	$(SRCROOT)atmos_param_am3/cosp/actsim/lmd_ipsl_stats.F90
load_hydrometeor_classes.o: $(SRCROOT)atmos_param_am3/cosp/quickbeam/load_hydrometeor_classes.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_defs.H radar_simulator_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/cosp	$(SRCROOT)atmos_param_am3/cosp/quickbeam/load_hydrometeor_classes.F90
load_mie_table.o: $(SRCROOT)atmos_param_am3/cosp/quickbeam/load_mie_table.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_defs.H radar_simulator_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/cosp	$(SRCROOT)atmos_param_am3/cosp/quickbeam/load_mie_table.F90
longwave_clouds.o: $(SRCROOT)atmos_param_am3/radiation/longwave/longwave_clouds.F90 longwave_utilities.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/radiation/longwave/longwave_clouds.F90
longwave_driver.o: $(SRCROOT)atmos_param_am3/radiation/sea_esf_rad/longwave_driver.F90 rad_utilities.o aerosolrad_types.o sealw99.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/radiation/sea_esf_rad/longwave_driver.F90
longwave_fluxes.o: $(SRCROOT)atmos_param_am3/radiation/longwave/longwave_fluxes.F90 longwave_utilities.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/radiation/longwave/longwave_fluxes.F90
longwave_params.o: $(SRCROOT)atmos_param_am3/radiation/longwave/longwave_params.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/radiation/longwave/longwave_params.F90
longwave_tables.o: $(SRCROOT)atmos_param_am3/radiation/longwave/longwave_tables.F90 longwave_utilities.o longwave_params.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/radiation/longwave/longwave_tables.F90
longwave_utilities.o: $(SRCROOT)atmos_param_am3/radiation/longwave/longwave_utilities.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/radiation/longwave/longwave_utilities.F90
lscale_cond.o: $(SRCROOT)atmos_param_am3/lscale_cond/lscale_cond.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/lscale_cond/lscale_cond.F90
lw_gases_stdtf.o: $(SRCROOT)atmos_param_am3/radiation/longwave/lw_gases_stdtf.F90 longwave_utilities.o gas_tf.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/radiation/longwave/lw_gases_stdtf.F90
m_tracname.o: $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/m_tracname.F90 moz.mods.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/m_tracname.F90
math_lib.o: $(SRCROOT)atmos_param_am3/cosp/quickbeam/math_lib.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_defs.H mrgrnk.o array_lib.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/cosp	$(SRCROOT)atmos_param_am3/cosp/quickbeam/math_lib.F90
matrix_operations.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/matrix_operations.F90 clubb_precision.o error_code.o constants_clubb.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/matrix_operations.F90
mean_adv.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/mean_adv.F90 grid_class.o model_flags.o clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/mean_adv.F90
mg_const.o: $(SRCROOT)atmos_param_am3/strat_cloud/mg_const.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/strat_cloud/mg_const.F90
mg_drag.o: $(SRCROOT)atmos_param_am3/mg_drag/mg_drag.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/mg_drag/mg_drag.F90
micro_mg.o: $(SRCROOT)atmos_param_am3/strat_cloud/micro_mg.F90 gamma_mg.o strat_cloud_utilities.o simple_pdf.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/strat_cloud/micro_mg.F90
microphys_rad.o: $(SRCROOT)atmos_param_am3/radiation/clouds/microphys_rad.F90 cloudrad_types.o sealw99.o esfsw_driver.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/radiation/clouds/microphys_rad.F90
microphysics.o: $(SRCROOT)atmos_param_am3/strat_cloud/microphysics.F90 rotstayn_klein_mp.o morrison_gettelman_microp.o strat_cloud_utilities.o cldwat2m_micro.o micro_mg.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/strat_cloud/microphysics.F90
mixing_length.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/mixing_length.F90 constants_clubb.o parameters_tunable.o parameters_model.o grid_class.o numerical_check.o saturation.o error_code.o model_flags.o clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/mixing_length.F90
mo_chem_utls.o: $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_chem_utls.F90 moz.mods.o m_tracname.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_chem_utls.F90
mo_chemdr.o: $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_chemdr.F90 moz.mods.o mo_photo.o mo_exp_slv.o mo_imp_slv.o mo_rodas_slv.o mo_usrrxt.o mo_setinv.o moz.subs.o mo_setsox.o mo_fphoto.o mo_chem_utls.o strat_chem_utilities.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_chemdr.F90
mo_chemini.o: $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_chemini.F90 mo_photo.o mo_chem_utls.o mo_usrrxt.o moz.mods.o mo_exp_slv.o mo_imp_slv.o mo_rodas_slv.o mo_read_sim_chm.o mo_fphoto.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_chemini.F90
mo_exp_slv.o: $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_exp_slv.F90 moz.mods.o mo_chem_utls.o moz.mat.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_exp_slv.F90
mo_fastjx.o: $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_fastjx.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_fastjx.F90
mo_fphoto.o: $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_fphoto.F90 mo_fastjx.o mo_chem_utls.o moz.mods.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_fphoto.F90
mo_hook.o: $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_hook.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_hook.F90
mo_imp_slv.o: $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_imp_slv.F90 moz.mods.o mo_chem_utls.o m_tracname.o moz.mat.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_imp_slv.F90
mo_jpl.o: $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_jpl.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_jpl.F90
mo_photo.o: $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_photo.F90 mo_chem_utls.o moz.mods.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_photo.F90
mo_read_sim_chm.o: $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_read_sim_chm.F90 moz.mods.o m_tracname.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_read_sim_chm.F90
mo_rodas_slv.o: $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_rodas_slv.F90 moz.mods.o mo_chem_utls.o moz.mat.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_rodas_slv.F90
mo_setinv.o: $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_setinv.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_setinv.F90
mo_setsox.o: $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_setsox.F90 mo_chem_utls.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_setsox.F90
mo_usrrxt.o: $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_usrrxt.F90 strat_chem_utilities.o mo_chem_utls.o moz.mods.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_usrrxt.F90
model_flags.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/model_flags.F90 constants_clubb.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/model_flags.F90
modis_simulator.o: $(SRCROOT)atmos_param_am3/cosp/MODIS_simulator/modis_simulator.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_defs.H cosp_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/cosp	$(SRCROOT)atmos_param_am3/cosp/MODIS_simulator/modis_simulator.F90
moist_conv.o: $(SRCROOT)atmos_param_am3/moist_conv/moist_conv.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/moist_conv/moist_conv.F90
moist_processes.o: $(SRCROOT)atmos_param_am3/moist_processes/moist_processes.F90 betts_miller.o bm_massflux.o bm_omp.o donner_deep.o moist_conv.o lscale_cond.o uw_conv.o lin_cloud_microphys.o ras.o dry_adj.o strat_cloud.o detr_ice_num.o MG_microp_3D.o CLUBB_driver_SCM.o rh_clouds.o diag_cloud.o diag_integral.o cu_mo_trans.o mo_hook.o aerosol_types.o moist_processes_utils.o moistproc_kernels.o atmos_tracer_utilities.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/moist_processes/moist_processes.F90
moist_processes_utils.o: $(SRCROOT)atmos_param_am3/moist_processes/moist_processes_utils.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/moist_processes/moist_processes_utils.F90
moistproc_kernels.o: $(SRCROOT)atmos_param_am3/moist_processes/moistproc_kernels.F90 betts_miller.o bm_massflux.o bm_omp.o diag_cloud.o donner_deep.o moist_conv.o lscale_cond.o uw_conv.o lin_cloud_microphys.o ras.o strat_cloud.o rh_clouds.o cu_mo_trans.o atmos_tracer_utilities.o mo_hook.o aerosol_types.o moist_processes_utils.o detr_ice_num.o MG_microp_3D.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/moist_processes/moistproc_kernels.F90
monin_obukhov.o: $(SRCROOT)atmos_param_am3/monin_obukhov/monin_obukhov.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/monin_obukhov/monin_obukhov.F90
monin_obukhov_kernel.o: $(SRCROOT)atmos_param_am3/monin_obukhov/monin_obukhov_kernel.F90 $(SRCROOT)fms/include/fms_platform.h $(SRCROOT)atmos_param_am3/monin_obukhov/monin_obukhov_interfaces.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include -I$(SRCROOT)atmos_param_am3/monin_obukhov	$(SRCROOT)atmos_param_am3/monin_obukhov/monin_obukhov_kernel.F90
mono_flux_limiter.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/mono_flux_limiter.F90 grid_class.o constants_clubb.o clubb_precision.o error_code.o fill_holes.o stats_type.o stats_variables.o mean_adv.o lapack_wrap.o anl_erf.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/mono_flux_limiter.F90
morrison_gettelman_microp.o: $(SRCROOT)atmos_param_am3/strat_cloud/morrison_gettelman_microp.F90 polysvp.o strat_cloud_utilities.o gamma_mg.o mg_const.o simple_pdf.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/strat_cloud/morrison_gettelman_microp.F90
moz.mat.o: $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/moz.mat.F90 moz.mods.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/moz.mat.F90
moz.mods.o: $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/moz.mods.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/moz.mods.F90
moz.subs.o: $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/moz.subs.F90 moz.mods.o mo_jpl.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/moz.subs.F90
mrgrnk.o: $(SRCROOT)atmos_param_am3/cosp/quickbeam/mrgrnk.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/cosp	$(SRCROOT)atmos_param_am3/cosp/quickbeam/mrgrnk.F90
mt95.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/mt95.f90
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/mt95.f90
my25_turb.o: $(SRCROOT)atmos_param_am3/my25_turb/my25_turb.F90 monin_obukhov.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/my25_turb/my25_turb.F90
nc_cond.o: $(SRCROOT)atmos_param_am3/strat_cloud/nc_cond.F90 betaDistribution.o strat_cloud_utilities.o polysvp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/strat_cloud/nc_cond.F90
nonfms_donner.o: $(SRCROOT)atmos_param_am3/donner_deep/nonfms_donner.F90 $(SRCROOT)atmos_param_am3/donner_deep/donner_nml.h donner_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/donner_deep	$(SRCROOT)atmos_param_am3/donner_deep/nonfms_donner.F90
numerical_check.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/numerical_check.F90 grid_class.o clubb_precision.o parameters_model.o pdf_parameter_module.o stats_variables.o constants_clubb.o variables_diagnostic_module.o variables_prognostic_module.o parameters_microphys.o error_code.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/numerical_check.F90
optical_path.o: $(SRCROOT)atmos_param_am3/radiation/longwave/optical_path.F90 longwave_utilities.o longwave_params.o lw_gases_stdtf.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/radiation/longwave/optical_path.F90
optics_lib.o: $(SRCROOT)atmos_param_am3/cosp/quickbeam/optics_lib.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/cosp	$(SRCROOT)atmos_param_am3/cosp/quickbeam/optics_lib.F90
output_grads.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/output_grads.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/recl.inc clubb_precision.o constants_clubb.o stat_file_module.o model_flags.o endian.o calendar.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/clubb/CLUBB_core	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/output_grads.F90
output_netcdf.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/output_netcdf.F90 stat_file_module.o clubb_precision.o constants_clubb.o parameters_model.o parameters_tunable.o parameter_indices.o model_flags.o parameters_microphys.o parameters_radiation.o calendar.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/output_netcdf.F90
ozone.o: $(SRCROOT)atmos_param_am3/radiation/radiative_gases/ozone.F90 rad_utilities.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/radiation/radiative_gases/ozone.F90
parameter_indices.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/parameter_indices.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/parameter_indices.F90
parameters_microphys.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/parameters_microphys.F90 clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/parameters_microphys.F90
parameters_model.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/parameters_model.F90 clubb_precision.o constants_clubb.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/parameters_model.F90
parameters_radiation.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/parameters_radiation.F90 clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/parameters_radiation.F90
parameters_tunable.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/parameters_tunable.F90 parameter_indices.o grid_class.o clubb_precision.o constants_clubb.o error_code.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/parameters_tunable.F90
pdf_closure_module.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/pdf_closure_module.F90 constants_clubb.o parameters_model.o parameters_tunable.o pdf_parameter_module.o anl_erf.o numerical_check.o saturation.o error_code.o stats_variables.o clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/pdf_closure_module.F90
pdf_parameter_module.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/pdf_parameter_module.F90 clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/pdf_parameter_module.F90
pf_to_mr.o: $(SRCROOT)atmos_param_am3/cosp/llnl/pf_to_mr.F $(SRCROOT)atmos_param_am3/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/cosp	$(SRCROOT)atmos_param_am3/cosp/llnl/pf_to_mr.F
physics_driver.o: $(SRCROOT)atmos_param_am3/physics_driver/physics_driver.F90 atmos_tracer_driver.o rad_utilities.o aerosol_types.o physics_radiation_exch.o physics_types.o aerosol.o cosp_driver.o moist_processes.o vert_turb_driver.o vert_diff_driver.o damping_driver.o grey_radiation.o CLUBB_driver_SCM.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/physics_driver/physics_driver.F90
physics_radiation_exch.o: $(SRCROOT)atmos_param_am3/physics_radiation_exch/physics_radiation_exch.F90 $(SRCROOT)fms/include/fms_platform.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)atmos_param_am3/physics_radiation_exch/physics_radiation_exch.F90
physics_types.o: $(SRCROOT)atmos_param_am3/physics_driver/physics_types.F90 $(SRCROOT)fms/include/fms_platform.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)atmos_param_am3/physics_driver/physics_types.F90
polysvp.o: $(SRCROOT)atmos_param_am3/strat_cloud/polysvp.F90 mg_const.o strat_cloud_utilities.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/strat_cloud/polysvp.F90
pos_definite_module.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/pos_definite_module.F90 grid_class.o constants_clubb.o clubb_precision.o error_code.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/pos_definite_module.F90
prec_scops.o: $(SRCROOT)atmos_param_am3/cosp/llnl/prec_scops.F $(SRCROOT)atmos_param_am3/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/cosp	$(SRCROOT)atmos_param_am3/cosp/llnl/prec_scops.F
qe_moist_convection.o: $(SRCROOT)atmos_param_am3/qe_moist_convection/qe_moist_convection.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/qe_moist_convection/qe_moist_convection.F90
rad_output_file.o: $(SRCROOT)atmos_param_am3/radiation/util/rad_output_file.F90 aerosol_types.o rad_utilities.o aerosolrad_types.o cloudrad_types.o esfsw_driver.o sealw99.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/radiation/util/rad_output_file.F90
rad_utilities.o: $(SRCROOT)atmos_param_am3/radiation/util/rad_utilities.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/radiation/util/rad_utilities.F90
radar_simulator.o: $(SRCROOT)atmos_param_am3/cosp/quickbeam/radar_simulator.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_defs.H mrgrnk.o array_lib.o math_lib.o optics_lib.o radar_simulator_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/cosp	$(SRCROOT)atmos_param_am3/cosp/quickbeam/radar_simulator.F90
radar_simulator_types.o: $(SRCROOT)atmos_param_am3/cosp/quickbeam/radar_simulator_types.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/cosp	$(SRCROOT)atmos_param_am3/cosp/quickbeam/radar_simulator_types.F90
radiation_diag.o: $(SRCROOT)atmos_param_am3/radiation/util/radiation_diag.F90 radiation_driver.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/radiation/util/radiation_diag.F90
radiation_driver.o: $(SRCROOT)atmos_param_am3/radiation/driver/radiation_driver.F90 radiation_types.o physics_radiation_exch.o rad_utilities.o aerosol_types.o aerosolrad_types.o cloudrad_types.o sea_esf_rad.o rad_output_file.o aerosolrad_driver.o cloudrad_driver.o cloudrad_diagnostics.o radiative_gases.o radiation_driver_diag.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/radiation/driver/radiation_driver.F90
radiation_driver_diag.o: $(SRCROOT)atmos_param_am3/radiation/driver/radiation_driver_diag.F90 rad_utilities.o aerosolrad_types.o esfsw_driver.o sealw99.o diag_integral.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/radiation/driver/radiation_driver_diag.F90
radiation_types.o: $(SRCROOT)atmos_param_am3/radiation/driver/radiation_types.F90 $(SRCROOT)fms/include/fms_platform.h atmos_co2.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)atmos_param_am3/radiation/driver/radiation_types.F90
radiative_gases.o: $(SRCROOT)atmos_param_am3/radiation/radiative_gases/radiative_gases.F90 rad_utilities.o sealw99.o ozone.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/radiation/radiative_gases/radiative_gases.F90
ras.o: $(SRCROOT)atmos_param_am3/ras/ras.F90 aerosol_types.o aer_ccn_act.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/ras/ras.F90
rh_clouds.o: $(SRCROOT)atmos_param_am3/rh_clouds/rh_clouds.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/rh_clouds/rh_clouds.F90
rotstayn_klein_mp.o: $(SRCROOT)atmos_param_am3/strat_cloud/rotstayn_klein_mp.F90 cloud_generator.o aer_in_act.o polysvp.o strat_cloud_utilities.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/strat_cloud/rotstayn_klein_mp.F90
saturation.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/saturation.F90 model_flags.o clubb_precision.o constants_clubb.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/saturation.F90
scops.o: $(SRCROOT)atmos_param_am3/cosp/icarus-scops-3.7/scops.F $(SRCROOT)atmos_param_am3/cosp/cosp_defs.H $(SRCROOT)atmos_param_am3/cosp/icarus-scops-3.7/congvec.H
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/cosp -I$(SRCROOT)atmos_param_am3/cosp/icarus-scops-3.7	$(SRCROOT)atmos_param_am3/cosp/icarus-scops-3.7/scops.F
sea_esf_rad.o: $(SRCROOT)atmos_param_am3/radiation/sea_esf_rad/sea_esf_rad.F90 rad_utilities.o aerosolrad_types.o longwave_driver.o shortwave_driver.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/radiation/sea_esf_rad/sea_esf_rad.F90
sealw99.o: $(SRCROOT)atmos_param_am3/radiation/longwave/sealw99.F90 longwave_utilities.o longwave_params.o longwave_clouds.o longwave_fluxes.o longwave_tables.o optical_path.o gas_tf.o lw_gases_stdtf.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/radiation/longwave/sealw99.F90
shallow_conv.o: $(SRCROOT)atmos_param_am3/shallow_conv/shallow_conv.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/shallow_conv/shallow_conv.F90
shallow_physics.o: $(SRCROOT)atmos_param_am3/shallow_physics/shallow_physics.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/shallow_physics/shallow_physics.F90
shortwave_driver.o: $(SRCROOT)atmos_param_am3/radiation/sea_esf_rad/shortwave_driver.F90 rad_utilities.o aerosolrad_types.o esfsw_driver.o solar_data.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/radiation/sea_esf_rad/shortwave_driver.F90
sigma_sqd_w_module.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/sigma_sqd_w_module.F90 constants_clubb.o clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/sigma_sqd_w_module.F90
simple_pdf.o: $(SRCROOT)atmos_param_am3/strat_cloud/simple_pdf.F90 betaDistribution.o strat_cloud_utilities.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/strat_cloud/simple_pdf.F90
solar_data.o: $(SRCROOT)atmos_param_am3/radiation/sea_esf_rad/solar_data.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/radiation/sea_esf_rad/solar_data.F90
sponge_layer_damping.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/sponge_layer_damping.F90 clubb_precision.o grid_class.o constants_clubb.o interpolation.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/sponge_layer_damping.F90
stable_bl_turb.o: $(SRCROOT)atmos_param_am3/stable_bl_turb/stable_bl_turb.F90 monin_obukhov.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/stable_bl_turb/stable_bl_turb.F90
stat_file_module.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/stat_file_module.F90 clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/stat_file_module.F90
stats_LH_zt.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/stats_LH_zt.F90 constants_clubb.o stats_variables.o stats_type.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/stats_LH_zt.F90
stats_rad_zm.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/stats_rad_zm.F90 constants_clubb.o stats_variables.o stats_type.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/stats_rad_zm.F90
stats_rad_zt.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/stats_rad_zt.F90 constants_clubb.o stats_variables.o stats_type.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/stats_rad_zt.F90
stats_sfc.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/stats_sfc.F90 constants_clubb.o stats_variables.o stats_type.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/stats_sfc.F90
stats_subs.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/stats_subs.F90 stats_variables.o clubb_precision.o output_grads.o output_netcdf.o stats_zm.o stats_zt.o stats_LH_zt.o stats_rad_zt.o stats_rad_zm.o stats_sfc.o error_code.o constants_clubb.o parameters_microphys.o grid_class.o variables_diagnostic_module.o pdf_parameter_module.o T_in_K_module.o parameters_model.o stats_type.o fill_holes.o interpolation.o saturation.o array_index.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/stats_subs.F90
stats_type.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/stats_type.F90 stat_file_module.o clubb_precision.o grid_class.o error_code.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/stats_type.F90
stats_variables.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/stats_variables.F90 stats_type.o clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/stats_variables.F90
stats_zm.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/stats_zm.F90 constants_clubb.o stats_variables.o stats_type.o parameters_model.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/stats_zm.F90
stats_zt.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/stats_zt.F90 constants_clubb.o stats_variables.o stats_type.o parameters_model.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/stats_zt.F90
strat_chem_driver.o: $(SRCROOT)atmos_shared_am3/tracer_driver/stratchem/strat_chem_driver.F90 strat_chem_model.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/stratchem/strat_chem_driver.F90
strat_chem_model.o: $(SRCROOT)atmos_shared_am3/tracer_driver/stratchem/strat_chem_model.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/stratchem/strat_chem_model.F90
strat_chem_utilities.o: $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/strat_chem_utilities.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/strat_chem_utilities.F90
strat_cloud.o: $(SRCROOT)atmos_param_am3/strat_cloud/strat_cloud.F90 $(SRCROOT)atmos_param_am3/strat_cloud/strat_nml.h cloud_rad.o aerosol_types.o microphysics.o nc_cond.o polysvp.o aerosol_cloud.o strat_cloud_utilities.o strat_netcdf.o strat_cloud_legacy.o check_nan.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/strat_cloud	$(SRCROOT)atmos_param_am3/strat_cloud/strat_cloud.F90
strat_cloud_legacy.o: $(SRCROOT)atmos_param_am3/strat_cloud/strat_cloud_legacy.F90 betaDistribution.o aer_ccn_act.o aer_in_act.o cloud_rad.o cloud_generator.o aerosol_types.o strat_cloud_utilities.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/strat_cloud/strat_cloud_legacy.F90
strat_cloud_utilities.o: $(SRCROOT)atmos_param_am3/strat_cloud/strat_cloud_utilities.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/strat_cloud/strat_cloud_utilities.F90
strat_clouds_W.o: $(SRCROOT)atmos_param_am3/radiation/clouds/strat_clouds_W.F90 aerosol_types.o cloud_rad.o cloudrad_types.o get_random_number_stream.o cloud_generator.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/radiation/clouds/strat_clouds_W.F90
strat_netcdf.o: $(SRCROOT)atmos_param_am3/strat_cloud/strat_netcdf.F90 strat_cloud_utilities.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/strat_cloud/strat_netcdf.F90
surface_varnce_module.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/surface_varnce_module.F90 parameters_model.o constants_clubb.o numerical_check.o error_code.o array_index.o stats_type.o clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/surface_varnce_module.F90
tke_turb.o: $(SRCROOT)atmos_param_am3/tke_turb/tke_turb.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/tke_turb/tke_turb.F90
topo_drag.o: $(SRCROOT)atmos_param_am3/topo_drag/topo_drag.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/topo_drag/topo_drag.F90
tropchem_driver.o: $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/tropchem_driver.F90 atmos_tracer_utilities.o mo_chemdr.o mo_chemini.o m_tracname.o moz.mods.o mo_hook.o strat_chem_utilities.o mo_chem_utls.o atmos_sulfate.o shortwave_driver.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/tropchem_driver.F90
two_stream_gray_rad.o: $(SRCROOT)atmos_param_am3/two_stream_gray_rad/two_stream_gray_rad.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/two_stream_gray_rad/two_stream_gray_rad.F90
uw_clouds_W.o: $(SRCROOT)atmos_param_am3/radiation/clouds/uw_clouds_W.F90 cloudrad_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/radiation/clouds/uw_clouds_W.F90
uw_conv.o: $(SRCROOT)atmos_param_am3/shallow_cu/uw_conv.F90 $(SRCROOT)fms/include/fms_platform.h atmos_tracer_utilities.o aerosol_types.o aer_ccn_act.o conv_utilities.o conv_utilities_k.o conv_plumes_k.o conv_closures.o deep_conv.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)atmos_param_am3/shallow_cu/uw_conv.F90
variables_diagnostic_module.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/variables_diagnostic_module.F90 pdf_parameter_module.o clubb_precision.o constants_clubb.o parameters_model.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/variables_diagnostic_module.F90
variables_prognostic_module.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/variables_prognostic_module.F90 pdf_parameter_module.o clubb_precision.o constants_clubb.o parameters_model.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/variables_prognostic_module.F90
variables_radiation_module.o: $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/variables_radiation_module.F90 clubb_precision.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/clubb/CLUBB_core/variables_radiation_module.F90
vert_advection.o: $(SRCROOT)atmos_shared_am3/vert_advection/vert_advection.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared_am3/vert_advection/vert_advection.F90
vert_diff.o: $(SRCROOT)atmos_param_am3/vert_diff/vert_diff.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/vert_diff/vert_diff.F90
vert_diff_driver.o: $(SRCROOT)atmos_param_am3/vert_diff_driver/vert_diff_driver.F90 vert_diff.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/vert_diff_driver/vert_diff_driver.F90
vert_turb_driver.o: $(SRCROOT)atmos_param_am3/vert_turb_driver/vert_turb_driver.F90 my25_turb.o tke_turb.o diffusivity.o edt.o strat_cloud.o shallow_conv.o stable_bl_turb.o entrain.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/vert_turb_driver/vert_turb_driver.F90
wet_deposition_0D.o: $(SRCROOT)atmos_param_am3/donner_deep/wet_deposition_0D.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param_am3/donner_deep/wet_deposition_0D.F90
zeff.o: $(SRCROOT)atmos_param_am3/cosp/quickbeam/zeff.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_defs.H math_lib.o optics_lib.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param_am3/cosp	$(SRCROOT)atmos_param_am3/cosp/quickbeam/zeff.F90

SRC = $(SRCROOT)atmos_param_am3/radiation/sea_esf_rad/longwave_driver.F90 $(SRCROOT)atmos_param_am3/radiation/shortwave/esfsw_parameters.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_photo.F90 $(SRCROOT)atmos_param_am3/radiation/shortwave/esfsw_bands.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_rttov_simulator.F90 $(SRCROOT)atmos_param_am3/cosp/quickbeam/radar_simulator.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/atmos_soa.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/model_flags.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/stats_type.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_types.F90 $(SRCROOT)atmos_param_am3/radiation/aerosols/aerosolrad_driver.F90 $(SRCROOT)atmos_param_am3/betts_miller/bm_massflux.F90 $(SRCROOT)atmos_param_am3/radiation/aerosols/aerosolrad_package.F90 $(SRCROOT)atmos_param_am3/cosp/icarus-scops-3.7/isccp_cloud_types.F $(SRCROOT)atmos_param_am3/cosp/quickbeam/atmos_lib.F90 $(SRCROOT)atmos_param_am3/radiation/longwave/optical_path.F90 $(SRCROOT)atmos_param_am3/strat_cloud/strat_cloud_utilities.F90 $(SRCROOT)atmos_shared_am3/vert_advection/vert_advection.F90 $(SRCROOT)atmos_param_am3/cosp/llnl/llnl_stats.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/advance_xm_wpxp_module.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/gmres_cache.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_isccp_simulator.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_hook.F90 $(SRCROOT)atmos_param_am3/stable_bl_turb/stable_bl_turb.F90 $(SRCROOT)atmos_param_am3/donner_deep/donner_rad_k.F90 $(SRCROOT)atmos_param_am3/grey_radiation/grey_radiation.F90 $(SRCROOT)atmos_param_am3/entrain/entrain.F90 $(SRCROOT)atmos_param_am3/radiation/radiative_gases/radiative_gases.F90 $(SRCROOT)atmos_param_am3/topo_drag/topo_drag.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/array_index.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/atmos_dust.F90 $(SRCROOT)atmos_param_am3/cosp/MODIS_simulator/modis_simulator.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_setsox.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/variables_diagnostic_module.F90 $(SRCROOT)atmos_param_am3/betts_miller/betts_miller.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_jpl.F90 $(SRCROOT)atmos_param_am3/cosp/cosp.F90 $(SRCROOT)atmos_param_am3/shallow_conv/shallow_conv.F90 $(SRCROOT)atmos_param_am3/cosp/llnl/cosp_radar.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/mixing_length.F90 $(SRCROOT)atmos_param_am3/strat_cloud/gamma_mg.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/saturation.F90 $(SRCROOT)atmos_param_am3/clubb/MG_microp_3D.F90 $(SRCROOT)atmos_param_am3/monin_obukhov/monin_obukhov_kernel.F90 $(SRCROOT)atmos_param_am3/moist_processes/detr_ice_num.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/mt95.f90 $(SRCROOT)atmos_param_am3/radiation/clouds/bulkphys_rad.F90 $(SRCROOT)atmos_param_am3/radiation/clouds/cloudrad_driver.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/atmos_sulfate.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/atmos_convection_tracer.F90 $(SRCROOT)atmos_param_am3/cosp/MISR_simulator/MISR_simulator.F $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_setinv.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_chemdr.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_driver.F90 $(SRCROOT)atmos_param_am3/donner_deep/fms_donner.F90 $(SRCROOT)atmos_param_am3/radiation/shortwave/esfsw_driver.F90 $(SRCROOT)atmos_param_am3/cosp/quickbeam/mrgrnk.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/atmos_regional_tracer_driver.F90 $(SRCROOT)atmos_param_am3/cosp/actsim/lmd_ipsl_stats.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/calendar.F90 $(SRCROOT)atmos_param_am3/physics_driver/physics_types.F90 $(SRCROOT)atmos_param_am3/radiation/clouds/uw_clouds_W.F90 $(SRCROOT)atmos_param_am3/radiation/driver/radiation_driver_diag.F90 $(SRCROOT)atmos_param_am3/radiation/sea_esf_rad/solar_data.F90 $(SRCROOT)atmos_param_am3/clubb/alt_cloud.F90 $(SRCROOT)atmos_param_am3/diag_integral/diag_integral.F90 $(SRCROOT)atmos_param_am3/cosp/quickbeam/load_mie_table.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/sigma_sqd_w_module.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/atmos_carbon_aerosol.F90 $(SRCROOT)atmos_param_am3/radiation/sea_esf_rad/shortwave_driver.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_3D_var.F90 $(SRCROOT)atmos_param_am3/shallow_cu/conv_utilities_k.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_constants.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/parameters_radiation.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/parameter_indices.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/hydrostatic_module.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_diagnostics.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/aer_ccn_act/ice_nucl.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/stats_rad_zm.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/stats_subs.F90 $(SRCROOT)atmos_param_am3/strat_cloud/mg_const.F90 $(SRCROOT)atmos_param_am3/radiation/clouds/cloudrad_diagnostics.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/output_grads.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/m_tracname.F90 $(SRCROOT)atmos_param_am3/strat_cloud/strat_netcdf.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_driver_SCM.F90 $(SRCROOT)atmos_param_am3/ras/ras.F90 $(SRCROOT)atmos_param_am3/radiation/aerosols/aerosolrad_types.F90 $(SRCROOT)atmos_param_am3/moist_conv/moist_conv.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/gmres_wrap.F90 $(SRCROOT)atmos_param_am3/cosp/quickbeam/optics_lib.F90 $(SRCROOT)atmos_param_am3/strat_cloud/microphysics.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/interpolation.F90 $(SRCROOT)atmos_param_am3/lin_cloud_microphys/lin_cloud_microphys.F90 $(SRCROOT)atmos_param_am3/shallow_cu/conv_utilities.F90 $(SRCROOT)atmos_param_am3/radiation/sea_esf_rad/sea_esf_rad.F90 $(SRCROOT)atmos_param_am3/cosp/quickbeam/radar_simulator_types.F90 $(SRCROOT)atmos_param_am3/strat_cloud/morrison_gettelman_microp.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/sponge_layer_damping.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_lidar.F90 $(SRCROOT)atmos_param_am3/radiation/clouds/cloudrad_package.F90 $(SRCROOT)atmos_param_am3/cloud_zonal/cloud_zonal.F90 $(SRCROOT)atmos_param_am3/donner_deep/donner_meso_k.F90 $(SRCROOT)atmos_param_am3/dry_adj/dry_adj.F90 $(SRCROOT)atmos_param_am3/shallow_cu/conv_closures.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/file_functions.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/mean_adv.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/atmos_co2.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/error_code.F90 $(SRCROOT)atmos_shared_am3/atmos_nudge/atmos_nudge.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/output_netcdf.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/stratchem/strat_chem_model.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/advance_helper_module.F90 $(SRCROOT)atmos_param_am3/cloud_rad/cloud_rad.F90 $(SRCROOT)atmos_param_am3/cosp/quickbeam/math_lib.F90 $(SRCROOT)atmos_param_am3/radiation/longwave/longwave_tables.F90 $(SRCROOT)atmos_param_am3/physics_radiation_exch/physics_radiation_exch.F90 $(SRCROOT)atmos_param_am3/donner_deep/donner_lite_k.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/atmos_ch3i.F90 $(SRCROOT)atmos_param_am3/cosp/quickbeam/gases.F90 $(SRCROOT)atmos_param_am3/monin_obukhov/monin_obukhov.F90 $(SRCROOT)atmos_param_am3/donner_deep/nonfms_donner.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/input_names.F90 $(SRCROOT)atmos_param_am3/edt/edt.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/pos_definite_module.F90 $(SRCROOT)atmos_param_am3/strat_cloud/strat_cloud_legacy.F90 $(SRCROOT)atmos_param_am3/strat_cloud/polysvp.F90 $(SRCROOT)atmos_param_am3/donner_deep/donner_cloud_model_k.F90 $(SRCROOT)atmos_param_am3/radiation/radiative_gases/ozone.F90 $(SRCROOT)atmos_param_am3/shallow_cu/deep_conv.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/lapack_wrap.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/atmos_age_tracer.F90 $(SRCROOT)atmos_param_am3/shallow_cu/conv_plumes.F90 $(SRCROOT)atmos_param_am3/donner_deep/donner_types.F90 $(SRCROOT)atmos_param_am3/radiation/clouds/cloudrad_types.F90 $(SRCROOT)atmos_param_am3/radiation/clouds/cloud_spec.F90 $(SRCROOT)atmos_param_am3/radiation/util/rad_output_file.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/mono_flux_limiter.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/clip_explicit.F90 $(SRCROOT)atmos_param_am3/strat_cloud/strat_cloud.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_imp_slv.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_chem_utls.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_modis_simulator.F90 $(SRCROOT)atmos_param_am3/radiation/longwave/longwave_clouds.F90 $(SRCROOT)atmos_param_am3/shallow_cu/conv_plumes_k.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/stats_variables.F90 $(SRCROOT)atmos_param_am3/cosp/quickbeam/format_input.F90 $(SRCROOT)atmos_param_am3/qe_moist_convection/qe_moist_convection.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/atmos_sulfur_hex.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/stats_zm.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/pdf_parameter_module.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/anl_erf.F90 $(SRCROOT)atmos_param_am3/diag_cloud/diag_cloud.F90 $(SRCROOT)atmos_param_am3/donner_deep/donner_cape_k.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/moz.mods.F90 $(SRCROOT)atmos_param_am3/donner_deep/donner_deep_miz.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/parameters_microphys.F90 $(SRCROOT)atmos_param_am3/tke_turb/tke_turb.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/extrapolation.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_stats.F90 $(SRCROOT)atmos_param_am3/radiation/longwave/longwave_utilities.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/stat_file_module.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/pdf_closure_module.F90 $(SRCROOT)atmos_param_am3/radiation/longwave/longwave_params.F90 $(SRCROOT)atmos_param_am3/donner_deep/donner_deep_k.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_fastjx.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/stats_rad_zt.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/variables_prognostic_module.F90 $(SRCROOT)atmos_param_am3/radiation/clouds/donner_deep_clouds_W.F90 $(SRCROOT)atmos_param_am3/radiation/shortwave/esfsw_utilities.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/surface_varnce_module.F90 $(SRCROOT)atmos_param_am3/cosp/quickbeam/array_lib.F90 $(SRCROOT)atmos_param_am3/vert_turb_driver/vert_turb_driver.F90 $(SRCROOT)atmos_param_am3/radiation/longwave/lw_gases_stdtf.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/stats_zt.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/aer_ccn_act/aer_ccn_act.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/moz.mat.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/input_reader.F90 $(SRCROOT)atmos_param_am3/vert_diff_driver/vert_diff_driver.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_simulator.F90 $(SRCROOT)atmos_param_am3/moist_processes/moist_processes.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/clubb_precision.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/diffusion.F90 $(SRCROOT)atmos_param_am3/radiation/util/radiation_diag.F90 $(SRCROOT)atmos_param_am3/strat_cloud/rotstayn_klein_mp.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/aer_ccn_act/aerosol_params.F90 $(SRCROOT)atmos_param_am3/cosp/llnl/pf_to_mr.F $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/clubb_core.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_usrrxt.F90 $(SRCROOT)atmos_param_am3/radiation/clouds/microphys_rad.F90 $(SRCROOT)atmos_param_am3/radiation/longwave/longwave_fluxes.F90 $(SRCROOT)atmos_param_am3/lscale_cond/lscale_cond.F90 $(SRCROOT)atmos_param_am3/cloud_generator/betaDistribution.F90 $(SRCROOT)atmos_param_am3/cosp/icarus-scops-3.7/scops.F $(SRCROOT)atmos_param_am3/moist_processes/moistproc_kernels.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/moz.subs.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_read_sim_chm.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_fphoto.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/parameters_tunable.F90 $(SRCROOT)atmos_param_am3/rh_clouds/rh_clouds.F90 $(SRCROOT)atmos_shared_am3/aerosol/aerosol.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/T_in_K_module.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/variables_radiation_module.F90 $(SRCROOT)atmos_param_am3/moist_processes/moist_processes_utils.F90 $(SRCROOT)atmos_param_am3/clubb/advance_sclrm_Nd_module.F90 $(SRCROOT)atmos_param_am3/cloud_obs/cloud_obs.F90 $(SRCROOT)atmos_param_am3/radiation/util/rad_utilities.F90 $(SRCROOT)atmos_param_am3/diffusivity/diffusivity.F90 $(SRCROOT)atmos_param_am3/cloud_generator/cloud_generator.F90 $(SRCROOT)atmos_param_am3/cosp/quickbeam/load_hydrometeor_classes.F90 $(SRCROOT)atmos_shared_am3/aerosol/aerosol_types.F90 $(SRCROOT)atmos_param_am3/donner_deep/cumulus_closure_k.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_utils.F90 $(SRCROOT)atmos_param_am3/radiation/longwave/sealw99.F90 $(SRCROOT)atmos_param_am3/radiation/clouds/isccp_clouds.F90 $(SRCROOT)atmos_param_am3/clouds/clouds.F90 $(SRCROOT)atmos_param_am3/my25_turb/my25_turb.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_chemini.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/atmos_radon.F90 $(SRCROOT)atmos_param_am3/cosp/llnl/prec_scops.F $(SRCROOT)atmos_param_am3/damping_driver/damping_driver.F90 $(SRCROOT)atmos_param_am3/cu_mo_trans/cu_mo_trans.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/grid_class.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/tropchem_driver.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/stats_sfc.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_rodas_slv.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/csr_matrix_class_3array.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/fill_holes.F90 $(SRCROOT)atmos_param_am3/radiation/longwave/gas_tf.F90 $(SRCROOT)atmos_param_am3/donner_deep/donner_utilities_k.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/parameters_model.F90 $(SRCROOT)atmos_param_am3/donner_deep/donner_deep.F90 $(SRCROOT)atmos_param_am3/radiation/driver/radiation_driver.F90 $(SRCROOT)atmos_param_am3/strat_cloud/aerosol_cloud.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/aer_ccn_act/aer_in_act.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_io.F90 $(SRCROOT)atmos_param_am3/shallow_cu/uw_conv.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/corr_matrix_module.F90 $(SRCROOT)atmos_param_am3/betts_miller/bm_omp.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/aer_ccn_act/aer_ccn_act_k.F90 $(SRCROOT)atmos_param_am3/diag_cloud_rad/diag_cloud_rad.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/constants_clubb.F90 $(SRCROOT)atmos_param_am3/radiation/clouds/get_random_number_stream.F90 $(SRCROOT)atmos_param_am3/shallow_physics/shallow_physics.F90 $(SRCROOT)atmos_param_am3/strat_cloud/check_nan.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/hyper_diffusion_4th_ord.F90 $(SRCROOT)atmos_param_am3/donner_deep/donner_lscloud_k.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/atmos_tracer_utilities.F90 $(SRCROOT)atmos_param_am3/physics_driver/physics_driver.F90 $(SRCROOT)atmos_param_am3/mg_drag/mg_drag.F90 $(SRCROOT)atmos_param_am3/cosp/icarus-scops-3.7/icarus.F90 $(SRCROOT)atmos_param_am3/donner_deep/wet_deposition_0D.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_misr_simulator.F90 $(SRCROOT)atmos_param_am3/two_stream_gray_rad/two_stream_gray_rad.F90 $(SRCROOT)atmos_param_am3/cosp/quickbeam/dsd.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/mo_exp_slv.F90 $(SRCROOT)atmos_shared_am3/atmos_cmip_diag/atmos_cmip_diag.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/advance_xp2_xpyp_module.F90 $(SRCROOT)atmos_param_am3/cosp/actsim/lidar_simulator.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/numerical_check.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/stratchem/strat_chem_driver.F90 $(SRCROOT)atmos_param_am3/cg_drag/cg_drag.F90 $(SRCROOT)atmos_param_am3/strat_cloud/nc_cond.F90 $(SRCROOT)atmos_param_am3/radiation/driver/radiation_types.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/stats_LH_zt.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/Skw_module.F90 $(SRCROOT)atmos_param_am3/strat_cloud/cldwat2m_micro.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/advance_wp2_wp3_module.F90 $(SRCROOT)atmos_param_am3/radiation/clouds/strat_clouds_W.F90 $(SRCROOT)atmos_param_am3/strat_cloud/simple_pdf.F90 $(SRCROOT)atmos_param_am3/vert_diff/vert_diff.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/atmos_tracer_driver.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/endian.F90 $(SRCROOT)atmos_param_am3/cosp/quickbeam/zeff.F90 $(SRCROOT)atmos_param_am3/strat_cloud/micro_mg.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/tropchem/strat_chem_utilities.F90 $(SRCROOT)atmos_shared_am3/tracer_driver/atmos_sea_salt.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/matrix_operations.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/clip_semi_implicit.F90 $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/advance_windm_edsclrm_module.F90 $(SRCROOT)atmos_param_am3/cosp/cosp_defs.H $(SRCROOT)atmos_param_am3/clubb/CLUBB_core/recl.inc $(SRCROOT)fms/include/fms_platform.h $(SRCROOT)atmos_param_am3/donner_deep/donner_types.h $(SRCROOT)atmos_param_am3/monin_obukhov/monin_obukhov_interfaces.h $(SRCROOT)atmos_param_am3/cosp/icarus-scops-3.7/congvec.H $(SRCROOT)atmos_param_am3/strat_cloud/strat_nml.h $(SRCROOT)atmos_shared_am3/atmos_cmip_diag/atmos_cmip_interp.inc $(SRCROOT)atmos_param_am3/donner_deep/donner_nml.h
OBJ = longwave_driver.o esfsw_parameters.o mo_photo.o esfsw_bands.o cosp_rttov_simulator.o radar_simulator.o atmos_soa.o model_flags.o stats_type.o cosp_types.o aerosolrad_driver.o bm_massflux.o aerosolrad_package.o isccp_cloud_types.o atmos_lib.o optical_path.o strat_cloud_utilities.o vert_advection.o llnl_stats.o advance_xm_wpxp_module.o gmres_cache.o cosp_isccp_simulator.o mo_hook.o stable_bl_turb.o donner_rad_k.o grey_radiation.o entrain.o radiative_gases.o topo_drag.o array_index.o atmos_dust.o modis_simulator.o mo_setsox.o variables_diagnostic_module.o betts_miller.o mo_jpl.o cosp.o shallow_conv.o cosp_radar.o mixing_length.o gamma_mg.o saturation.o MG_microp_3D.o monin_obukhov_kernel.o detr_ice_num.o mt95.o bulkphys_rad.o cloudrad_driver.o atmos_sulfate.o atmos_convection_tracer.o MISR_simulator.o mo_setinv.o mo_chemdr.o cosp_driver.o fms_donner.o esfsw_driver.o mrgrnk.o atmos_regional_tracer_driver.o lmd_ipsl_stats.o calendar.o physics_types.o uw_clouds_W.o radiation_driver_diag.o solar_data.o alt_cloud.o diag_integral.o load_mie_table.o sigma_sqd_w_module.o atmos_carbon_aerosol.o shortwave_driver.o CLUBB_3D_var.o conv_utilities_k.o cosp_constants.o parameters_radiation.o parameter_indices.o hydrostatic_module.o cosp_diagnostics.o ice_nucl.o stats_rad_zm.o stats_subs.o mg_const.o cloudrad_diagnostics.o output_grads.o m_tracname.o strat_netcdf.o CLUBB_driver_SCM.o ras.o aerosolrad_types.o moist_conv.o gmres_wrap.o optics_lib.o microphysics.o interpolation.o lin_cloud_microphys.o conv_utilities.o sea_esf_rad.o radar_simulator_types.o morrison_gettelman_microp.o sponge_layer_damping.o cosp_lidar.o cloudrad_package.o cloud_zonal.o donner_meso_k.o dry_adj.o conv_closures.o file_functions.o mean_adv.o atmos_co2.o error_code.o atmos_nudge.o output_netcdf.o strat_chem_model.o advance_helper_module.o cloud_rad.o math_lib.o longwave_tables.o physics_radiation_exch.o donner_lite_k.o atmos_ch3i.o gases.o monin_obukhov.o nonfms_donner.o input_names.o edt.o pos_definite_module.o strat_cloud_legacy.o polysvp.o donner_cloud_model_k.o ozone.o deep_conv.o lapack_wrap.o atmos_age_tracer.o conv_plumes.o donner_types.o cloudrad_types.o cloud_spec.o rad_output_file.o mono_flux_limiter.o clip_explicit.o strat_cloud.o mo_imp_slv.o mo_chem_utls.o cosp_modis_simulator.o longwave_clouds.o conv_plumes_k.o stats_variables.o format_input.o qe_moist_convection.o atmos_sulfur_hex.o stats_zm.o pdf_parameter_module.o anl_erf.o diag_cloud.o donner_cape_k.o moz.mods.o donner_deep_miz.o parameters_microphys.o tke_turb.o extrapolation.o cosp_stats.o longwave_utilities.o stat_file_module.o pdf_closure_module.o longwave_params.o donner_deep_k.o mo_fastjx.o stats_rad_zt.o variables_prognostic_module.o donner_deep_clouds_W.o esfsw_utilities.o surface_varnce_module.o array_lib.o vert_turb_driver.o lw_gases_stdtf.o stats_zt.o aer_ccn_act.o moz.mat.o input_reader.o vert_diff_driver.o cosp_simulator.o moist_processes.o clubb_precision.o diffusion.o radiation_diag.o rotstayn_klein_mp.o aerosol_params.o pf_to_mr.o clubb_core.o mo_usrrxt.o microphys_rad.o longwave_fluxes.o lscale_cond.o betaDistribution.o scops.o moistproc_kernels.o moz.subs.o mo_read_sim_chm.o mo_fphoto.o parameters_tunable.o rh_clouds.o aerosol.o T_in_K_module.o variables_radiation_module.o moist_processes_utils.o advance_sclrm_Nd_module.o cloud_obs.o rad_utilities.o diffusivity.o cloud_generator.o load_hydrometeor_classes.o aerosol_types.o cumulus_closure_k.o cosp_utils.o sealw99.o isccp_clouds.o clouds.o my25_turb.o mo_chemini.o atmos_radon.o prec_scops.o damping_driver.o cu_mo_trans.o grid_class.o tropchem_driver.o stats_sfc.o mo_rodas_slv.o csr_matrix_class_3array.o fill_holes.o gas_tf.o donner_utilities_k.o parameters_model.o donner_deep.o radiation_driver.o aerosol_cloud.o aer_in_act.o cosp_io.o uw_conv.o corr_matrix_module.o bm_omp.o aer_ccn_act_k.o diag_cloud_rad.o constants_clubb.o get_random_number_stream.o shallow_physics.o check_nan.o hyper_diffusion_4th_ord.o donner_lscloud_k.o atmos_tracer_utilities.o physics_driver.o mg_drag.o icarus.o wet_deposition_0D.o cosp_misr_simulator.o two_stream_gray_rad.o dsd.o mo_exp_slv.o atmos_cmip_diag.o advance_xp2_xpyp_module.o lidar_simulator.o numerical_check.o strat_chem_driver.o cg_drag.o nc_cond.o radiation_types.o stats_LH_zt.o Skw_module.o cldwat2m_micro.o advance_wp2_wp3_module.o strat_clouds_W.o simple_pdf.o vert_diff.o atmos_tracer_driver.o endian.o zeff.o micro_mg.o strat_chem_utilities.o atmos_sea_salt.o matrix_operations.o clip_semi_implicit.o advance_windm_edsclrm_module.o

clean: neat
	$(RM) $(OBJ) libatmos_phys.a
	$(RM) *.mod

neat:
	$(RM) $(TMPFILES)

TAGS: $(SRC)
	etags $(SRC)

tags: $(SRC)
	ctags $(SRC)

libatmos_phys.a: $(OBJ)
	$(AR) $(ARFLAGS) libatmos_phys.a $(OBJ)
