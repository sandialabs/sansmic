// © 2024 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
// Government retains certain rights in this software.

// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file libsansmic.cpp
 * @author See AUTHORS.md file
 * @brief Python bindings for SANSMIC
 * @version 1.0.0
 */

#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "sansmic.hpp"

namespace py = pybind11;

PYBIND11_MODULE (libsansmic, m)
{
  m.doc () = "SANSMIC model plugin";

  py::enum_<sansmic::Mode> (m, "CRunMode",
                            "See :py:class:`sansmic.model.SimulationMode`.")
      .value ("ORDINARY", sansmic::Mode::Ordinary,
              "ordinary leaching (water in, brine out)")
      .value ("WITHDRAWAL", sansmic::Mode::Withdrawal,
              "oil withdrawal (water in, oil out)")
      .value ("LEACH_FILL", sansmic::Mode::LeachFill,
              "leach and fill (water and oil in, brine out)")
      .value ("OIL_FILL", sansmic::Mode::OilFill,
              "oil fill (oil in, brine out)")
      .export_values ();

  py::enum_<sansmic::IData> (m, "CGeometryFormat",
                             "See :py:class:`sansmic.model.GeometryFormat`.")
      .value ("RADIUS_LIST", sansmic::IData::RadiusList, "radius list")
      .value ("VOLUME_LIST", sansmic::IData::VolumeList,
              "volumes follow depths")
      .value ("VOLUME_TABLE", sansmic::IData::VolumeTable,
              "volume-depth pairs")
      .value ("RADIUS_TABLE", sansmic::IData::RadiusTable,
              "radius-depth pairs")
      .export_values ();

  py::class_<sansmic::Stage> (m, "CStage",
                              "See :py:class:`sansmic.model.StageDefinition`")
      .def (py::init<> ())
      .def_readwrite ("title", &sansmic::Stage::title, "str : stage title")
      .def_readwrite ("mode", &sansmic::Stage::leach_mode,
                      ":class:`RunMode` : leaching mode")
      .def_readwrite ("print_interval", &sansmic::Stage::print_freq,
                      "int : steps between report output")
      .def_readwrite ("subsequent", &sansmic::Stage::is_subsequent,
                      "int : is this a subsequent stage?")
      .def_readwrite (
          "rest_duration", &sansmic::Stage::t_rest,
          "double : time cavern sits static before next stage (in hours)")
      .def_readwrite ("num_coallescing", &sansmic::Stage::coalesce_num_wells,
                      "int : number of coallescing wells (typically 1)")
      .def_readwrite ("geometry_format", &sansmic::Stage::data_format,
                      "int : cavern geometry format")
      .def_readwrite (
          "stop_value", &sansmic::Stage::stop_cond_val,
          "double : stop criteria; >0 = volume; <0 = OBI; 0 = time")
      .def_readwrite ("injection_depth", &sansmic::Stage::h_inj,
                      "double : injection height above floor (in feet)")
      .def_readwrite ("production_depth", &sansmic::Stage::h_prod,
                      "double : production height above floor (in feet)")
      .def_readwrite (
          "interface_depth", &sansmic::Stage::h_obi,
          "double : oil-brine interface height above floor (in feet)")
      .def_readwrite ("ullage_standoff", &sansmic::Stage::h_ullSO,
                      "double : ullage reference height above floor (in feet)")
      .def_readwrite ("injection_rate", &sansmic::Stage::Q_raw,
                      "double : water input rate (in bbl/d)")
      .def_readwrite ("inn_tbg_inside_radius", &sansmic::Stage::r_tbgID,
                      "double : inner radius of inner tubing (in inches)")
      .def_readwrite ("inn_tbg_outside_radius", &sansmic::Stage::r_tbgOD,
                      "double : outer radius of inner tubing (in inches)")
      .def_readwrite ("out_csg_inside_radius", &sansmic::Stage::r_csgID,
                      "double : inner radius of outer tubing (in inches)")
      .def_readwrite ("out_csg_outside_radius", &sansmic::Stage::r_csgOD,
                      "double : outer radius of outer tubing (in inches)")
      .def_readwrite ("injection_fluid_sg", &sansmic::Stage::sg_raw,
                      "double : injected water specific gravity")
      .def_readwrite ("cavern_sg", &sansmic::Stage::sg_init,
                      "double : cavern brine specific gravity")
      .def_readwrite ("timestep", &sansmic::Stage::dt,
                      "double : simulation timestep size (in hours)")
      .def_readwrite ("injection_duration", &sansmic::Stage::t_stage,
                      "double : injection duration (in hours)")
      .def_readwrite ("fill_rate", &sansmic::Stage::Q_oil,
                      "double : oil input rate (in bbl/d)")
      .def_readwrite ("coallescing_well_separation",
                      &sansmic::Stage::coalesce_well_sep,
                      "double : coallescing wells separation (in feet)")
      .def_readwrite ("radius_vector", &sansmic::Stage::radii,
                      "list[num_cells+1] : radius at each node "
                      "(num_cells+1)")
      .def_readwrite ("depth_vector", &sansmic::Stage::depths,
                      "list[num_cells+1] : depth at each node (num_cells+1)")
      .def_readwrite ("volume_vector", &sansmic::Stage::volumes,
                      "list[num_cells+1] : volume at each node (num_cells+1)");

  py::class_<sansmic::Salt> (m, "CSalt",
                             "A model for salt density and saturation.")
      .def (py::init<> ())
      .def ("set_saturated_sg", &sansmic::Salt::set_sg_max, py::arg ("sg"),
            "Set the maximum sg of saturated brine")
      .def ("get_saturated_sg", &sansmic::Salt::get_sg_max,
            "Get the maximum sg of saturated brine")
      .def ("set_solid_density", &sansmic::Salt::set_solid_density,
            py::arg ("rho"), "Set the density (g/cc) of the rock salt")
      .def ("get_solid_density", &sansmic::Salt::get_solid_density,
            "Get the density (g/cc) of the rock salt")
      .def ("get_recession_rate_coeff",
            &sansmic::Salt::get_recession_rate_coeff,
            "**a** for the recession rate equation")
      .def (
          "get_wt_pct_sg_coeff", &sansmic::Salt::get_sg_wt_pct_convert_coeff,
          "Get coefficients **c** for the wt-pct vs specific gravity equation")
      .def ("set_recession_rate_coeff",
            &sansmic::Salt::set_recession_rate_coeff, py::arg ("a"),
            "Set coefficients **a**[n=6] for the recession rate equation")
      .def ("set_wt_pct_sg_coeff",
            &sansmic::Salt::set_density_conversion_coeff, py::arg ("c"),
            "Set coefficients **c**[n=3] for the wt-pct vs specific gravity "
            "equation")
      .def ("wt_pct_to_sg", &sansmic::Salt::get_sg, py::arg ("wt"),
            py::arg ("T") = 75.0,
            "Convert weight-percent NaCl to specific gravity")
      .def ("recession_rate", &sansmic::Salt::get_recession_rate,
            py::arg ("sg"),
            "Calculate the recession rate in ft/s for wall dissolution")
      .def ("sg_to_wt_pct", &sansmic::Salt::get_wt_pct, py::arg ("sg"),
            py::arg ("T") = 75.0,
            "Convert specific gravity to weight-percent NaCl");

  py::class_<sansmic::Model> (m, "CModel",
                              "See :py:class:`sansmic.model.Scenario`.")
      .def (py::init<std::string> ())
      .def ("add_stage", &sansmic::Model::add_stage, py::arg ("stage"),
            "int : Add `stage` to the model and return new number of stages.")
      .def ("open_outfiles", &sansmic::Model::open_outfiles,
            py::arg ("append"),
            "Open output files, optionally in `append` mode, by default "
            "append=false.")
      .def ("initialize", &sansmic::Model::init_model,
            "Initialize the model object.")
      .def ("run_sim", &sansmic::Model::run_sim,
            py::call_guard<py::scoped_ostream_redirect,
                           py::scoped_estream_redirect> (),
            "Run complete sansmic simulation.")
      // .def ("run_next_stage", &sansmic::Model::run_stage,
      //       py::call_guard<py::scoped_ostream_redirect,
      //                      py::scoped_estream_redirect> (),
      //       "Run the stage `istage` of the model.")
      .def ("run_next_step", &sansmic::Model::run_step,
            py::call_guard<py::scoped_ostream_redirect,
                           py::scoped_estream_redirect> (),
            "Run the next step in stage `istage` of the model.")
      .def ("num_stages", &sansmic::Model::num_stages,
            "Get the number of model stages.")
      .def ("is_running", &sansmic::Model::get_running_status,
            "int: Is the model in the middle of a stage.")
      .def ("close_outfiles", &sansmic::Model::close_outfiles,
            "Close the C-library model output files.")
      .def ("_get_results", &sansmic::Model::get_results)
      .def ("_get_current_stage", &sansmic::Model::get_current_stage)
      .def ("_get_current_time", &sansmic::Model::get_current_time)
      .def ("_get_current_volume", &sansmic::Model::get_current_volume)
      .def ("_get_current_state", &sansmic::Model::get_current_state)
      .def ("_set_jet_model_version", &sansmic::Model::set_jet_model_version,
            py::arg ("ver"), "Set the jet model version to `ver`.")
      .def ("_set_plume_model_version",
            &sansmic::Model::set_plume_model_version, py::arg ("ver"),
            "Set the plume model version to `ver`.")
      .def ("_set_temperature_model_version",
            &sansmic::Model::set_temperature_model_version, py::arg ("ver"),
            "Set the temperature model version to `ver`.")
      .def ("_set_num_cells", &sansmic::Model::set_num_cells,
            py::arg ("n_cells"))
      .def ("_set_cavern_height", &sansmic::Model::set_cavern_height,
            py::arg ("height"))
      .def ("_set_floor_depth", &sansmic::Model::set_initial_floor_depth,
            py::arg ("measured_depth"))
      .def ("_set_entrainment_alpha", &sansmic::Model::set_entrainment_coeff,
            py::arg ("alpha"))
      .def ("_set_diffusion_beta", &sansmic::Model::set_diffusion_beta_coeff,
            py::arg ("beta"))
      .def ("_set_diffusion_D_mol",
            &sansmic::Model::set_molecular_diffusion_coeff, py::arg ("D_mol"))
      .def ("_set_diffusion_D_0", &sansmic::Model::set_eddy_coeff,
            py::arg ("D_0"))
      .def ("_set_fraction_insol", &sansmic::Model::set_fraction_insolubles,
            py::arg ("f_insol"))
      .def ("_set_dissolution_factor", &sansmic::Model::set_dissolution_factor,
            py::arg ("f_diss"))
      .def ("_set_solver_abserr",
            &sansmic::Model::set_ODESolver_absolute_tolerance, py::arg ("tol"))
      .def ("_set_solver_relerr",
            &sansmic::Model::set_ODESolver_relative_tolerance, py::arg ("tol"))
      //  .def_readonly ("_salt", &sansmic::Model::salt, "CSalt : salt
      //  properties")
      .def_readwrite ("_geometry_format", &sansmic::Model::geom_data_format,
                      "int : cavern geometry format")
      .def_readwrite ("_radius_vector", &sansmic::Model::geom_radii,
                      "list[num_cells+1] : radius at each node "
                      "(num_cells+1)")
      .def_readwrite ("_depth_vector", &sansmic::Model::geom_depths,
                      "list[num_cells+1] : depth at each node (num_cells+1)")
      .def_readwrite ("_volume_vector", &sansmic::Model::geom_volumes,
                      "list[num_cells+1] : volume at each node (num_cells+1)")
      .def ("_get_stages", &sansmic::Model::get_stages);

  m.def ("get_error_text", &sansmic::error::get_error_text,
         "Convert error code to text", py::arg ("code"));

  py::class_<sansmic::Results> (m, "CResults",
                                "See :class:`sansmic.model.Results`.")
      .def (py::init<> ())
      .def_readwrite ("r_0", &sansmic::Results::r_0)
      .def_readwrite ("z_0", &sansmic::Results::z_0)
      .def_readwrite ("h_0", &sansmic::Results::h_0)
      .def_readwrite ("step", &sansmic::Results::step)
      .def_readwrite ("stage", &sansmic::Results::stage)
      .def_readwrite ("phase", &sansmic::Results::phase)
      .def_readwrite ("injCell", &sansmic::Results::injCell)
      .def_readwrite ("prodCell", &sansmic::Results::prodCell)
      .def_readwrite ("obiCell", &sansmic::Results::obiCell)
      .def_readwrite ("plumCell", &sansmic::Results::plmCell)
      .def_readwrite ("t", &sansmic::Results::t)
      .def_readwrite ("err", &sansmic::Results::err)
      .def_readwrite ("z_obi", &sansmic::Results::z_obi)
      .def_readwrite ("z_inj", &sansmic::Results::z_inj)
      .def_readwrite ("z_prod", &sansmic::Results::z_prod)
      .def_readwrite ("z_plm", &sansmic::Results::z_plm)
      .def_readwrite ("z_insol", &sansmic::Results::z_insol)
      .def_readwrite ("h_insol", &sansmic::Results::h_insol)
      .def_readwrite ("l_jet", &sansmic::Results::l_jet)
      .def_readwrite ("r_jet", &sansmic::Results::r_jet)
      .def_readwrite ("u_jet", &sansmic::Results::u_jet)
      .def_readwrite ("V_injTot", &sansmic::Results::V_injTot)
      .def_readwrite ("V_fillTot", &sansmic::Results::V_fillTot)
      .def_readwrite ("V_cavTot", &sansmic::Results::V_cavTot)
      .def_readwrite ("V_insolTot", &sansmic::Results::V_insolTot)
      .def_readwrite ("V_insolVent", &sansmic::Results::V_insolVent)
      .def_readwrite ("Q_out", &sansmic::Results::Q_out)
      .def_readwrite ("sg_out", &sansmic::Results::sg_out)
      .def_readwrite ("sg_cavAve", &sansmic::Results::sg_cavAve)
      .def_readwrite ("dt", &sansmic::Results::dt)
      .def_readwrite ("r", &sansmic::Results::r_cav)
      .def_readwrite ("dr_0", &sansmic::Results::dr_cav)
      .def_readwrite ("sg", &sansmic::Results::sg)
      .def_readwrite ("theta", &sansmic::Results::theta)
      .def_readwrite ("Q_inj", &sansmic::Results::Q_inj)
      .def_readwrite ("V", &sansmic::Results::V)
      .def_readwrite ("f_dis", &sansmic::Results::f_dis)
      .def_readwrite ("f_flag", &sansmic::Results::f_flag)
      .def_readwrite ("xincl", &sansmic::Results::xincl)
      .def_readwrite ("amd", &sansmic::Results::amd)
      .def_readwrite ("D_coeff", &sansmic::Results::D_coeff)
      .def_readwrite ("dC_dz", &sansmic::Results::dC_dz)
      .def_readwrite ("C_old", &sansmic::Results::C_old)
      .def_readwrite ("C_new", &sansmic::Results::C_new)
      .def_readwrite ("dC", &sansmic::Results::dC_dt)
      .def_readwrite ("dr", &sansmic::Results::dr_dt)
      .def_readwrite ("C_plm", &sansmic::Results::C_plm)
      .def_readwrite ("u_plm", &sansmic::Results::u_plm)
      .def_readwrite ("r_plm", &sansmic::Results::r_plm);
}
