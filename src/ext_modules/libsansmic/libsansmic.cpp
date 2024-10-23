// © 2024 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
// Government retains certain rights in this software.

// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file libsansmic.cpp
 * @brief Python bindings for sansmic.libsansmic
 */

#include "libsansmic.hpp"

#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(libsansmic, m) {
  m.doc() = "SANSMIC model plugin";

  py::enum_<sansmic::LeachMode>(m, "CRunMode",
                                "See :py:class:`sansmic.model.SimulationMode`.")
      .value("ORDINARY", sansmic::LeachMode::Ordinary,
             "ordinary leaching (water in, brine out)")
      .value("WITHDRAWAL", sansmic::LeachMode::Withdrawal,
             "oil withdrawal (water in, oil out)")
      .value("LEACH_FILL", sansmic::LeachMode::LeachFill,
             "leach and fill (water and oil in, brine out)")
      .value("OIL_FILL", sansmic::LeachMode::OilFill,
             "oil fill (oil in, brine out)")
      .export_values();

  py::enum_<sansmic::GeomFormat>(
      m, "CGeometryFormat", "See :py:class:`sansmic.model.GeometryFormat`.")
      .value("RADIUS_LIST", sansmic::GeomFormat::RadiusList, "radius list")
      .value("VOLUME_LIST", sansmic::GeomFormat::VolumeList,
             "volumes follow depths")
      .value("VOLUME_TABLE", sansmic::GeomFormat::VolumeTable,
             "volume-depth pairs")
      .value("RADIUS_TABLE", sansmic::GeomFormat::RadiusTable,
             "radius-depth pairs")
      .export_values();

  py::class_<sansmic::Scenario>(m, "CScenario",
                                "See :py:class:`sansmic.model.Scenario`")
      .def(py::init<>())
      .def("add_stage", &sansmic::Scenario::add_stage, py::arg("stage"),
           "int : Add `stage` to the model and return new number of stages.")
      .def_readwrite("title", &sansmic::Scenario::title, "str : scenario title")
      .def_readwrite("num_cells", &sansmic::Scenario::num_cells,
                     "int : number of cells")
      .def_readwrite("cavern_height", &sansmic::Scenario::cavern_height,
                     "float : height of cavern in feet")
      .def_readwrite("ullage_standoff", &sansmic::Scenario::ullage_standoff,
                     "float : ullage standoff distance in feet")
      .def_readwrite("floor_depth", &sansmic::Scenario::floor_depth,
                     "float : floor depth as in feet measured depth")
      .def_readwrite("geometry_format", &sansmic::Scenario::geometry_format,
                     "int : cavern geometry format")
      .def_readwrite("fraction_insolubles",
                     &sansmic::Scenario::fraction_insolubles,
                     "float : volume fraction of insolubles in salt")
      .def_readwrite("jet_model_version", &sansmic::Scenario::jet_model_version,
                     "int : Jet model version, default is 1")
      .def_readwrite("plume_model_version",
                     &sansmic::Scenario::plume_model_version,
                     "int : Plume model version, default is 1")
      .def_readwrite("temperature_model_version",
                     &sansmic::Scenario::temperature_model_version,
                     "int : Temperature model version, default is 0 (off)")
      .def_readwrite(
          "entrainment_coeff", &sansmic::Scenario::entrainment_coeff,
          "float : Dissolution entrainment coefficient, default is 0.09")
      .def_readwrite("diffusion_beta", &sansmic::Scenario::diffusion_beta,
                     "float : Diffusion beta coefficient, default is 0.147")
      .def_readwrite(
          "molecular_diffusion", &sansmic::Scenario::molecular_diffusion,
          "float : Molecular diffusion coefficient, default is 5.03e-5")
      .def_readwrite("eddy_coefficient", &sansmic::Scenario::eddy_coefficient,
                     "float : Eddy coefficient, default is 1.142e5")
      .def_readwrite(
          "relative_error", &sansmic::Scenario::relative_error,
          "float : ODE solver absolute error tolerance, default is 0.01")
      .def_readwrite(
          "absolute_error", &sansmic::Scenario::absolute_error,
          "float : ODE solver absolute error tolerance, default is 0.0001")
      .def_readwrite("dissolution_factor",
                     &sansmic::Scenario::dissolution_factor,
                     "float : Dissolution factor, default is 1.0")
      .def_readwrite(
          "max_brine_sg", &sansmic::Scenario::max_brine_sg,
          "float : Maximum brine specific gravity, default is 1.2019 (NaCl)")
      .def_readwrite(
          "solid_density", &sansmic::Scenario::solid_density,
          "float : Salt (NaCl) density in solid form, default is 2.16 g/cc")
      .def_readwrite("coallescing_wells", &sansmic::Scenario::coallescing_wells,
                     "int : Number of coallescing wells, default is 1")
      .def_readwrite(
          "well_separation", &sansmic::Scenario::well_separation,
          "float : Separation between coallescing wells, default is 0.0 feet")

      .def_readwrite("radius_vector", &sansmic::Scenario::geom_radii,
                     "list[float] : radius at each node, len="
                     "(num_cells+1)")
      .def_readwrite("depth_vector", &sansmic::Scenario::geom_depths,
                     "list[float] : depth at each node, len=(num_cells+1)")
      .def_readwrite("volume_vector", &sansmic::Scenario::geom_volumes,
                     "list[float] : volume at each node, len=(num_cells+1)");

  py::class_<sansmic::Stage>(m, "CStage",
                             "See :py:class:`sansmic.model.StageDefinition`")
      .def(py::init<>())
      .def_readwrite("title", &sansmic::Stage::title, "str : stage title")
      .def_readwrite("mode", &sansmic::Stage::leach_mode,
                     ":class:`RunMode` : leaching mode")
      .def_readwrite("print_interval", &sansmic::Stage::print_freq,
                     "int : steps between report output")
      .def_readwrite("subsequent", &sansmic::Stage::is_subsequent,
                     "int : is this a subsequent stage?")
      .def_readwrite(
          "rest_duration", &sansmic::Stage::t_rest,
          "double : time cavern sits static before next stage (in hours)")
      .def_readwrite("stop_value", &sansmic::Stage::stop_cond_val,
                     "double : stop criteria; >0 = volume; <0 = OBI; 0 = time")
      .def_readwrite("injection_depth", &sansmic::Stage::h_inj,
                     "double : injection height above floor (in feet)")
      .def_readwrite("production_depth", &sansmic::Stage::h_prod,
                     "double : production height above floor (in feet)")
      .def_readwrite(
          "interface_depth", &sansmic::Stage::h_obi,
          "double : oil-brine interface height above floor (in feet)")
      .def_readwrite("injection_rate", &sansmic::Stage::Q_raw,
                     "double : water input rate (in bbl/d)")
      .def_readwrite("inn_tbg_inside_radius", &sansmic::Stage::r_tbgID,
                     "double : inner radius of inner tubing (in inches)")
      .def_readwrite("inn_tbg_outside_radius", &sansmic::Stage::r_tbgOD,
                     "double : outer radius of inner tubing (in inches)")
      .def_readwrite("out_csg_inside_radius", &sansmic::Stage::r_csgID,
                     "double : inner radius of outer tubing (in inches)")
      .def_readwrite("out_csg_outside_radius", &sansmic::Stage::r_csgOD,
                     "double : outer radius of outer tubing (in inches)")
      .def_readwrite("injection_fluid_sg", &sansmic::Stage::sg_raw,
                     "double : injected water specific gravity")
      .def_readwrite("cavern_sg", &sansmic::Stage::sg_init,
                     "double : cavern brine specific gravity")
      .def_readwrite("timestep", &sansmic::Stage::dt,
                     "double : simulation timestep size (in hours)")
      .def_readwrite("injection_duration", &sansmic::Stage::t_stage,
                     "double : injection duration (in hours)")
      .def_readwrite("fill_rate", &sansmic::Stage::Q_oil,
                     "double : oil input rate (in bbl/d)");

  py::class_<sansmic::Salt>(m, "CSalt",
                            "A model for salt density and saturation.")
      .def(py::init<>())
      .def("set_saturated_sg", &sansmic::Salt::set_sg_max, py::arg("sg"),
           "Set the maximum sg of saturated brine")
      .def("get_saturated_sg", &sansmic::Salt::get_sg_max,
           "Get the maximum sg of saturated brine")
      .def("set_solid_density", &sansmic::Salt::set_solid_density,
           py::arg("rho"), "Set the density (g/cc) of the rock salt")
      .def("get_solid_density", &sansmic::Salt::get_solid_density,
           "Get the density (g/cc) of the rock salt")
      .def("get_recession_rate_coeff", &sansmic::Salt::get_recession_rate_coeff,
           "**a** for the recession rate equation")
      .def("get_wt_pct_sg_coeff", &sansmic::Salt::get_sg_wt_pct_convert_coeff,
           "Get coefficients **c** for the wt-pct vs specific gravity equation")
      .def("set_recession_rate_coeff", &sansmic::Salt::set_recession_rate_coeff,
           py::arg("a"),
           "Set coefficients **a** [n=6] for the recession rate equation")
      .def("set_wt_pct_sg_coeff", &sansmic::Salt::set_density_conversion_coeff,
           py::arg("c"),
           "Set coefficients **c** [n=3] for the wt-pct vs specific gravity "
           "equation")
      .def("wt_pct_to_sg", &sansmic::Salt::sg, py::arg("wt"),
           py::arg("T") = 75.0,
           "Convert weight-percent NaCl to specific gravity")
      .def("recession_rate", &sansmic::Salt::recession_rate, py::arg("sg"),
           "Calculate the recession rate in ft/s for wall dissolution")
      .def("sg_to_wt_pct", &sansmic::Salt::wt_pct, py::arg("sg"),
           py::arg("T") = 75.0,
           "Convert specific gravity to weight-percent NaCl");

  py::class_<sansmic::Model>(m, "CModel",
                             "See :py:class:`sansmic.model.Scenario`.")
      .def(py::init<std::string>())
      .def("open_outfiles", &sansmic::Model::open_outfiles, py::arg("append"),
           "Open output files, optionally in `append` mode, by default "
           "append=false.")
      .def("configure", &sansmic::Model::configure, py::arg("scenario"),
           "Initialize the model object.")
      .def("run_sim", &sansmic::Model::run_sim,
           py::call_guard<py::scoped_ostream_redirect,
                          py::scoped_estream_redirect>(),
           "Run complete sansmic simulation.")
      .def("run_next_step", &sansmic::Model::run_step,
           py::call_guard<py::scoped_ostream_redirect,
                          py::scoped_estream_redirect>(),
           "Run the next step in stage `istage` of the model.")
      .def("num_stages", &sansmic::Model::get_num_stages,
           "int : Get the number of model stages.")
      .def("is_running", &sansmic::Model::get_is_running,
           "int: Is the model in the middle of a stage.")
      .def("close_outfiles", &sansmic::Model::close_outfiles,
           "Close the C-library model output files.")
      .def("generate_tst_file", &sansmic::Model::set_use_tstfile,
           py::arg("use_file"), "Generate and output a .TST file")
      .def("generate_out_file", &sansmic::Model::set_use_outfile,
           py::arg("use_file"), "Generate and output an .OUT file")
      .def("set_verbosity_level", &sansmic::Model::set_verbosity_level,
           py::arg("verb"), "Set the console verbosity level")
      .def("_get_results", &sansmic::Model::get_results,
           "CResults : get the full results object")
      .def("_get_current_stage", &sansmic::Model::get_current_stage,
           "int : get the current stage number")
      .def("_get_current_time", &sansmic::Model::get_current_time,
           "float : get the current time")
      .def("_get_current_volume", &sansmic::Model::get_current_volume,
           "float : get the current cavern volume")
      .def("_get_current_state", &sansmic::Model::get_current_state,
           "CResults : get the model status results right now")
      .def("_get_stages", &sansmic::Model::get_stages,
           "list[CStage] : get a list of all the stages");

  py::class_<sansmic::Results>(m, "CResults",
                               "See :class:`sansmic.model.Results`.")
      .def(py::init<>())
      .def_readwrite("r_0", &sansmic::Results::r_0, "list[float] : initial radius (ft)")
      .def_readwrite("z_0", &sansmic::Results::z_0, "list[float] : depth (ft MD)")
      .def_readwrite("h_0", &sansmic::Results::h_0, "list[float] : height (ft above original TD)")
      .def_readwrite("step", &sansmic::Results::step, "list[int] : step number")
      .def_readwrite("stage", &sansmic::Results::stage,"list[int] : stage number")
      .def_readwrite("phase", &sansmic::Results::phase, "list[int] : active/passive phase")
      .def_readwrite("injCell", &sansmic::Results::injCell,"list[int] : injection cell index")
      .def_readwrite("prodCell", &sansmic::Results::prodCell, "list[int] : production cell index")
      .def_readwrite("obiCell", &sansmic::Results::obiCell, "list[int] : obi cell index")
      .def_readwrite("plumCell", &sansmic::Results::plmCell, "list[int] : plume cell index")
      .def_readwrite("t", &sansmic::Results::t, "list[float] : time (h)")
      .def_readwrite("err", &sansmic::Results::err, "list[float] : convergence factor")
      .def_readwrite("z_obi", &sansmic::Results::z_obi, "list[float] : obi depth (ft MD)")
      .def_readwrite("z_inj", &sansmic::Results::z_inj, "list[float] : injection depth (ft MD)")
      .def_readwrite("z_prod", &sansmic::Results::z_prod, "list[float] : production depth (ft MD)")
      .def_readwrite("z_plm", &sansmic::Results::z_plm, "list[float] : plume stagnation depth (ft MD)")
      .def_readwrite("z_insol", &sansmic::Results::z_insol, "list[float] : insoluble-based current TD (ft MD)")
      .def_readwrite("h_insol", &sansmic::Results::h_insol, "list[float] : insoluble material height (ft)")
      .def_readwrite("l_jet", &sansmic::Results::l_jet, "list[float] : jet length (ft)")
      .def_readwrite("r_jet", &sansmic::Results::r_jet, "list[float] : jet radius (ft)")
      .def_readwrite("u_jet", &sansmic::Results::u_jet, "list[float] : jet velocity (ft/s)")
      .def_readwrite("V_injTot", &sansmic::Results::V_injTot, "list[float] :total injection volume (bbl)")
      .def_readwrite("V_fillTot", &sansmic::Results::V_fillTot, "list[float] : total fill volume (bbl)")
      .def_readwrite("V_cavTot", &sansmic::Results::V_cavTot, "list[float] : total cavern volume (bbl)")
      .def_readwrite("V_insolTot", &sansmic::Results::V_insolTot, "list[float] : total insolubles deposited volume (bbl)")
      .def_readwrite("V_insolVent", &sansmic::Results::V_insolVent, "list[float] : total insolubles vented volume (bbl)")
      .def_readwrite("Q_out", &sansmic::Results::Q_out, "list[float] : instantaneous outflow rate (bbl/h)")
      .def_readwrite("sg_out", &sansmic::Results::sg_out, "list[float] : instantaneous outflow concentration (sg)")
      .def_readwrite("sg_cavAve", &sansmic::Results::sg_cavAve, "list[float] : instantaneous average cavern brine concentration (sg)")
      .def_readwrite("dt", &sansmic::Results::dt, "list[float] : solver timestep (h)")
      .def_readwrite("r", &sansmic::Results::r_cav, "list[list[float]] : cavern radius (ft)")
      .def_readwrite("dr_0", &sansmic::Results::dr_cav, "list[list[float]] : change in cavern radius (ft)")
      .def_readwrite("sg", &sansmic::Results::sg, "list[list[float]] : brine concentration (sg)")
      .def_readwrite("theta", &sansmic::Results::theta, "list[list[float]] : wall angle (deg)")
      .def_readwrite("Q_inj", &sansmic::Results::Q_inj, "list[list[float]] : injection rate (bbl/h)")
      .def_readwrite("V", &sansmic::Results::V, "list[list[float]] : cell volume (bbl)")
      .def_readwrite("f_dis", &sansmic::Results::f_dis, "list[list[float]] : dissolution factor (:1)")
      .def_readwrite("f_flag", &sansmic::Results::f_flag, "list[list[float]] : dissolution flag (status)")
      .def_readwrite("xincl", &sansmic::Results::xincl, "list[list[float]] : debug")
      .def_readwrite("amd", &sansmic::Results::amd, "list[list[float]] : debug")
      .def_readwrite("D_coeff", &sansmic::Results::D_coeff, "list[list[float]] : diffusion coefficient")
      .def_readwrite("dC_dz", &sansmic::Results::dC_dz, "list[list[float]] : vertical diffusion rate (sg/ )")
      .def_readwrite("C_old", &sansmic::Results::C_old, "list[list[float]] : previous timestep concentration (sg)")
      .def_readwrite("C_new", &sansmic::Results::C_new, "list[list[float]] : current timestep concentration (sg)")
      .def_readwrite("dC", &sansmic::Results::dC_dt, "list[list[float]] : rate of change in concentration (sg/h)")
      .def_readwrite("dr", &sansmic::Results::dr_dt, "list[list[float]] : rate of change in cavern radius (ft/h)")
      .def_readwrite("C_plm", &sansmic::Results::C_plm, "list[list[float]] : plume concentration (sg)")
      .def_readwrite("u_plm", &sansmic::Results::u_plm, "list[list[float]] : plume velocity (ft/s)")
      .def_readwrite("r_plm", &sansmic::Results::r_plm, "list[list[float]] : plume radius (ft)");
}
