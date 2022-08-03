// Copyright (c) Stanford University, The Regents of the University of
//               California, and others.
//
// All Rights Reserved.
//
// See Copyright-SimVascular.txt for additional details.
//
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject
// to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
// OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/**
 * @file interface.cpp
 * @brief svZeroDSolver callinle interace.
 */
#include "interface.h"

typedef double T;

template <typename TT>
using S = ALGEBRA::SparseSystem<TT>;

// Static member data.
//
int SolverInterface::problem_id_count = 0;
std::map<int,SolverInterface*> SolverInterface::interface_list;

//-----------------
// SolverInterface
//-----------------
//
SolverInterface::SolverInterface(const std::string& input_file_name) : input_file_name_(input_file_name)
{
  problem_id = problem_id_count++;
  SolverInterface::interface_list[problem_id] = this;
}

//////////////////////////////////////////////////////////
//            Callible interface functions              //
//////////////////////////////////////////////////////////

extern "C" void initialize(const char* input_file, int& problem_id);
extern "C" void increment_time(const int problem_id, const double time);

//------------
// initialize
//------------
//
void initialize(const char* input_file_arg, int& problem_id)
{
  std::cout << "========== svzero initialize ==========" << std::endl;
  std::string input_file(input_file_arg);
  std::cout << "[initialize] input_file: " << input_file << std::endl;
  std::string output_file = "svzerod.json";

  auto interface = new SolverInterface(input_file);
  problem_id = interface->problem_id;
  std::cout << "[initialize] problem_id: " << problem_id << std::endl;

  // Create configuration reader.
  IO::ConfigReader<T> config(input_file);

  // Create a model.
  auto model = config.get_model();
  //interface->model = config.get_model();
  std::cout << "[initialize] System size: " << model.dofhandler.size() << std::endl;

  // Get simulation parameters.
  //
  T time_step_size = config.get_time_step_size();
  int num_time_steps = config.get_num_time_steps();
  T absolute_tolerance = config.get_scalar_simulation_parameter("absolute_tolerance", 1e-8);
  int max_nliter = config.get_int_simulation_parameter("maximum_nonlinear_iterations", 30);
  int output_interval = config.get_int_simulation_parameter("output_interval", 1);
  bool steady_initial = config.get_bool_simulation_parameter("steady_initial", true);
  bool output_mean_only = config.get_bool_simulation_parameter("output_mean_only", false);

  // Setup system
  ALGEBRA::State<T> state = ALGEBRA::State<T>::Zero(model.dofhandler.size());

  // Create steady initial
  //
  if (steady_initial) {
    std::cout << "[initialize] " << std::endl;
    std::cout << "[initialize] ----- Calculating steady initial condition ----- " << std::endl;
    T time_step_size_steady = config.cardiac_cycle_period / 10.0;
    std::cout << "[initialize] Create steady model ... " << std::endl;
    auto model_steady = config.get_model();

    std::cout << "[initialize] model_steady.to_steady() ... " << std::endl;
    model_steady.to_steady();

    std::cout << "[initialize] integrator_steady ... " << std::endl;
    ALGEBRA::Integrator<T, S> integrator_steady(model_steady, time_step_size_steady, 0.1, absolute_tolerance, max_nliter);

    for (size_t i = 0; i < 31; i++) {
      state = integrator_steady.step(state, time_step_size_steady * T(i), model_steady);
    }
  }
  std::cout << "[initialize] ... done Calculating steady initial condition" << std::endl;

  #define n_do_original_steppping
  #ifdef do_original_steppping
  std::cout << "[initialize] ----- do_original_steppping ----- " << std::endl;
  ALGEBRA::Integrator<T, S> integrator(model, time_step_size, 0.1, absolute_tolerance, max_nliter);
  std::vector<ALGEBRA::State<T>> states;
  std::vector<T> times;
  states.reserve(num_time_steps);
  times.reserve(num_time_steps);

  T time = 0.0;

  states.push_back(state);
  times.push_back(time);

  int interval_counter = 0;
  for (size_t i = 1; i < num_time_steps; i++) {
    state = integrator.step(state, time, interface->model);
    interval_counter += 1;
    time = time_step_size * T(i);
    /*
    std::cout << "[initialize] " << std::endl;
    std::cout << "[initialize] time: " << time << std::endl;
    std::cout << "[initialize] state.y[0]: " << state.y[0] << std::endl;
    std::cout << "[initialize] state.y[1]: " << state.y[1] << std::endl;
    std::cout << "[initialize] state.y[2]: " << state.y[2] << std::endl;
    */

    if (interval_counter == output_interval) {
      times.push_back(time);
      states.push_back(std::move(state));
      interval_counter = 0;
    }
  }

  #endif

  //interface->model = config.get_model();
  //interface->model = model;

  interface->num_time_steps_ = num_time_steps;
  interface->time_step_size_ = time_step_size;
  interface->max_nliter_ = max_nliter;
  interface->absolute_tolerance_ = absolute_tolerance;

  interface->states.reserve(num_time_steps);
  interface->states.push_back(state);

  interface->times.reserve(num_time_steps);
  interface->times.push_back(0.0);

  interface->time_step_ = 0;
  interface->save_interval_counter_ = 0;
  interface->output_interval_ = output_interval;

  /*
  std::cout << "[initialize] " << std::endl;
  std::cout << "[initialize] Examine model data ..." << std::endl;

  for (auto &&elem : interface->model.blocks) {
    std::cout << "[initialize] ----- Block " << elem.first << " -----" << std::endl;
    std::visit([&](auto &&block) { 
      std::cout << "[initialize] inlet_nodes[0]: " << block.inlet_nodes.size() << std::endl;
      if ( block.inlet_nodes.size() == 1) { 
        std::cout << "[initialize] inlet_nodes[0]: " << block.inlet_nodes[0] << std::endl;
      }
    }, elem.second);
  }
  */

  std::cout << "[initialize] Done " << std::endl;
}

//----------------
// increment_time
//----------------
//
void increment_time(const int problem_id, const double time)
{
  std::cout << "========== svzero increment_time ==========" << std::endl;
  std::cout << "[increment_time] Step in time " << std::endl;
  std::cout << "[increment_time] time: " << time << std::endl;
  auto interface = SolverInterface::interface_list[problem_id];

  // Create a model.
  if (interface->first_time_step_) {
    IO::ConfigReader<T> config(interface->input_file_name_);
    interface->model = config.get_model();
    interface->first_time_step_ = false;
  }
  IO::ConfigReader<T> config(interface->input_file_name_);
  interface->model = config.get_model();
  //auto model = config.get_model();
  //std::cout << "[increment_time] System size: " << interface->model.dofhandler.size() << std::endl;

  auto time_step_size = interface->time_step_size_;
  auto absolute_tolerance = interface->absolute_tolerance_;
  auto max_nliter = interface->max_nliter_;
  std::cout << "[initialize] time_step_size: " << time_step_size << std::endl;
  std::cout << "[initialize] absolute_tolerance: " << absolute_tolerance << std::endl;
  std::cout << "[initialize] max_nliter: " << max_nliter << std::endl;
  std::cout << "[initialize] System size: " << interface->model.dofhandler.size() << std::endl;

  ALGEBRA::Integrator<T, S> integrator(interface->model, time_step_size, 0.1, absolute_tolerance, max_nliter);
  //ALGEBRA::Integrator<T, S> integrator(interface->model, time_step_size, 0.1, absolute_tolerance, max_nliter);
  //interface->integrator = &integrator; 

  return;


  //std::vector<ALGEBRA::State<T>> states;
  //std::vector<T> times;
  //interface->states.reserve(num_time_steps);
  //interface->times.reserve(num_time_steps);

  //interface->states.push_back(state);
  //interface->times.push_back(time);


  auto state = interface->states.back();
  //state = interface->integrator->step(state, time, interface->model);
  int time_step = interface->time_step_;
  //double time = time_step_size * T(time_step);
  std::cout << "[initialize] time: " << time << std::endl;
  interface->time_step_ += 1;

  if (interface->save_interval_counter_ == interface->output_interval_) {
    interface->times.push_back(time);
    interface->states.push_back(std::move(state));
    interface->save_interval_counter_ = 0;
  }

/*
  std::string output;

  if (HELPERS::endswith(output_file, ".csv")) {
    DEBUG_MSG("Saving csv result file to " << output_file);
    output = IO::write_csv<T>( interface->times,  interface->states,  interface->model, output_mean_only);
  } else if (HELPERS::endswith(output_file, ".json")) {
    DEBUG_MSG("Saving json result file to " << output_file);
    output = IO::write_json<T>( interface->times,  interface->states,  interface->model);
  } else {
    throw std::runtime_error("Unsupported outfile file format.");
  }

  std::ofstream ofs(output_file);
  ofs << output;
  ofs.close();
*/
}
