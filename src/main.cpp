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
 * @file main.cpp
 * @brief Main routine of svZeroDSolver
 */
#include "algebra/densesystem.hpp"
#include "algebra/integrator.hpp"
#include "algebra/sparsesystem.hpp"
#include "algebra/state.hpp"
#include "helpers/debug.hpp"
#include "helpers/endswith.hpp"
#include "io/configreader.hpp"
#include "io/csvwriter.hpp"
#include "io/jsonwriter.hpp"
#include "model/model.hpp"

//typedef float T;
typedef double T;

template <typename TT>
using S = ALGEBRA::SparseSystem<TT>;

template <typename T>
class MainTest {
  public:
    MODEL::Model<T>* model;
};

/**
 *
 * @brief svZeroDSolver main routine
 *
 * This is the main routine of the svZeroDSolver. It exectutes the following
 * steps:
 *
 * 1. Read the input file
 * 2. Create the 0D model
 * 3. (Optional) Solve for steady initial condition
 * 4. Run simulation
 * 5. Write output to file
 *
 * @param argc Number of command line arguments
 * @param argv Command line arguments
 * @return Return code
 */
int main(int argc, char *argv[]) 
{
  DEBUG_MSG("Starting svZeroDSolver");

  // Get input and output file name
  if (argc != 3) {
    std::cout << "Usage: svzerodsolver path/to/config.json path/to/output.json"
              << std::endl;
    exit(1);
  }
  std::string input_file = argv[1];
  std::string output_file = argv[2];
  DEBUG_MSG("Reading configuration from " << input_file);

  // Create configuration reader
  std::cout << "Create configuration reader ..." << std::endl;
  IO::ConfigReader<T> config(input_file);

  // Create model
  std::cout << std::endl;
  std::cout << "Creating model ..." << std::endl;
  MainTest<T> mtest;
  mtest.model = config.get_model();
  auto model = mtest.model;
  //auto model = config.get_model();
  std::cout << "Size of system: " << model->dofhandler.size() << std::endl;

  // Get simulation parameters
  std::cout << std::endl;
  std::cout << "Setup simulutation ..." << std::endl;
  T time_step_size = config.get_time_step_size();
  DEBUG_MSG("Time step size:      " << time_step_size);
  int num_time_steps = config.get_num_time_steps();
  DEBUG_MSG("Number of timesteps: " << num_time_steps);

  T absolute_tolerance = config.get_scalar_simulation_parameter("absolute_tolerance", 1e-8);
  int max_nliter = config.get_int_simulation_parameter("maximum_nonlinear_iterations", 30);
  int output_interval = config.get_int_simulation_parameter("output_interval", 1);
  bool steady_initial = config.get_bool_simulation_parameter("steady_initial", true);
  bool output_mean_only = config.get_bool_simulation_parameter("output_mean_only", false);

  // Setup system
  std::cout << "Starting simulutation ..." << std::endl;
  ALGEBRA::State<T> state = ALGEBRA::State<T>::Zero(model->dofhandler.size());

  // Create steady initial
  if (steady_initial) {
    std::cout << std::endl;
    std::cout << "------ Calculating steady initial condition ------ " << std::endl;
    T time_step_size_steady = config.cardiac_cycle_period / 10.0;
    auto model_steady = config.get_model();
    model_steady->to_steady();
    ALGEBRA::Integrator<T, S> integrator_steady(*model_steady, time_step_size_steady, 0.1, absolute_tolerance, max_nliter);

    for (size_t i = 0; i < 31; i++) {
      state = integrator_steady.step(state, time_step_size_steady * T(i), *model_steady);
    }
    std::cout << "... done Calculating steady initial condition" << std::endl;
  }
  std::cout << std::endl;
  std::cout << "[main] ----- model inlet_nodes -----" << std::endl;
  for (auto& node : model->nodes) {
    std::cout << "[main] node.flow_dof: " << node.flow_dof << std::endl;
  }

  std::cout << std::endl;
  std::cout << "Create integrator ..." << std::endl;
  ALGEBRA::Integrator<T, S> integrator(*model, time_step_size, 0.1, absolute_tolerance, max_nliter);

  std::cout << std::endl;
  std::cout << "[main] ----- model inlet_nodes -----" << std::endl;
  for (auto& node : model->nodes) {
    std::cout << "[main] node.flow_dof: " << node.flow_dof << std::endl;
  }

  std::vector<ALGEBRA::State<T>> states;
  std::vector<T> times;
  states.reserve(num_time_steps);
  times.reserve(num_time_steps);

  T time = 0.0;

  states.push_back(state);
  times.push_back(time);

  std::cout << std::endl;
  std::cout << "------ Do simulation ------ " << std::endl;
  int interval_counter = 0;
  for (size_t i = 1; i < num_time_steps; i++) {
    state = integrator.step(state, time, *model);
    interval_counter += 1;
    time = time_step_size * T(i);
    /*
    std::cout << "[main] " << std::endl;
    std::cout << "[main] time: " << time << std::endl;
    std::cout << "[main] state.y[0]: " << state.y[0] << std::endl;
    std::cout << "[main] state.y[1]: " << state.y[1] << std::endl;
    std::cout << "[main] state.y[2]: " << state.y[2] << std::endl;
    */

    if (interval_counter == output_interval) {
      times.push_back(time);
      states.push_back(std::move(state));
      interval_counter = 0;
    }
  }
  DEBUG_MSG("Simulation completed");

  std::string output;
  if (HELPERS::endswith(output_file, ".csv")) {
    DEBUG_MSG("Saving csv result file to " << output_file);
    output = IO::write_csv<T>(times, states, *model, output_mean_only);
  } else if (HELPERS::endswith(output_file, ".json")) {
    DEBUG_MSG("Saving json result file to " << output_file);
    output = IO::write_json<T>(times, states, *model);
  } else {
    throw std::runtime_error("Unsupported outfile file format.");
  }

  std::ofstream ofs(output_file);
  ofs << output;
  ofs.close();

  return 0;
}
