
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

#include <map>
#include <string>
#include <vector>

typedef double T;
template <typename TT>
using S = ALGEBRA::SparseSystem<TT>;

//-----------------
// SolverInterface
//-----------------
//
class SolverInterface 
{
  public:
    SolverInterface(const std::string& input_file_name);
    ~SolverInterface(); 

    static int problem_id_count_;
    static std::map<int,SolverInterface*> interface_list_; 

    int problem_id_;
    std::string input_file_name_;

    // Solver parameters.
    //
    // These are set via the solver interface.
    //
    double solver_time_step_ = 0.1;

    // 0D solver parameters. 
    //
    // These are read in from the input JSON solver configuration file.
    //
    double time_step_size_;
    int num_time_steps_;
    double absolute_tolerance_;
    int max_nliter_;

    MODEL::Model<T>* model_;
    //ALGEBRA::Integrator<double,S>* integrator;

    int time_step_;
    std::vector<ALGEBRA::State<double>> states_;
    std::vector<double> times_;
    int save_interval_counter_;
    int output_interval_;

};

