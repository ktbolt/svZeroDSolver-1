
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
    ~SolverInterface() { 
      std::cout << "########## ~SolverInterface ##########" << std::endl; 
      std::cout << "~SolverInterface" << std::endl; 
    };

    static int problem_id_count;
    static std::map<int,SolverInterface*> interface_list; 

    int problem_id;
    std::string input_file_name_;

    double time_step_size_;
    int num_time_steps_;
    double absolute_tolerance_;
    int max_nliter_;

    MODEL::Model<T> model;
    //ALGEBRA::Integrator<double,S>* integrator;

    int time_step_;
    std::vector<ALGEBRA::State<double>> states;
    std::vector<double> times;
    int save_interval_counter_;
    int output_interval_;

    bool first_time_step_ = true;

};

