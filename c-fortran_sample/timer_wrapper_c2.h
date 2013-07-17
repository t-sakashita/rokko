//#include <rokko/solver.hpp>
//#include <rokko/grid.hpp>
//#include <rokko/distributed_matrix.hpp>
//#include <rokko/localized_vector.hpp>

#ifdef __cplusplus
namespace rokko {
  
  extern "C"{
#endif

    void* initialize_timer();
    void delete_timer(void*);
    void timer_start(void*, unsigned int);
    void timer_stop(void*, unsigned int);
    void timer_registrate(void*, unsigned int, char*);
    double timer_get_count(void*, unsigned int);
    int timer_get_average(void*, unsigned int);
#ifdef __cplusplus
  }
}
#endif
