#include "timer_wrapper_c.h";
void initialize_timer_(void** timer){ 
  *timer = initialize_timer();
}
void delete_timer_(void** timer){ 
  
  delete_timer(*timer);
}
void timer_start_(void** timer, int* id){ 
  timer_start(*timer, *id);
}
void timer_stop_(void** timer, int* id){ 
  
  timer_stop(*timer, *id);
}
void timer_registrate_(void** timer, int* id, char* label, long length_label){   
  timer_registrate(*timer, *id, label);
}
double timer_get_count_(void** timer, int* id){
  
  return timer_get_count(*timer, *id);
}    
int timer_get_average_(void** timer, int* id){
  
  return timer_get_average(*timer, *id);
}    
