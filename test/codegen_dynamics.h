/* This file was automatically generated by CasADi.
   The CasADi copyright holders make no ownership claim of its contents. */
#ifdef __cplusplus
extern "C" {
#endif

#ifndef casadi_real
#define casadi_real double
#endif

int dynamics(const casadi_real** arg, casadi_real** res, int* iw, casadi_real* w, void* mem);
void dynamics_incref(void);
void dynamics_decref(void);
int dynamics_n_in(void);
int dynamics_n_out(void);
const char* dynamics_name_in(int i);
const char* dynamics_name_out(int i);
const int* dynamics_sparsity_in(int i);
const int* dynamics_sparsity_out(int i);
int dynamics_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
#ifdef __cplusplus
} /* extern "C" */
#endif
