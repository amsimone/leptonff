#ifndef LEPTONFF_FRAGMENTATION_C_H
#define LEPTONFF_FRAGMENTATION_C_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
#define LEPTONFF_FF_ABI_VERSION 1
enum { LEPTONFF_MUON_MINUS=0, LEPTONFF_MUON_PLUS=1, LEPTONFF_PHOTON=2 };
enum { LEPTONFF_FF_ORDER1=1, LEPTONFF_FF_COUPLED_NLL=2 };
enum { LEPTONFF_ALPHA_ON_SHELL=0, LEPTONFF_ALPHA_MSBAR=1 };
enum { LEPTONFF_ALPHA_FIXED=0, LEPTONFF_ALPHA_ONE_LOOP=1, LEPTONFF_ALPHA_TWO_LOOP=2 };
enum { LEPTONFF_ONE_FAMILY=0, LEPTONFF_LEPTON_THRESHOLDS=1 };
typedef struct leptonff_fragmentation_context leptonff_fragmentation_context;
typedef struct {
  size_t struct_size;
  unsigned abi_version;
  int model,alpha_scheme,running,flavours;
  /* On-shell alpha: reference_scale=0. MSbar: reference_scale=initial_scale.
   * Running and thresholds apply to coupled NLL, independently of the scheme. */
  double alpha,alpha_reference_scale,muon_mass,initial_scale,final_scale;
  /* bins=0: moments only for coupled NLL. Otherwise prepare the x-space
   * operator once. Weights are z*D(z) dz, with positive lower bin support. */
  size_t grid_bins;
  double log_span,absolute_tolerance,relative_tolerance;
} leptonff_fragmentation_config;
void leptonff_fragmentation_config_default(leptonff_fragmentation_config*);
int leptonff_fragmentation_create(const leptonff_fragmentation_config*,leptonff_fragmentation_context**);
void leptonff_fragmentation_destroy(leptonff_fragmentation_context*);
int leptonff_fragmentation_moment(const leptonff_fragmentation_context*,int parent,int observed,double N,double*);
/* [lower,upper), with the endpoint included for upper=1. FO uses adaptive
 * distribution integration; coupled NLL uses the prepared signed lattice. */
int leptonff_fragmentation_integral(const leptonff_fragmentation_context*,int parent,int observed,
                                   double lower,double upper,double*);
/* Fixed-order convolution only. Callback must return finite values and must
 * not throw across C. The callback/data are used only during this call. */
typedef double (*leptonff_function)(double z,void* data);
int leptonff_fragmentation_convolve(const leptonff_fragmentation_context*,int parent,int observed,
                                   double x,leptonff_function,void* data,double*);
/* Thread-local; valid until the next fallible API call on this thread.
 * Calls return 0 on success and nonzero on error. Exceptions are contained. */
const char* leptonff_last_error(void);
const char* leptonff_version(void);
#ifdef __cplusplus
}
#endif
#endif
