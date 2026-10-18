#ifndef LEPTONFF_C_API_H
#define LEPTONFF_C_API_H
#include <stddef.h>
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
#define LEPTONFF_ABI_VERSION 1
#include "leptonff/fragmentation_c.h"
typedef struct leptonff_decay_context leptonff_decay_context;
/* Generator-independent process identifiers (PDG boson IDs). Legacy 1,6,31
 * remain accepted by ABI v1 for existing clients. */
enum {LEPTONFF_PROCESS_W_PLUS=24,LEPTONFF_PROCESS_W_MINUS=-24,LEPTONFF_PROCESS_Z=23};
enum {LEPTONFF_WPLUS=1,LEPTONFF_WMINUS=6,LEPTONFF_Z=31};
enum {LEPTONFF_PHOTONIC_LL=1,LEPTONFF_PHOTONIC_NLL=2,LEPTONFF_LEPTONIC_NLL=3,LEPTONFF_QED_HVP=4};
typedef struct {
  size_t struct_size;
  unsigned abi_version;
  int process,model;
  /* ABI-v1 presets: alpha is fixed in photonic modes, on-shell in the
   * leptonic/HVP presets. Use the explicit C++ theory API for MSbar input. */
  double alpha,muon_mass,mu_frag;
  size_t grid_bins;
  double log_span;
  unsigned hard_samples;
  double soft_cutoff,collinear_cutoff,minimum_z;
  double muon_pt_min,muon_eta_max,neutrino_pt_min,transverse_mass_min;
  double mass_min,mass_max,hadronic_matching_scale;
} leptonff_decay_config;
typedef struct {const char* observable;const double* edges;size_t edge_count;} leptonff_observable;
typedef struct {double born,order1,matched,order2,higher;} leptonff_prediction;
typedef struct {
  unsigned abi_version;
  int photonic_unitary_conditioning;
  double a_initial,source_number_moment,raw_inclusive_factor,wrap_bound,tail_bound;
} leptonff_decay_diagnostics;
/* All momenta below are (energy, px, py, pz), in GeV. first is the charged W
 * lepton or negative Z lepton; second is the neutrino or positive Z lepton.
 * Each subtraction group must share a radiation seed. Weights remain external.
 * Contexts are immutable and may be used concurrently; last_error is per thread.
 * Return 0 on success, nonzero on failure. Exceptions never cross this ABI. */
void leptonff_decay_config_default(leptonff_decay_config*);
int leptonff_decay_create(const leptonff_decay_config*,const leptonff_observable*,size_t,
                         leptonff_decay_context**);
void leptonff_decay_destroy(leptonff_decay_context*);
/* Output index zero is the fiducial rate; following entries are physical bins. */
size_t leptonff_decay_output_count(const leptonff_decay_context*);
/* Photonic NLL is explicitly conditioned to the inclusive branching ratio.
 * The removed source-sector factor is recorded here; it is not a partial width. */
int leptonff_decay_get_diagnostics(const leptonff_decay_context*,leptonff_decay_diagnostics*);
int leptonff_decay_evaluate(const leptonff_decay_context*,const double first[4],
                           const double second[4],uint64_t radiation_seed,
                           leptonff_prediction*,size_t capacity);
/* Direct response tensors, row-major 3x3 complex entries as real/imag pairs.
 * Component order: born, order1, matched, order2, higher; 90 doubles/output.
 * Cartesian axes: rotation-free boost from LAB. Z is charge symmetric only.
 * capacity is a count of doubles, not a count of output bins. */
int leptonff_decay_response(const leptonff_decay_context*,const double parent[4],
                          uint64_t seed,unsigned cosine_nodes,unsigned azimuth_nodes,
                          double* output,size_t capacity);
const char* leptonff_last_error(void);
#ifdef __cplusplus
}
#endif
#endif
