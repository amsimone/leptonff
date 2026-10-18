#include "leptonff/c_api.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

int main(void) {
  leptonff_decay_config c;leptonff_decay_config_default(&c);
  c.alpha=0;c.grid_bins=128;c.hard_samples=2;
  const double edges[]={25,30,35,40,45,60,INFINITY};
  leptonff_observable obs={"pt",edges,7};leptonff_decay_context* context=NULL;
  assert(leptonff_decay_create(&c,&obs,1,&context)==0);
  const size_t n=leptonff_decay_output_count(context);assert(n==7);
  leptonff_decay_diagnostics diagnostics;
  assert(leptonff_decay_get_diagnostics(context,&diagnostics)==0);
  assert(diagnostics.photonic_unitary_conditioning==1 && diagnostics.raw_inclusive_factor==1);
  assert(fabs(diagnostics.source_number_moment+.719311605774509)<1e-12);
  leptonff_prediction* out=calloc(n,sizeof(*out));
  const double mu[]={50,40,0,30},nu[]={50,-40,0,-30};
  assert(leptonff_decay_evaluate(context,mu,nu,42,out,n)==0);
  assert(out[0].born==1 && out[0].matched==1 && out[0].order1==1);
  double sum=0;for(size_t i=1;i<n;i++)sum+=out[i].matched;
  assert(sum==1);
  assert(leptonff_decay_evaluate(context,mu,nu,42,out,n-1)!=0);
  assert(leptonff_last_error()[0]!='\0');
  double* tensors=calloc(90*n,sizeof(double));const double q[]={100,0,0,0};
  assert(leptonff_decay_response(context,q,42,8,8,tensors,90*n)==0);
  for(size_t i=0;i<90*n;i++)assert(isfinite(tensors[i]));
  free(tensors);free(out);leptonff_decay_destroy(context);
  c.abi_version=999;context=NULL;
  assert(leptonff_decay_create(&c,NULL,0,&context)!=0 && context==NULL);
  puts("C ABI lifecycle, custom bins, identity, tensor layout and error containment passed");
  return 0;
}
