#include "leptonff/c_api.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

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
  size_t label_size=0;assert(leptonff_decay_output_label(context,0,NULL,&label_size)==0);
  char* label=malloc(label_size);assert(leptonff_decay_output_label(context,0,label,&label_size)==0 && strcmp(label,"total")==0);free(label);
  size_t json_size=0;assert(leptonff_decay_configuration_json(context,NULL,&json_size)==0);
  char* resolved=malloc(json_size);assert(leptonff_decay_configuration_json(context,resolved,&json_size)==0 && strstr(resolved,"alpha_msbar_mu0"));free(resolved);
  const double first_batch[]={50,40,0,30,50,40,0,30},second_batch[]={50,-40,0,-30,50,-40,0,-30};
  const uint64_t streams[]={42,137};leptonff_prediction* batch=calloc(2*n,sizeof(*batch));
  assert(leptonff_decay_evaluate_batch(context,first_batch,second_batch,streams,2,2,batch,2*n)==0);
  assert(batch[0].matched==out[0].matched && batch[n].matched==out[0].matched);
  assert(leptonff_decay_evaluate_batch(context,first_batch,second_batch,streams,2,2,batch,2*n-1)!=0);
  for(size_t i=0;i<2*n;++i)batch[i].matched=-123.4;
  const double invalid_second[]={50,-40,0,-30,50,40,0,30};
  assert(leptonff_decay_evaluate_batch(context,first_batch,invalid_second,streams,2,2,batch,2*n)!=0);
  for(size_t i=0;i<2*n;++i)assert(batch[i].matched==-123.4);
  free(batch);
  assert(leptonff_decay_evaluate(context,mu,nu,42,out,n-1)!=0);
  assert(leptonff_last_error()[0]!='\0');
  double* tensors=calloc(90*n,sizeof(double));const double q[]={100,0,0,0};
  assert(leptonff_decay_response(context,q,42,8,8,tensors,90*n)==0);
  for(size_t i=0;i<90*n;i++)assert(isfinite(tensors[i]));
  const double coarse_edges[]={25,40,INFINITY};const leptonff_observable coarse={"pt",coarse_edges,3};
  leptonff_decay_context* view=NULL;assert(leptonff_decay_with_observables(context,NULL,&coarse,1,&view)==0);
  assert(leptonff_decay_output_count(view)==3);
  leptonff_decay_config invalid=c;invalid.alpha=.01;
  leptonff_decay_context* bad=NULL;assert(leptonff_decay_with_observables(context,&invalid,&coarse,1,&bad)!=0 && bad==NULL);
  free(tensors);free(out);leptonff_decay_destroy(context);
  leptonff_prediction viewed[3];assert(leptonff_decay_evaluate(view,mu,nu,42,viewed,3)==0 && viewed[0].matched==1);
  leptonff_decay_destroy(view);
  leptonff_decay_theory t;leptonff_decay_theory_default(&t);
  assert(leptonff_decay_create_with_theory(&c,&t,&obs,1,&context)==0);leptonff_decay_destroy(context);
  t.running=LEPTONFF_ALPHA_TWO_LOOP;
  assert(leptonff_decay_create_with_theory(&c,&t,&obs,1,&context)!=0 && context==NULL);
  c.model=LEPTONFF_LEPTONIC_NLL;t.alpha_reference_scale=.2;t.initial_scale=.2;
  assert(leptonff_decay_create_with_theory(&c,&t,&obs,1,&context)==0);leptonff_decay_destroy(context);
  c.abi_version=999;context=NULL;
  assert(leptonff_decay_create(&c,NULL,0,&context)!=0 && context==NULL);
  puts("C ABI lifecycle, custom bins, identity, tensor layout and error containment passed");
  return 0;
}
