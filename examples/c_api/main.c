#include <leptonff/c_api.h>
#include <stdio.h>
#include <stdlib.h>

int main(void) {
  leptonff_decay_config config;leptonff_decay_config_default(&config);
  config.grid_bins=128;config.hard_samples=4;
  leptonff_decay_context* context=NULL;
  if(leptonff_decay_create(&config,NULL,0,&context)) {
    fprintf(stderr,"%s\n",leptonff_last_error());return 1;
  }
  const double mu[4]={50,40,0,30},nu[4]={50,-40,0,-30};
  const size_t size=leptonff_decay_output_count(context);
  leptonff_prediction* values=calloc(size,sizeof(*values));
  if(!values || leptonff_decay_evaluate(context,mu,nu,42,values,size)) {
    fprintf(stderr,"%s\n",leptonff_last_error());free(values);leptonff_decay_destroy(context);return 1;
  }
  printf("Born acceptance %.12g; matched response %.12g; %zu bins\n",values[0].born,values[0].matched,size-1);
  free(values);leptonff_decay_destroy(context);return 0;
}
