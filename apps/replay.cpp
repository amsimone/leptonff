#include "leptonff/production_records.hpp"
#include "leptonff/decay_setup.hpp"
#include "../src/parallel.hpp"
#include <filesystem>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>

namespace {
using namespace leptonff;
std::string quoted(const std::string& value) {
  std::ostringstream out;out<<'"';
  for(unsigned char c:value) {
    if(c=='"' || c=='\\')out<<'\\'<<char(c);
    else if(c<32)out<<"\\u"<<std::hex<<std::setw(4)<<std::setfill('0')<<unsigned(c)<<std::dec;
    else out<<char(c);
  }
  return out.str()+'"';
}
DecayConfiguration configuration(SpinDecay kind,const std::string& model,Real alpha,std::size_t bins,unsigned hard,Real scale,Real soft,Real collinear) {
  auto c=make_decay_configuration(kind,decay_preset(model),AlphaInput::on_shell(alpha),.1056583755L,scale);
  c.grid.bins=bins;c.hard_samples=hard;c.soft_cutoff=soft;c.collinear_cutoff=collinear;
  return c;
}
}
int main(int argc,char** argv) {
  try {
    using namespace leptonff;
    if(argc==2 && std::string(argv[1])=="--help") {
      std::cout<<"Usage: leptonff_replay production.lff --output report.json [--model born|photonic_ll|photonic_nll|leptonic_nll|qed_hvp] [--subdivide 4|--binning edges.conf] [--grid-bins N] [--hard-samples N] [--threads N] [--batch-groups N] [--covariance full|per-weight|diagonal]\n"
        <<"Repeat --model for paired predictions and full cross-model covariance. Use leptonff-replay manifest.json for checksum-verified replay.\n";return 0;
    }
    if(argc<2)throw std::invalid_argument("Usage: leptonff_replay production.lff --output report.json [--model born|photonic_ll|photonic_nll|leptonic_nll|qed_hvp] [--subdivide 4|--binning edges.conf] [--grid-bins N] [--hard-samples N] [--threads N] [--batch-groups N] [--covariance full|per-weight|diagonal]");
    std::string input=argv[1],output,binning;std::vector<std::string> models;
    std::size_t grid=4096,batch_groups=128;unsigned hard=16,subdivide=1,threads=1;std::string covariance="full";
    Real alpha=1/137.035999084L,scale=0,soft=1e-5L,collinear=1e-5L;
    WLabCuts w;FiducialCuts z;
    for(int i=2;i<argc;++i) {
      const std::string key=argv[i];if(i+1>=argc)throw std::invalid_argument("Missing option value: "+key);const std::string value=argv[++i];
      if(key=="--output")output=value;else if(key=="--model")models.push_back(value);else if(key=="--binning")binning=value;
      else if(key=="--threads") {const auto n=std::stoul(value);if(!n || n>1024)throw std::invalid_argument("Invalid worker count");threads=unsigned(n);}
      else if(key=="--batch-groups")batch_groups=std::stoull(value);
      else if(key=="--covariance")covariance=value;
      else if(key=="--subdivide")subdivide=std::stoul(value);else if(key=="--grid-bins")grid=std::stoull(value);
      else if(key=="--hard-samples")hard=std::stoul(value);else if(key=="--alpha")alpha=std::stold(value);
      else if(key=="--mu-frag")scale=std::stold(value);else if(key=="--soft-cutoff")soft=std::stold(value);
      else if(key=="--collinear-cutoff")collinear=std::stold(value);
      else if(key=="--pt-min")w.muon_pt_min=z.pt_min=std::stold(value);else if(key=="--eta-max")w.muon_eta_max=z.eta_max=std::stold(value);
      else if(key=="--mass-min")z.mass_min=std::stold(value);else if(key=="--mass-max")z.mass_max=std::stold(value);
      else if(key=="--nu-pt-min")w.neutrino_pt_min=std::stold(value);else if(key=="--mt-min")w.transverse_mass_min=std::stold(value);
      else throw std::invalid_argument("Unknown option: "+key);
    }
    if(output.empty() || std::filesystem::exists(output))throw std::invalid_argument("A new --output path is required");
    if(!binning.empty() && subdivide!=1)throw std::invalid_argument("Choose explicit binning or subdivision");
    if(models.empty())models={"photonic_nll"};
    if(std::set<std::string>(models.begin(),models.end()).size()!=models.size())throw std::invalid_argument("Duplicate replay model");
    ProductionReader reader(input);const auto& meta=reader.metadata();const bool neutral=meta.process==SpinDecay::neutral_vector;
    if(scale==0)scale=neutral?91.1876L:80.385L;
    const auto spectrum=LabSpectrum::configured(neutral,binning.empty()?LabSpectrum::reference_binning(neutral,subdivide):LabSpectrum::read_binning(binning),w,z);
    std::vector<DecayCalculator> calculators;
    for(const auto& model:models)calculators.emplace_back(configuration(meta.process,model,alpha,grid,hard,scale,soft,collinear));
    const std::size_t channels=meta.weight_names.size()*models.size()*(1+spectrum.bins.size())*5;
    CovarianceRequest request;
    if(covariance=="diagonal")request.mode=CovarianceMode::diagonal;
    else if(covariance=="per-weight") {request.mode=CovarianceMode::blocks;request.block_size=models.size()*(1+spectrum.bins.size())*5;}
    else if(covariance!="full")throw std::invalid_argument("Covariance must be full, per-weight or diagonal");
    if(!batch_groups || batch_groups>100000)throw std::invalid_argument("Invalid replay batch size");
    batch_groups=std::min(batch_groups,std::max(std::size_t{1},(64*1024*1024)/(channels*sizeof(Real))));
    const MeasurementPlan measurement(neutral,spectrum);
    GroupedIntegration integral(channels,request);ProductionBlock block;std::size_t groups=0,terms=0,blocks=0;
    std::vector<ProductionBlock> block_metadata;
    std::vector<ProductionGroup> pending;pending.reserve(batch_groups);
    const auto flush=[&]{
      std::vector<std::vector<Real>> values(pending.size(),std::vector<Real>(channels));
      detail::parallel_for(pending.size(),threads,[&](std::size_t i){
        for(const auto& t:pending[i].terms)for(std::size_t model=0;model<models.size();++model) {
          const auto result=calculators[model].evaluate(t.point,measurement);
          for(std::size_t weight=0;weight<meta.weight_names.size();++weight)for(std::size_t bin=0;bin<=spectrum.bins.size();++bin) {
            const auto p=bin?result.bins[bin-1].values():result.total.values();
            const auto offset=(((weight*models.size()+model)*(1+spectrum.bins.size()))+bin)*5;
            for(int j=0;j<5;++j)values[i][offset+j]+=t.weights[weight]*p[j];
          }
        }
      });
      // Stable input order, independent of worker scheduling. Complete signed
      // groups are reduced before products; no splitting across host blocks.
      for(const auto& group:values)integral.add_group(group);
      pending.clear();
    };
    while(reader.next_block(block,[&](const ProductionGroup& group) {
      ++groups;terms+=group.terms.size();pending.push_back(group);
      if(pending.size()==batch_groups)flush();
    })) {flush();integral.end_block(block);block_metadata.push_back(block);++blocks;}
    const auto result=integral.estimate();const Real factor=meta.units=="fb"?1e-3L:1;
    std::ofstream out(output);if(!out)throw std::runtime_error("Cannot create replay output");out<<std::setprecision(18);
    out<<"{\n\"schema\":"<<(covariance=="full"?1:2)<<",\"scope\":\"grouped pre-FSR replay; integration covariance, no theory uncertainty claim\",\n\"process\":"
      <<quoted(neutral?"z":meta.process==SpinDecay::wplus?"wplus":"wminus")<<",\"input\":"<<quoted(input)
      <<",\"producer\":"<<quoted(meta.producer)<<",\"physics_id\":"<<quoted(meta.physics_id)<<",\"run_id\":"<<quoted(meta.run_id)
      <<",\"random_seed\":"<<meta.random_seed<<",\"groups\":"<<groups<<",\"terms\":"<<terms<<",\"blocks\":"<<blocks<<",\n\"configuration\":{\"alpha_input\":"<<alpha
      <<",\"mu_frag\":"<<scale<<",\"grid_bins\":"<<grid<<",\"hard_samples\":"<<hard<<",\"soft_cutoff\":"<<soft<<",\"collinear_cutoff\":"<<collinear
      <<",\"pt_min\":"<<w.muon_pt_min<<",\"eta_max\":"<<w.muon_eta_max<<",\"neutrino_pt_min\":"<<w.neutrino_pt_min<<",\"transverse_mass_min\":"<<w.transverse_mass_min
      <<",\"mass_min\":"<<z.mass_min<<",\"mass_max\":"<<z.mass_max<<"},\n\"models\":[";
    for(std::size_t i=0;i<models.size();++i)out<<(i?",":"")<<quoted(models[i]);
    out<<"],\n\"model_configuration\":[";
    for(std::size_t i=0;i<models.size();++i) {
      const auto& c=calculators[i].configuration();if(i)out<<',';
      out<<"{\"name\":"<<quoted(models[i])<<",\"alpha_msbar_mu0\":"<<c.evolution.alpha_msbar_at_mu0
        <<",\"muon_mass\":"<<c.evolution.muon_mass<<",\"mu0\":"<<c.evolution.mu_0<<",\"log_span\":"<<c.grid.log_span
        <<",\"minimum_z\":"<<c.minimum_z<<",\"fraction_measure\":"<<quoted(c.fraction_measure==FractionMeasure::atoms?"atoms":"conservative_cells")
        <<",\"coupled_flavours\":"<<(c.evolution.flavour_scheme?"true":"false")<<",\"pointlike_hvp\":"<<(c.evolution.hadronic_matching?"true":"false");
      const auto* p=calculators[i].photonic();
      if(p && p->nll_evolution()) {
        const auto& d=p->nll_evolution()->diagnostics();
        out<<",\"normalisation\":\"unitary-conditioned\",\"source_number_moment\":"<<d.source_number_moment
          <<",\"raw_inclusive_factor\":"<<d.raw_inclusive_factor<<",\"wrap_bound\":"<<d.wrap_bound<<",\"tail_bound\":"<<d.omitted_tail_bound;
      }
      out<<'}';
    }
    out<<"],\n\"integration_blocks\":[";
    for(std::size_t i=0;i<block_metadata.size();++i) {
      const auto& b=block_metadata[i];if(i)out<<',';
      out<<"{\"part\":"<<b.part<<",\"channel\":"<<b.channel<<",\"iteration\":"<<b.iteration<<",\"trials\":"<<b.trials<<",\"stored_groups\":"<<b.stored_groups<<'}';
    }
    out<<"],\n\"bins\":[";
    for(std::size_t k=0;k<spectrum.bins.size();++k) {
      const auto& b=spectrum.bins[k];if(k)out<<',';out<<"{\"observable\":"<<quoted(b.observable)<<",\"lower\":"<<b.lower<<",\"upper\":";
      if(std::isfinite(b.upper))out<<b.upper;else out<<"null";out<<'}';
    }
    const char* names[]={"born","order1","matched","order2","higher"};out<<"],\n\"labels\":[";bool first=true;
    for(const auto& weight:meta.weight_names)for(const auto& model:models)for(std::size_t k=0;k<=spectrum.bins.size();++k)for(int j=0;j<5;++j) {
      if(!first)out<<',';first=false;const std::string bin=k?spectrum.bins[k-1].observable+":"+std::to_string(k-1):"total";
      out<<quoted(weight+"/"+model+"/"+bin+"/"+names[j]);
    }
    out<<"],\n\"mean_pb\":[";for(std::size_t i=0;i<channels;++i)out<<(i?",":"")<<factor*result.mean[i];
    out<<"],\n\"covariance_mode\":"<<quoted(covariance)<<",\"worker_threads\":"<<threads<<",\"batch_groups\":"<<batch_groups;
    out<<",\n\"variance_pb2\":[";for(std::size_t i=0;i<channels;++i)out<<(i?",":"")<<factor*factor*result.variance[i];
    if(covariance!="full") {
      out<<"],\n\"covariance_pairs\":[";
      for(std::size_t k=0;k<result.covariance_pairs.size();++k) {if(k)out<<',';auto pair=result.covariance_pairs[k];out<<'['<<pair.first<<','<<pair.second<<']';}
      out<<"],\n\"covariance_values_pb2\":[";
      for(std::size_t k=0;k<result.covariance_values.size();++k)out<<(k?",":"")<<factor*factor*result.covariance_values[k];
    }else {
    out<<"],\n\"covariance_pb2\":[";
    for(std::size_t i=0;i<channels;++i) {if(i)out<<',';out<<'[';for(std::size_t j=0;j<channels;++j)out<<(j?",":"")<<factor*factor*result.covariance[i][j];out<<']';}
    }
    out<<"]\n}\n";if(!out)throw std::runtime_error("Replay output write failed");
    std::cout<<"Replayed "<<groups<<" groups / "<<terms<<" terms into "<<spectrum.bins.size()<<" bins for "<<models.size()<<" models\n";
  }catch(const std::exception& e){std::cerr<<"leptonff_replay: "<<e.what()<<'\n';return 1;}
}
