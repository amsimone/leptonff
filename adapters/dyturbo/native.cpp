// SPDX-License-Identifier: MIT
// Native DYTurbo client. The separately installed backend retains GPLv3.
#include "dyturbo.h"
#include "settings.h"
#include "coupling.h"
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace {
std::string quoted(const std::string& s) {
  std::ostringstream out;out<<'"';
  for(unsigned char c:s) {
    if(c=='"' || c=='\\')out<<'\\'<<char(c);
    else if(c<32)throw std::runtime_error("Control character in backend metadata");
    else out<<char(c);
  }
  return out.str()+'"';
}
double finite(double x) {
  if(!std::isfinite(x))throw std::runtime_error("Non-finite backend result");
  return x;
}
}
int main(int argc,char** argv) {
  try {
    const auto* filename=std::getenv("LEPTONFF_DYTURBO_RESULT");
    if(!filename || !*filename)throw std::runtime_error("Set LEPTONFF_DYTURBO_RESULT to a new JSON filename");
    if(std::ifstream(filename))throw std::runtime_error("Refusing to overwrite backend results");
    const auto start=std::chrono::steady_clock::now();
    DYTurbo::Init(argc,argv);
    // Results are integrated cross sections in fb, never densities per bin.
    if(opts.ptbinwidth || opts.ybinwidth || opts.mbinwidth || opts.dsdqt2 || opts.edsdp3)
      throw std::runtime_error("Native client requires integrated cross sections");
    if(opts.qed)throw std::runtime_error("Disable DYTurbo QED resummation; this interface benchmarks QCD production");
    DYTurbo::WarmUp();
    std::ostringstream out;out.precision(17);
    out<<"{\"schema\":\"leptonff-dyturbo-native-1\",\"backend_version\":\"1.4.2\","
       <<"\"scope\":\"QCD production; no leptonff decay applied\",\"units\":\"fb\","
       <<"\"process\":"<<opts.nproc<<",\"qcd_order\":"<<opts.order
       <<",\"fixed_order\":"<<(opts.fixedorder?"true":"false")
       <<",\"fiducial\":"<<(opts.makecuts?"true":"false")
       <<",\"pdf\":"<<quoted(opts.LHAPDFset)<<",\"pdf_member\":"<<opts.LHAPDFmember
       <<",\"seed\":"<<opts.rseed<<",\"resolved_orders\":{\"sudakov\":"<<opts.order_sudak
       <<",\"hard\":"<<opts.order_hcoef<<",\"evolution\":"<<opts.order_evol
       <<",\"exponentiation\":"<<opts.order_expc<<",\"counterterm\":"<<opts.order_ct
       <<",\"vjet\":"<<opts.order_vjet<<"},\"dimensions\":{\"born\":"<<opts.intDimBorn
       <<",\"resummed\":"<<opts.intDimRes<<",\"counterterm\":"<<opts.intDimCT
       <<",\"vjet\":"<<opts.intDimVJ<<"},\"physics\":{\"sqrts\":"<<opts.sroot
       <<",\"ew_scheme\":"<<opts.ewscheme<<",\"gf\":"<<coupling::Gf<<",\"alpha_ew\":"<<coupling::aemmz
       <<",\"mw\":"<<coupling::wmass<<",\"mz\":"<<coupling::zmass<<",\"gamma_w\":"<<coupling::wwidth
       <<",\"gamma_z\":"<<coupling::zwidth<<",\"muR_function\":"<<opts.fmuren
       <<",\"muF_function\":"<<opts.fmufac<<",\"muRes_function\":"<<opts.fmures
       <<",\"muR_factor\":"<<opts.kmuren<<",\"muF_factor\":"<<opts.kmufac
       <<",\"muRes_factor\":"<<opts.kmures<<",\"qt_cut\":"<<opts.qtcut<<",\"qt_over_m_cut\":"<<opts.xqtcut
       <<",\"b_prescription\":"<<opts.bprescription<<",\"np_form_factor\":"<<opts.npff
       <<"},\"bins\":[";
    bool first_bin=true;
    for(DYTurbo::BoundIterator bounds;!bounds.IsEnd();++bounds) {
      DYTurbo::SetBounds(bounds);
      if(!first_bin)out<<',';first_bin=false;
      out<<"{\"mass\":["<<bounds.loBound(DYTurbo::b_M)<<','<<bounds.hiBound(DYTurbo::b_M)
         <<"],\"qt\":["<<bounds.loBound(DYTurbo::b_QT)<<','<<bounds.hiBound(DYTurbo::b_QT)
         <<"],\"rapidity\":["<<bounds.loBound(DYTurbo::b_Y)<<','<<bounds.hiBound(DYTurbo::b_Y)
         <<"],\"terms\":[";
      bool first_term=true;
      for(DYTurbo::TermIterator term;!term.IsEnd();++term) {
        auto& t=*term;t.RunIntegration();
        if(!first_term)out<<',';first_term=false;
        out<<"{\"name\":"<<quoted(t.name)<<",\"method\":"<<quoted(t.isVegas?"vegas":"quadrature")
           <<",\"value\":"<<finite(t.last_int.at(0))<<",\"error\":"<<finite(std::sqrt(t.last_err2))
           <<",\"seconds\":"<<finite(t.last_time)<<'}';
      }
      if(first_term)throw std::runtime_error("No active QCD integration terms");
      out<<"],\"value\":"<<finite(DYTurbo::subtotal.last_int.at(0))
         <<",\"error\":"<<finite(std::sqrt(DYTurbo::subtotal.last_err2))<<'}';
      DYTurbo::subtotal.last_reset();
    }
    DYTurbo::Terminate();
    if(first_bin)throw std::runtime_error("No production bins");
    out<<"],\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";
    std::ofstream file(filename);file<<out.str();file.close();
    if(!file)throw std::runtime_error("Cannot write native results");
    return 0;
  }catch(const std::exception& e) {std::cerr<<"leptonff DYTurbo: "<<e.what()<<'\n';return 2;}
}
