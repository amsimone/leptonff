#pragma once
#include "leptonff/flavours.hpp"

namespace leptonff { namespace detail {
enum FlavourBasis { Pff,Pfg,Pgf,Delta,Photonic,SameInterference,OppositeInterference,
  NonSingletPair,SingletPair,Pfg1,Pgf1Direct,Pgf1Vacuum,Pgg1,Dff0,Dfg0,QCDgg,flavour_basis_count };
using FlavourExpansion=std::array<Real,flavour_basis_count>;
Distribution flavour_basis_distribution(std::size_t);
FlavourExpansion flavour_kernel_expansion(const QEDFlavourScheme&,const std::vector<bool>&,
    std::size_t i,std::size_t j,int order); // order 0/1
FlavourExpansion flavour_boundary_expansion(const QEDFlavourScheme&,const std::vector<bool>&,
    std::size_t i,std::size_t j,Real mu);
FlavourExpansion flavour_matching_expansion(const QEDFlavourScheme&,std::size_t heavy,
    std::size_t i,std::size_t j,Real mu);
FlavourExpansion qcd_kernel_expansion(const QEDFlavourScheme&,const std::vector<bool>&,
    std::size_t i,std::size_t j);
}} // namespace leptonff::detail
