#pragma once
#include "leptonff/spin_decay.hpp"
void leptonff_photos_comparison(int,const double*,const double*,const double*,const double*,double*);
std::vector<std::array<leptonff::Real,4>> leptonff_photos_batch(int,const double*,const double*,
    const double*,const double*,const std::vector<leptonff::DecayMeasurement>&);
