#pragma once

#include <TCanvas.h>
#include <TTree.h>

struct RankingParams {
    double tan_low;
    double tan_up;
    int vph_standard_plus;
    double xy_lin_max;
    double lat_lin_max;
};

void position(TCanvas *c1, TTree *tree, int pl, const double *AreaParam, const std::vector<double> &TDRange) noexcept;

void position_projection(
    TCanvas *c1, TTree *tree, const size_t entries, int pl, const std::vector<double> &TDRange, const double *AreaParam
) noexcept;

void angle(TCanvas *c1, TTree *tree, int pl, const double angle_max, const double angle_resolution) noexcept;

void angle_projection(TCanvas *c1, TTree *tree, int pl, const double angle_max, const double angle_resolution) noexcept;

void d_angle(TCanvas *c1, TTree *tree) noexcept;

void d_angle_Ncut(
    TCanvas *c1, TTree *tree, const TString da_cutX, const TString da_cutY, const int da_cutPH
) noexcept;

void d_angle_rl(
    TCanvas *c1, TTree *tree, const double angle_max, const double dlat_range, const double drad_range, const int face
) noexcept;

void phvph_2D(
    TCanvas *c1, TTree *tree, const int vph_range, const double angle_max, const double angle_resolution
) noexcept;

void phvph_1D(TCanvas *c1, TTree *tree, const int vph_range, const uint8_t i, const double interval) noexcept;

void ranking(TCanvas *c1, TTree *tree, const std::vector<RankingParams> &params, uint32_t vph_standard) noexcept;

void dxdydz(TCanvas *c1, TTree *tree, const double *AreaParam) noexcept;

void dxdy(TCanvas *c1, TTree *tree, const double *AreaParam) noexcept;

void sensor_not(TCanvas *c1, const int fieldsOfView[2], const uint32_t NoTcount[2][72]) noexcept;

void sensor_drdz(
    TCanvas *c1, const uint32_t NoTcount_same[72], const double dr_sum[72], const double dz_sum[72]
) noexcept;
