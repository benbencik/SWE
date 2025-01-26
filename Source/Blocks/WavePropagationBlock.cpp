/**
 * @file
 * This file is part of SWE.
 *
 * @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Sebastian Rettenberger (rettenbs AT in.tum.de,
 * http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
 * @author Michael Bader (bader AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Michael_Bader)
 *
 * @section LICENSE
 *
 * SWE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SWE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SWE.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * @section DESCRIPTION
 *
 * Implementation of Blocks::Block that uses solvers in the wave propagation formulation.
 */

#include "WavePropagationBlock.hpp"

#include <iostream>
#include <omp.h>

Blocks::WavePropagationBlock::WavePropagationBlock(int nx, int ny, RealType dx, RealType dy):
  Block(nx, ny, dx, dy),
  hNetUpdatesLeft_(nx + 1, ny),
  hNetUpdatesRight_(nx + 1, ny),
  huNetUpdatesLeft_(nx + 1, ny),
  huNetUpdatesRight_(nx + 1, ny),
  hNetUpdatesBelow_(nx, ny + 1),
  hNetUpdatesAbove_(nx, ny + 1),
  hvNetUpdatesBelow_(nx, ny + 1),
  hvNetUpdatesAbove_(nx, ny + 1) {}

Blocks::WavePropagationBlock::WavePropagationBlock(
  int nx, int ny, RealType dx, RealType dy, Tools::Float2D<RealType>& h, Tools::Float2D<RealType>& hu, Tools::Float2D<RealType>& hv
):
  Block(nx, ny, dx, dy, h, hu, hv),
  hNetUpdatesLeft_(nx + 1, ny),
  hNetUpdatesRight_(nx + 1, ny),
  huNetUpdatesLeft_(nx + 1, ny),
  huNetUpdatesRight_(nx + 1, ny),
  hNetUpdatesBelow_(nx, ny + 1),
  hNetUpdatesAbove_(nx, ny + 1),
  hvNetUpdatesBelow_(nx, ny + 1),
  hvNetUpdatesAbove_(nx, ny + 1) {}

void Blocks::WavePropagationBlock::computeNumericalFluxes() {
  // Maximum (linearized) wave speed within one iteration
  RealType maxWaveSpeed = RealType(0.0);

#pragma omp parallel
  {
    // Initialize two local maximum wave speeds for vertical and horizontal edges
    RealType maxWaveSpeedLocal_u = RealType(0.0);
    RealType maxWaveSpeedLocal_v = RealType(0.0);

#pragma omp for
    for (int i = 1; i < nx_ + 1; i++) {
#pragma omp simd reduction(max : maxWaveSpeedLocal_u) reduction(max : maxWaveSpeedLocal_v)
      for (int j = 1; j < ny_ + 1; ++j) {
        // Compute the net-updates for the vertical edges
        RealType maxEdgeSpeed = RealType(0.0);
        wavePropagationSolver_.computeNetUpdates(
          h_[i - 1][j],
          h_[i][j],
          hu_[i - 1][j],
          hu_[i][j],
          b_[i - 1][j],
          b_[i][j],
          hNetUpdatesLeft_[i - 1][j - 1],
          hNetUpdatesRight_[i - 1][j - 1],
          huNetUpdatesLeft_[i - 1][j - 1],
          huNetUpdatesRight_[i - 1][j - 1],
          maxEdgeSpeed
        );

        // Update the thread-local maximum wave speed
        maxWaveSpeedLocal_u = std::max(maxWaveSpeed, maxEdgeSpeed);

        // Compute the net-updates for the horizontal edges
        maxEdgeSpeed = RealType(0.0);
        wavePropagationSolver_.computeNetUpdates(
          h_[i][j - 1],
          h_[i][j],
          hv_[i][j - 1],
          hv_[i][j],
          b_[i][j - 1],
          b_[i][j],
          hNetUpdatesBelow_[i - 1][j - 1],
          hNetUpdatesAbove_[i - 1][j - 1],
          hvNetUpdatesBelow_[i - 1][j - 1],
          hvNetUpdatesAbove_[i - 1][j - 1],
          maxEdgeSpeed
        );

        // Update the thread-local maximum wave speed
        maxWaveSpeedLocal_v = std::max(maxWaveSpeed, maxEdgeSpeed);
      }
    }

// Dealing with the vertical edges for cells on the nx_+1 boundary
#pragma omp for
    for (int j = 1; j < ny_ + 1; ++j) {
      RealType maxEdgeSpeed = RealType(0.0);
      wavePropagationSolver_.computeNetUpdates(
        h_[nx_][j],
        h_[nx_ + 1][j],
        hu_[nx_][j],
        hu_[nx_ + 1][j],
        b_[nx_][j],
        b_[nx_ + 1][j],
        hNetUpdatesLeft_[nx_][j - 1],
        hNetUpdatesRight_[nx_][j - 1],
        huNetUpdatesLeft_[nx_][j - 1],
        huNetUpdatesRight_[nx_][j - 1],
        maxEdgeSpeed
      );
      maxWaveSpeedLocal_u = std::max(maxWaveSpeed, maxEdgeSpeed);
    }

// Dealing with the horizontal edges for cells on the ny_+1 boundary
#pragma omp for
    for (int i = 1; i < nx_ + 1; ++i) {
      RealType maxEdgeSpeed = RealType(0.0);
      wavePropagationSolver_.computeNetUpdates(
        h_[i][ny_],
        h_[i][ny_ + 1],
        hv_[i][ny_],
        hv_[i][ny_ + 1],
        b_[i][ny_],
        b_[i][ny_ + 1],
        hNetUpdatesBelow_[i - 1][ny_],
        hNetUpdatesAbove_[i - 1][ny_],
        hvNetUpdatesBelow_[i - 1][ny_],
        hvNetUpdatesAbove_[i - 1][ny_],
        maxEdgeSpeed
      );
      maxWaveSpeedLocal_v = std::max(maxWaveSpeed, maxEdgeSpeed);
    }

#pragma omp critical
    { maxWaveSpeed = std::max(std::max(maxWaveSpeedLocal_u, maxWaveSpeedLocal_v), maxWaveSpeed); }
  }


  if (maxWaveSpeed > 0.00001) {
    // Compute the time step width
    maxTimeStep_ = std::min(dx_ / maxWaveSpeed, dy_ / maxWaveSpeed);

    // Reduce maximum time step size by "safety factor"
    maxTimeStep_ *= RealType(0.4); // CFL-number = 0.5
  } else {
    // Might happen in dry cells
    maxTimeStep_ = std::numeric_limits<RealType>::max();
  }
}

#include <immintrin.h> // For SIMD intrinsics

void Blocks::WavePropagationBlock::updateUnknowns(RealType dt) {
  __m256d dt_dx  = _mm256_set1_pd(dt / dx_);
  __m256d dt_dy  = _mm256_set1_pd(dt / dy_);
  __m256d zero   = _mm256_set1_pd(0.0);
  __m256d dryTol = _mm256_set1_pd(0.1);

  // Update cell averages with the net-updates
  for (int i = 1; i < nx_ + 1; i++) {
    for (int j = 1; j < ny_ + 1; j += 4) { // Process 4 elements at a time for double precision
      __m256d h  = _mm256_loadu_pd(&h_[i][j]);
      __m256d hu = _mm256_loadu_pd(&hu_[i][j]);
      __m256d hv = _mm256_loadu_pd(&hv_[i][j]);

      __m256d hNetRight = _mm256_loadu_pd(&hNetUpdatesRight_[i - 1][j - 1]);
      __m256d hNetLeft  = _mm256_loadu_pd(&hNetUpdatesLeft_[i][j - 1]);
      __m256d hNetAbove = _mm256_loadu_pd(&hNetUpdatesAbove_[i - 1][j - 1]);
      __m256d hNetBelow = _mm256_loadu_pd(&hNetUpdatesBelow_[i - 1][j]);

      __m256d huNetRight = _mm256_loadu_pd(&huNetUpdatesRight_[i - 1][j - 1]);
      __m256d huNetLeft  = _mm256_loadu_pd(&huNetUpdatesLeft_[i][j - 1]);

      __m256d hvNetAbove = _mm256_loadu_pd(&hvNetUpdatesAbove_[i - 1][j - 1]);
      __m256d hvNetBelow = _mm256_loadu_pd(&hvNetUpdatesBelow_[i - 1][j]);

      __m256d hUpdate  = _mm256_add_pd(_mm256_mul_pd(dt_dx, _mm256_add_pd(hNetRight, hNetLeft)), _mm256_mul_pd(dt_dy, _mm256_add_pd(hNetAbove, hNetBelow)));
      __m256d huUpdate = _mm256_mul_pd(dt_dx, _mm256_add_pd(huNetRight, huNetLeft));
      __m256d hvUpdate = _mm256_mul_pd(dt_dy, _mm256_add_pd(hvNetAbove, hvNetBelow));

      h  = _mm256_sub_pd(h, hUpdate);
      hu = _mm256_sub_pd(hu, huUpdate);
      hv = _mm256_sub_pd(hv, hvUpdate);

      // Handle negative heights and dry cells
      __m256d maskNegative = _mm256_cmp_pd(h, zero, _CMP_LT_OQ);
      __m256d maskDry      = _mm256_cmp_pd(h, dryTol, _CMP_LT_OQ);

      h  = _mm256_blendv_pd(h, zero, maskNegative);
      hu = _mm256_blendv_pd(hu, zero, maskNegative);
      hv = _mm256_blendv_pd(hv, zero, maskNegative);

      hu = _mm256_blendv_pd(hu, zero, maskDry);
      hv = _mm256_blendv_pd(hv, zero, maskDry);

      _mm256_storeu_pd(&h_[i][j], h);
      _mm256_storeu_pd(&hu_[i][j], hu);
      _mm256_storeu_pd(&hv_[i][j], hv);
    }
  }
}
