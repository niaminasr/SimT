/**
* @brief Header file for the Matrix class.
*
* @author Dimitri Walther
* @date 11/03/2023
*/

#ifndef SIMT_MATRIX_HPP
#define SIMT_MATRIX_HPP

#include "ParEITConfig.hpp"
#include "geometry/Form.hpp"
#include "geometry/Grid.hpp"
#include "tools/EitBiVector.hpp"

class Matrix {
private:
  Grid *m_grid;
  Form *m_form;
  EitBiVector<real_t> *m_sigma;

  uint_t m_NxBeg, m_NxEnd;

public:
  // Constructor
  Matrix(Grid *grid, Form *form, EitBiVector<real_t> *sigma);
  // Copy Constructor
  Matrix(const Matrix &that);
  // Destructor
  ~Matrix(){};

  Vertexint getIJKfromglobalidx(uint_t idxglob);

  /*
      The neighbors stored in class Point are determined by their global index,
      so we need to ransform them into local indices, when we work locally.
      The function Computelocalindex does so by comunicationg via a breadcast the number of nodes managed on each processor
      and then the local index is computed by substracting the effective of nodes on all previos processors from the global index.
   */

  uint_t computeLocalIndex(uint_t idxglob);

  bool getMyProc(uint_t idx);

  // Scalar product definition
  EitBiVector<real_t> operator*(const EitBiVector<real_t> &u);
};

// EitBiVector<real_t> operator* (const Matrix &A, const EitBiVector<real_t> &u);

#endif /* SRC_SOLVER_MATRIX */
