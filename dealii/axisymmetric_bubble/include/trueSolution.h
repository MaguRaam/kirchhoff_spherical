#ifndef TRUESOLUTION_H_
#define TRUESOLUTION_H_

#include <deal.II/base/function.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/lac/vector.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include "Headers.h"

Point<2> grid_transform (Point<2>& );

class SineTransform : public ChartManifold<2,2>
{
public:
	virtual
  Point<2>
  pull_back(const Point<2> &space_point) const override;
  virtual
  Point<2>
  push_forward(const Point<2> &chart_point) const override;

  virtual std::unique_ptr<Manifold<2,2> > clone() const override;

};

class TrueSolution : public Function<2>
{
public:
  TrueSolution() : Function<2>(1)
  {}

  void soln(const Point<2> &, Vector<double>&, const double &) const;

  void soln_list(const std::vector<Point<2> >&,
                             std::vector< Vector<double> >&, const double&) const;

};

#endif
