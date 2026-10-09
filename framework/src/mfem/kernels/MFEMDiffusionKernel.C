//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMDiffusionKernel.h"
#include "MFEMProblem.h"
#include "MFEMKernel.h"
#include "MooseStringUtils.h"
#include "MooseUtils.h"

#include <functional>

registerMooseObject("MooseApp", MFEMDiffusionKernel);

InputParameters
MFEMDiffusionKernel::validParams()
{
  InputParameters params = MFEMKernel::validParams();
  params.addClassDescription("Adds the domain integrator to an MFEM problem for the bilinear form "
                             "$(k\\vec\\nabla u, \\vec\\nabla v)_\\Omega$ "
                             "arising from the weak form of the Laplacian operator "
                             "$- \\vec\\nabla \\cdot \\left( k \\vec \\nabla u \\right)$.");
  params.addParam<MFEMScalarCoefficientName>(
      MFEMKernel::COEFFICIENT_PARAM, "1.", "Name of property for diffusion coefficient k.");
  params.addParam<MFEMMatrixCoefficientName>(MFEMKernel::MATRIX_COEFFICIENT_PARAM,
                                             "Name of matrix coefficient for property k. Mutually "
                                             "exclusive with parameter 'coefficient'.");
  return params;
}

MFEMDiffusionKernel::MFEMDiffusionKernel(const InputParameters & parameters)
  : MFEMKernel(parameters)
{
  /*
   MatrixFunctionCoefficient(int dim,
                             std::function<void(const Vector &, DenseMatrix &)> F,
                             Coefficient *q = nullptr)
   */
  // 1. Parse the matrix_coefficient terms, if present.
  // 2. If a term is a number, convert to a new lambda that returns the constant and add to
  //    function list
  // 3. If term is a variable, look for a function with that name.
  // 4a. If function does not exist, error.
  // 4b. If function exists, add it to function list (may need to be wrapped?).
  // 5. Create MFEM MatrixFunctionCoefficient compatible functions that populates a DenseMatrix
  //    using the function list.

  if (!_pars.isParamSetByUser(MFEMKernel::MATRIX_COEFFICIENT_PARAM))
  {
    return;
  }
  const auto & mat_coeff_str =
      _pars.get<MFEMMatrixCoefficientName>(MFEMKernel::MATRIX_COEFFICIENT_PARAM);
  const auto mat_coeffs = MooseUtils::split(mat_coeff_str, " ");
  using FunctionType =
      std::function<double(/* time */ double, /* x */ double, /* y */ double, /* z */ double)>;
  std::vector<FunctionType> functions;

  std::size_t cols = 0;
  std::size_t current_row_size = 0;
  std::size_t rows = 0;
  for (auto item : mat_coeffs)
  {
    item = MooseUtils::trim(item);
    if (item.empty())
    {
      continue;
    }

    bool end_of_row = false;
    if (item.back() == ';')
    {
      end_of_row = true;
      item.pop_back();
    }

    FunctionType fn;
    double value;
    if (MooseUtils::convert<double>(item, value, false))
    {
      fn = [value](double, double, double, double) { return value; };
    }
    else if (getMFEMProblem().hasFunction(item))
    {
      auto * func = &getMFEMProblem().getFunction(item);
      fn = [func](double t, double x, double y, double z) { return func->value(t, {x, y, z}); };
    }
    else
    {
      mooseError("matrix_coefficient ",
                 item,
                 " cannot be converted to a scalar and does not name a function.");
    }

    functions.push_back(std::move(fn));
    ++current_row_size;

    if (end_of_row)
    {
      ++rows;
      if (cols == 0)
      {
        cols = current_row_size;
      }
      else if (current_row_size != cols)
      {
        mooseError("matrix_coefficient contains rows with different sizes.");
      }
      current_row_size = 0;
    }
  }

  // Make the matrix coefficient's function.
  std::function<void(const mfem::Vector &, mfem::real_t, mfem::DenseMatrix &)> matrix_function =
      [functions = std::move(functions),
       cols](const mfem::Vector & v, mfem::real_t t, mfem::DenseMatrix & mat)
  {
    mat.SetSize(cols);
    for (auto i = 0U; i < cols; ++i)
    {
      for (auto j = 0U; j < cols; ++j)
      {
        mat(i, j) = functions[j + cols * i](t, v(0), v(1), v(2));
      }
    }
  };

  getMFEMProblem().getCoefficients().declareMatrix<mfem::MatrixFunctionCoefficient>(
      MFEMKernel::MATRIX_COEFFICIENT_PARAM, cols, matrix_function);
}

mfem::BilinearFormIntegrator *
MFEMDiffusionKernel::createBFIntegrator()
{
  auto coeffs = getMFEMProblem().getCoefficients().resolveCoefficientVariant(
      _pars, MFEMKernel::COEFFICIENT_PARAM, MFEMKernel::MATRIX_COEFFICIENT_PARAM);
  return std::visit([](auto & c) { return new mfem::DiffusionIntegrator(c); }, coeffs);
}

#endif
