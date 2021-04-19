# Engy-5310 Problem 1: Poisson 1D FEM
# UMass Lowell Nuclear Chemical Engineering
# Prof. Valmor F. de Almeida
# 29Mar21 20:24:18

# Parameters
xmin = 0.00000e+00
xmax = 2.50000e+01
diff_coeff = 1.00000e-01
source_s = 1.00000e-03
rho_v = 1
rho_l = 1000
[Mesh]
  [1d]
    type = GeneratedMeshGenerator
    dim = 1
    xmin = ${replace xmin}
    xmax = ${replace xmax}
    nx = 3
  []
[]

[Variables]
  [velocityMixture]
    order = first
    family = lagrange
  []
  [fractionVapor]
    order = first
    family = lagrange
  []
[]

[Kernels]
  [mixture-mass-balance-divergence]
    type = MixtureMassBalDivergence
    variable = velocityMixture     # produced quantity
    diffCoeff = ${replace diff_coeff}
  []
  [vapor-drift-flux]
    type = VaporDriftDiffusion
    variable = fractionVapor     # produced quantity
    diffCoeff = ${replace diff_coeff}
  []
  [vapor-mass-transfer-source]
    type = VaporMassTransferSource
    variable = fractionVapor     # add to produced quantity
    sourceS = ${replace source_s}
  []
[]

[BCs]
  [entry-vapor-fraction]
    type = DirichletBC
    variable = fractionVapor
    boundary = left
    value = 0
  []
  [entry-mixture-velocity]
    type = DirichletBC
    variable = velocityMixture
    boundary = left
    value = 1
  []
  [exit]
    type = NeumannBC
    variable = fractionVapor
    boundary = right
    value = 0.00000e+00
  []
[]

[Executioner]
  type = Steady
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]


[VectorPostprocessors]
  [x-data]
    type = LineValueSampler
    execute_on = 'timestep_end final'
    variable = 'velocityMixture fractionVapor'    # output data
    start_point = '${replace xmin} 0 0'
    end_point = '${replace xmax} 0 0'
    num_points = 20
    sort_by = id
  []
[]

[Outputs]
  console = true
  [file-x-data]
    type = CSV
    file_base = 'output'
    execute_on = 'final'
    show = 'x-data'
  []
[]
