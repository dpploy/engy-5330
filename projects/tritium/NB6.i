
#  Engy-5310 Problem: Poisson 2D FEM
# UMass Lowell Nuclear Chemical Engineering
# Prof. Valmor F. de Almeida
# 11Apr21 15:21:36

# Parameters
xmin = 0.00000e+00
xmax = 2.50000e+01
ymin = -6.25000e+00
ymax = 6.25000e+00
diff_coeff = 1.00000e-01
source_s = 0
u_left = 3.00000e+00
u_right = 0.00000e+00

[Problem]
  type = FEProblem
  coord_type = XYZ
[]

[Mesh]
  [2d]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = ${replace xmin}
    xmax = ${replace xmax}
    ymin = ${replace ymin}
    ymax = ${replace ymax}
    nx = 3
    ny = 1
	elem_type= QUAD9
  []
[]

[Variables]
  [u]
    order = second
    family = lagrange
  []
[]

[Kernels]
  [diffusion-term]
    type = DiffusionTerm
    variable = u     # produced quantity
    diffCoeff = ${replace diff_coeff}
  []
  [source-term]
    type = SourceTerm
    variable = u     # add to produced quantity
    sourceS = ${replace source_s}
  []
[]

[BCs]
  [east]
    type = DirichletBC
    variable = u
    boundary = left
    value = ${replace u_left}
  []
  [west]
    type = NeumannBC
    variable = u
    boundary = right
    value = ${replace u_right}
  []
  [south]
    type = NeumannBC
    variable = u
    boundary = bottom
  []
  [north]
    type = NeumannBC
    variable = u
    boundary = top
  []
[]

[Executioner]
  type = Steady
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[VectorPostprocessors]
  [x-line]
    type = LineValueSampler
    execute_on = 'timestep_end final'
    variable = 'u'    # output data
    start_point = '${replace xmin} ${fparse (ymax+ymin)/2} 0'
    end_point = '${replace xmax} ${fparse (ymax+ymin)/2} 0'
    num_points = 20
    sort_by = id
  []
  [y-line]
    type = LineValueSampler
    execute_on = 'timestep_end final'
    variable = 'u'    # output data
    start_point = '${fparse (xmax+xmin)/2} ${replace ymin} 0'
    end_point = '${fparse (xmax+xmin)/2} ${replace ymax} 0'
    num_points = 2
    sort_by = id
  []
[]

[Outputs]
  console = true
  [vtk]
    type = VTK
    execute_on = final
    file_base = out
  []
  [x]
    type = CSV
    execute_on = 'final'
    show = 'x-line'
    file_base = out-x
  []
  [y]
    type = CSV
    execute_on = 'final'
    show = 'y-line'
    file_base = out-y
  []
[]