
xmin = -160
xmax = 230

[Mesh]
  [1d]
    type = GeneratedMeshGenerator
    dim = 1
    xmin = -160
    xmax = 230
    nx = 3
    elem_type = edge3
  []
[]
[Variables]
  [u]
    order = second
    family = lagrange
  []
[]

[AuxVariables]
  [diffFluxU_x]
    order = FIRST
    family = MONOMIAL
  []
[]      

[Kernels]
  [diffusion-term]
    type = DiffusionTerm
    variable = u
    diffCoeff = 1800
  []
  [source-term]
    type = SourceTerm
    variable = u
    sourceS = -1.13
  []
[]

[AuxKernels]
  [diffusion-flux-x]
    execute_on = timestep_end
    type = DiffusionFluxComponent
    diffusivity = 1800
    variable = diffFluxU_x
    field = u
    component = x
    param1 = 18
    param2 = 1
  []
[]

[BCs]
  [entry]
    type = DirichletBC
    variable = u
    boundary = left
    value = 3.1
  []
  [exit]
    type = DirichletBC
    variable = u
    boundary = right
    value = 27
  []
[]

[Executioner]
  type = Steady
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[VectorPostprocessors]
  [u]
    type = LineValueSampler
    execute_on = 'timestep_end final'
    variable = 'u'
    start_point = '-160 0 0'
    end_point = '230 0 0'
    num_points = 20
    sort_by = id
  []
  [x-data]
    type = LineValueSampler
    execute_on = 'timestep_end final'
    variable = 'u diffFluxU_x'
    start_point = '${replace xmin} 0 0'
    end_point = '${replace xmax} 0 0'
    num_points = 20
    sort_by = id
  []
[]

[Postprocessors]
  [bulk-energy]
    type = BulkEnergy
    execute_on = 'timestep_end final'
    variable2 = 'u'
    variable1 = 'diffFluxU_x'
    param1 = 1800
    param2 = -1.13
  []
[]

[Outputs]
  [csv]
    type = CSV
    file_base = 'output'
    execute_on = 'final'   
  []
[]
