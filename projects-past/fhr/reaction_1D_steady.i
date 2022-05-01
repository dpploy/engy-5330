# Steady-state test for the InterfaceReaction kernel.
#
# Specie M transport from domain 1 (0<=x<=1) to domain 2 (1<x<=2),
# u and v are concentrations in domain 1 and domain 2.
#
# Diffusion in both domains can be described by Ficks law and diffusion
# kernel is applied.
#
# Specie M has different diffusity in different domains, here set as D1=4, D2=2.
#
# Dirichlet boundary conditions are applied, i.e., u(0)=1, v(2)=0
#
# At the interface consider the following
#
# (a) Fluxes are matched from both domains (InterfaceDiffusion kernel)
#
# (b) First-order reaction is R = kf*u - kb*v
#
# Analytical solution is
# u = -0.2*u+1,    0<=u<=1
# v = -0.4*v+0.8,  1<v<=2
[Problem]
  type = FEProblem
  coord_type = RZ
  rz_coord_axis = Y
[]

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 10
    xmax = 2
    elem_type = edge3
  []
  [subdomain1]
    input = gen
    type = SubdomainBoundingBoxGenerator
    bottom_left = '1.0 0 0'
    block_id = 1
    top_right = '2.0 1.0 0'
  []
  [interface]
    type = SideSetsBetweenSubdomainsGenerator
    input = 'subdomain1'
    primary_block = '0'
    paired_block = '1'
    new_boundary = 'primary0_interface'
  []
[]

[Variables]
  [u]
    order = second
    family = LAGRANGE
    block = '0'
  []
  [v]
    order = second
    family = LAGRANGE
    block = '1'
  []
[]

[Kernels]
  [diff_u]
    type = MatDiffusion
    variable = u
    block = '0'
    diffusivity = D
  []
  #[sourceu]
  #  type = BodyForce
  #  variable = u
  #  function = 2
  #[]
  [sourceu]
    type = NuclearHeat
    block = 0
    variable = u
    sourceS = 2 
  []
  [diff_v]
    type = MatDiffusion
    variable = v
    block = '1'
    diffusivity = D
  []
  [sourcev]
    type = NuclearHeat
    block = 1
    variable = v
    sourceS = -2
  []
[]

[InterfaceKernels]
  [./interface]
    type = InterfaceDiffusion
    variable = u
    neighbor_var = 'v'
    boundary = 'primary0_interface'
    D = D
    D_neighbor = D
  [../]
  [./interface_reaction]
    type = InterfaceReaction
    variable = u
    neighbor_var = 'v'
    boundary = 'primary0_interface'
    kf = 1 # Forward reaction rate coefficient
    kb = 1 # Backward reaction rate coefficient
  [../]
[]

[BCs]
  #[left]
  #  type = DirichletBC
  #  variable = u
  #  boundary = 'left'
  #  value = 1
  #[]
  [left]
    type = NeumannBC
    variable = u
    boundary = 'left'
    value = 0
  []
  #[right]
  #  type = DirichletBC
  #  variable = v
  #  boundary = 'right'
  #  value = 0
  #[]
  [ro]
    type = NormalHeatFluxBC
    variable = v
    boundary = 'right'
    refTempFunc = refTempFunc
    transferCoeff = 1
  []
[]
[Functions]
  [refTempFunc]
    type = ParsedFunction
    value = temp_ref
    vars = temp_ref
    vals = 5
  []
[]

[Materials]
  [./block0]
    type = GenericConstantMaterial
    block = '0'
    prop_names = 'D'
    prop_values = '2'
  [../]
  [./block1]
    type = GenericConstantMaterial
    block = '1'
    prop_names = 'D'
    prop_values = '1'
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  nl_rel_tol = 1e-10
[]

[Outputs]
  print_linear_residuals = true
  execute_on = 'FINAL'
  exodus = true
  csv = true
[]

[Debug]
  show_var_residual_norms = true
[]

[Postprocessors]
  [elemental_error_u]
    type = ElementL2Error
    function = -0.2*x+1
    variable = 'u'
    block = '0'
  []
  [elemental_error_v]
    type = ElementL2Error
    function = -0.4*x+0.8
    variable = 'v'
    block = '1'
  []
[]
[VectorPostprocessors]
  [omega_1]
    type = LineValueSampler
    execute_on = 'timestep_end final'
    variable = 'u'
    start_point = '0 0 0'
    end_point = '${fparse 1.0-1e-3} 0 0'
    num_points = 11
    sort_by = id
  []
  [omega_2]
    type = LineValueSampler
    execute_on = 'timestep_end final'
    variable = 'v'
    start_point = '${fparse 1.0+1e-3} 0 0'
    end_point = '2.0 0 0'
    num_points = 11
    sort_by = id
  []
[]
