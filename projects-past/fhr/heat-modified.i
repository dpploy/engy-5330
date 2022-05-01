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
    nx = 40
    xmax = 2
    #elem_type = edge3
  []
  [subdomain1]
    input = gen
    type = SubdomainBoundingBoxGenerator
    block_id = 1
    bottom_left = '1 0 0'
    top_right = '2 1 0'
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
    order = first
    family = LAGRANGE
    block = '0'
    initial_condition = 1
  []
  [v]
    order = first
    family = LAGRANGE
    block = '1'
    initial_condition = 1
  []
[]

[Kernels]
  [diff_u]
    type = MatDiffusion
    variable = u
    block = '0'
    diffusivity = D
  []
  [sourceu]
    type = BodyForce
    block = 0
    variable = u
    function = 100
  []
  [diff_v]
    type = MatDiffusion
    variable = v
    block = '1'
    diffusivity = D
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
  [left]
    type = DirichletBC
    variable = u
    boundary = 'left'
    value = 10
  []
  [right]
    type = DirichletBC
    variable = v
    boundary = 'right'
    value = 5
  []
[]

[Materials]
  [./block0]
    type = GenericConstantMaterial
    block = '0'
    prop_names = 'D'
    prop_values = '1'
  [../]
  [./block1]
    type = GenericConstantMaterial
    block = '1'
    prop_names = 'D'
    prop_values = '5'
  [../]
[]
[Preconditioning]
  active = 'fdp-newt-full'
  [fdp-newt-full]
    type = FDP
    full = true
    solve_type = 'NEWTON'
    petsc_options_iname = '-pc_type -mat_fd_coloring_err -mat_fd_type'
    petsc_options_value = 'lu  1e-8          ds'
  []
[]
[Executioner]
  type = Steady
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

[VectorPostprocessors]
  [omega_1]
    type = LineValueSampler
    execute_on = 'timestep_end final'
    variable = 'u'
    start_point = '0 0 0'
    end_point = '${fparse 1-2e-3} 0 0'
    num_points = 21
    sort_by = id
  []
  [omega_2]
    type = LineValueSampler
    execute_on = 'timestep_end final'
    variable = 'v'
    start_point = '${fparse 1+1e-4} 0 0'
    end_point = '2 0 0'
    num_points = 21
    sort_by = id
  []
[]
