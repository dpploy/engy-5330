# Student: Anadi Mondal
# Mentor: Dr. Subash Sharma
# Instructor:Prof. Valmor F. de Almeida


# Parameters
xmin = 0
xmax = 6
diff_coeff = 5e-05
#source_s = 0.6883
#source_transfer_coeff = 5.00000e-03
#source_saturation = 1.00000e+00
u_left = 0.0459
u_right = 0
u2_left = 0
u2_right = 0
diff_coeff_2 = 9.2e-05
velocity = '14735.8 0.00000e+00 0.00000e+00'

[Problem]
  type = FEProblem
  coord_type = XYZ
[]

[Mesh]
  [1d]
    type = GeneratedMeshGenerator
    dim = 1
    xmin = ${replace xmin}
    xmax = ${replace xmax}
    nx = 20
	elem_type= edge3
  []
[]

[Variables]
  [u]
    order = second
    family = lagrange
    #initial_condition = ${fparse (u_right+u_left)/2}
  []
  [u2]
    order = second
    family = lagrange
    #initial_condition = ${fparse (u2_right+u2_left)/2}
  []
[]

[AuxVariables]
  [diffFluxU]
    order = CONSTANT
    family = MONOMIAL_VEC
  []
  [diffFluxU2]
    order = CONSTANT
    family = MONOMIAL_VEC
  []
  [diffFluxU_x]
    order = CONSTANT
    family = MONOMIAL
  []
  [diffFluxU2_x]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Kernels]
  [diffusion-term]
    type = DiffusionTerm
    variable = u     # produced quantity
    diffCoeff = ${replace diff_coeff}
	voidFraction =0.998
  []
  [source-term]
    type = SourceTerm
    variable = u     # add to produced quantity
    #sourceS = ${replace source_s}
	massArea = 15
	henryC =2.08e-04
    #transferCoeff = ${replace source_transfer_coeff}
    #saturation = ${replace source_saturation}
    coupledVariable = u2
  []
  [convection-term]
    type = ConvectiveTerm
    variable = u     # produced quantity
    velocity = ${replace velocity}
	voidFraction =0.998
  []
  [diffusion-term-2]
    type = DiffusionTermm
    variable = u2     # produced quantity
    diffCoeff = ${replace diff_coeff_2}
	voidFractionn =0.002
  []
  [source-term-2]
    type = SourceTerm
    variable = u2     # add to produced quantity
	#field = u
    #sourceSCoupled = ${replace source_s}
	massC= 15
	heC =2.08e-04
    #transferCoeffCoupled = ${replace source_transfer_coeff}
    #saturationCoupled = ${replace source_saturation}
	coupledVariable  = u
  []
  [convection-term-2]
    type = ConvectiveTermm
    variable = u2     # produced quantity
    velocity = ${replace velocity}
	voidFractionn =0.002
  []
[]

[AuxKernels]

  [diffusion-flux-x]
    execute_on = timestep_end
    type = DiffusionFluxComponent
    field = u
    diffCoeff = ${replace diff_coeff}
    component = x
    variable = diffFluxU_x     # produced quantity
  []

  [diffusion-flux-x-2]
    execute_on = timestep_end
    type = DiffusionFluxComponent
	field = u2
    diffCoeff = ${replace diff_coeff_2}   # produced quantity
    component = x
    variable = diffFluxU2_x   
  []
[]

[BCs]
  [entry-u]
    type = DirichletBC
    variable = u
    boundary = left
    value = ${replace u_left}

  []
  [entry-u2]
    type = DirichletBC
    variable = u2
    boundary = left
    value = ${replace u2_left}

  []
  [exit-u]
    type = NeumannBC
    variable = u
    boundary = right
    value = ${replace u_right}

  []
  [exit-u2]
    type = NeumannBC
    variable = u2
    boundary = right
    value = ${replace u2_right}
	
  []
[]

[Preconditioning]
  active = 'fdp-newt-full'
  [fdp-newt-full]
    type = FDP
    full = true
    solve_type = 'NEWTON'
    petsc_options_iname = '-pc_type -mat_fd_coloring_err -mat_fd_type'
    petsc_options_value = 'lu        9.9999999999999995474811182588626e-05              ds'
  []
[]

[Executioner]
  type = Steady
[]

[VectorPostprocessors]
  [x-data]
    type = LineValueSampler
    execute_on = 'timestep_end final'
    variable = 'u u2 diffFluxU_x diffFluxU2_x'    # output data
    start_point = '${replace xmin} 0 0'
    end_point = '${replace xmax} 0 0'
    num_points = 21
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
