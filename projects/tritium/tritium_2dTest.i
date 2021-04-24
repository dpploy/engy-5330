
# Engy-5310 Problem: Poisson 2D FEM
# UMass Lowell Nuclear Chemical Engineering
# Prof. Valmor F. de Almeida
# 11Apr21 15:22:32

# Parameters [cm]
xmin = 0
xmax = 150
ymin = 0
ymax = 300
diff_coeff =  1
source_s = 0
#u_right =0
u_left = 0
u_bottom = 0 
u_top = 10
velocity = '0 0 0'
PE = 135

[Problem]
  type = FEProblem
  #coord_type = RZ
  #rz_coord_axis = Y
    
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
    nx = 5
    ny = 10
	
	
  []
[]

[Variables]
  [u]
    order = first
    family = lagrange
  []
[]
[AuxVariables]
   [diffFluxU_x]
    order = FIRST
    family = MONOMIAL
	[]
	[diffFluxU_y]
    order = FIRST
    family = MONOMIAL
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
   [convection-term]
    type = ConvectionTerm
    variable = u     # produced quantity 
    velocity = ${replace velocity}
  []
  
[]
[AuxKernels]
  [diffusion-flux-x]
    execute_on = timestep_end
    type = DiffusionFluxComponent       # user-built auxiliary kernel
    field = u                           # user-defined parameter
    diffCoeff = ${replace diff_coeff}  # user-defined parameter
    component = x                       # user-defined parameter
    variable = diffFluxU_x              # produced quantity
  []
   [diffusion-flux-y]
    execute_on = timestep_end
    type = DiffusionFluxComponent       # user-built auxiliary kernel
    field = u                           # user-defined parameter
    diffCoeff = ${replace diff_coeff}   # user-defined parameter
    component = y                      # user-defined parameter
    variable = diffFluxU_y              # produced quantity
  []
[]

[BCs]
 [metal]                            # Liquid-metal interface
 type = NeumannBC          		# Diffusion of Tritium into SS316 is limited by the mass transfer through pipe
 variable = u
 boundary = right    
 value= 0      			        #inside pipe no tritium
 []
 [center]                           #center line of bulk fuel salt								
  type = NeumannBC				#DirichletBc was imposed by the assumption that at tritium is removed in the secondary fluid
  variable = u                      # our unknown variable from the [Variables] block
  boundary = left
  value = ${replace u_left}
 []
  
  [outflow]
    type = DirichletBC
    variable = u
    boundary = bottom
    value = ${replace u_bottom}
  []
  [inflow]
    type = DirichletBC
    variable = u
    boundary = top
    value = ${replace u_top}

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
    num_points = 20
    sort_by = id
  []
[]

[Outputs]
  console = true
  [tecplot]
    type = Tecplot
   tecplot = true
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