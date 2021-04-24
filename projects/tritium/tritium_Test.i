
# Engy-5310 Problem: Poisson 2D FEM
# UMass Lowell Nuclear Chemical Engineering
# Prof. Valmor F. de Almeida
# 11Apr21 15:21:36

# Parameters
xmin = 0.00000e+00
xmax = 1.50000e+02
ymin = -1.5000e+02
ymax = 1.5000e+02
diff_coeff =  .000266493
source_s = 0
velocity = '0 0 0'

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
	elem_type = QUAD9
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
  [convection-term]
    type = ConvectionTerm
    variable = u     # produced quantity 
    velocity = ${replace velocity}
  []
[]

[BCs]
 [metal]                            # Liquid-metal interface
 type = NormalFluxBC           		# Diffusion of Tritium into SS316 is limited by the mass transfer through pipe
 variable = u
 boundary = right    
 transferCoeff = 2.35e-04 			# estimated (see Excel Sheet) [cm/s]
 reference = 0      			    #inside pipe no tritium
 []
 [center]                           #center line of bulk fuel salt
  type = DirichletBC				#DirichletBc was imposed by the assumption that at tritium is removed in the secondary fluid
  variable = u                      # our unknown variable from the [Variables] block
  boundary = left
  value =  3.6e-07
 []
  
  [outflow]
    type = NeumannBC
    variable = u
    boundary = bottom
    value = 0
  []
  [inflow]
    type = DirichletBC
    variable = u
    boundary = top
    value = 3.6e-07

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