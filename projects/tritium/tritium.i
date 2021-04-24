


[Mesh]
 [1d]                               # any lower case name; this is our 1D block
  type = GeneratedMeshGenerator  
  dim = 1 
  xmin = 0 
  xmax = 2 							#[mm]
  nx = 10
  elem_type = edge3 
 []
[]

[Variables]
 [u]                                # our unknown variable to be related to below
  order = second
  family = lagrange
 []
[]

[Kernels]
 [diffusion-term]                   # our diffusion term kernel implemented in the app
  type = DiffusionTerm              # our name for the diffusion kernel C++ class
  variable = u                      # u in this case is the concentration of Tritium 
  diffCoeff = 63200            		# Diffusivity of Tritium into SS316[mm2/s] [Trident Source]
 []
 [source-term]                      # our source term kernel implemented in the app
  type = SourceTerm                 # our name for the source kernel C++ class
  variable = u                      # our unknown variable from the [Variables] block
  sourceS = 0               		# Currently No Source (Diffusion through SS316)
 []
[]

[BCs]
 [entry]                            # inlet boundary
 type = NormalFluxBC           		# Diffusion of Tritium into SS316 is limited by the mass transfer through pipe
 variable = u
 boundary = left     
 transferCoeff = 24e+03 				# estimated (see Excel Sheet) [mm/s]
 reference = 3.6e-16      			# Inital Concentration of Tritium in FliNak for an AHTR [g/mm] design [source Tritium in ATHR] [3]
 []
 [exit]                             # outlet boundary
  type = DirichletBC				#DirichletBc was imposed by the assumption that at tritium is removed in the secondary fluid
  variable = u                      # our unknown variable from the [Variables] block
  boundary = right
  value = 0.0
 []
[]

[Executioner]
 type = Steady
 solve_type = 'PJFNK'
 petsc_options_iname = '-pc_type -pc_hypre_type'
 petsc_options_value = 'hypre boomeramg'
[]

[VectorPostprocessors]
 [u]                                # our unknown variable
  type = LineValueSampler
  execute_on = 'timestep_end final'
  variable = 'u'                    # our unknown variable from the [Variables] block
  start_point = '0 0 0'
  end_point = '2 0 0'
  num_points = 7
  sort_by = id
 []
[]

[Outputs]

[tecplot]
    type = Tecplot
   tecplot = true
  []

 [csv]                             # our choice for data output: tabular
  type = CSV
  file_base = 'output'
  execute_on = 'final'
  
  
 []
[]