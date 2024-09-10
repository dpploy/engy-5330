#!/usr/bin/env python
# -*-  coding: utf-8 -*-
'''
This file is part of the Engy-5310 Computational Continuum Transport Phenomena
course at https://github.com/dpploy/engy-5310
'''

import os
import datetime
import subprocess

import math
import numpy as np

#*********************************************************************************
def engy5310_p1_exact_solution(x_a, x_b, u_a, u_b, diff_coeff, source_s):

    L = x_b - x_a

    a_hat = u_a/source_s*diff_coeff/L**2
    b_hat = u_b/source_s*diff_coeff/L**2

    u_hat = lambda x_hat, a_hat, b_hat: -x_hat**2/2 + (b_hat-a_hat)*x_hat + (a_hat+b_hat)/2 + 1/8

    flux_hat = lambda x_hat, a_hat, b_hat: -(-x_hat + b_hat - a_hat)

    # derive the correct energy
    energy_hat = (-b_hat + -a_hat)/2 + (-b_hat - -a_hat)**2/2 - 1/24

    return (u_hat, flux_hat, energy_hat)
#*********************************************************************************
def write_engy5310_p1_1d_input_file(x_left, x_right,
                                    u_left=None, u_right=None,
                                    qn_bias_left=None, transfer_coeff_left=None,
                                    u_reference_left=None,
                                    qn_bias_right=None, transfer_coeff_right=None,
                                    u_reference_right=None,
                                    diff_coeff=1.0,
                                    source_s=None,
                                    source_transfer_coeff=None,
                                    source_saturation=None,
                                    u2_left=None, u2_right=None,
                                    diff_coeff_2=None,
                                    source_s_2=None,
                                    source_transfer_coeff_2=None,
                                    source_saturation_2=None,
                                    velocity=None,
                                    n_felem=1, x_bias=None,
                                    order='second',
                                    n_plot_pts=10,
                                    compute_diffusion_flux=False,
                                    use_moose_diffusion_flux=False,
                                    compute_energy=False,
                                    use_moose_neumann_bc=False,
                                    solver='pjfnk-hypre-boom',
                                    fdp_dh=1e-6,
                                    file_name=None):

    if u_right is not None:
        assert qn_bias_right is None
        assert transfer_coeff_right is None
        assert u_reference_right is None
    else:
        assert qn_bias_right is not None or transfer_coeff_right is not None or \
                u_reference_right is not None

    if qn_bias_left is not None:
        assert u_right is None
    if transfer_coeff_left is not None:
        assert transfer_coeff_left >= 0.0
        assert u_right is None
    if u_reference_left is not None:
        assert u_right is None

    if velocity is not None:
        assert isinstance(velocity, tuple)
        assert len(velocity) == 3

    if x_bias is not None:
        assert 0.5 <= x_bias <= 2

    if file_name is None:
        fout = open('engy5310p1/input.hit', 'wt')
    else:
        fout = open(file_name, 'wt')

    fout.write('# Engy-5310 Problem 1: Poisson 1D FEM\n')
    fout.write('# UMass Lowell Nuclear Chemical Engineering\n')
    fout.write('# Prof. Valmor F. de Almeida\n')
    time_date_stamp = datetime.datetime.today().strftime('%d%b%y %H:%M:%S')
    fout.write('# %s\n'%time_date_stamp)

    fout.write('\n')

    fout.write('# Parameters\n')
    fout.write('xmin = %10.5e\n'%x_left)
    fout.write('xmax = %10.5e\n'%x_right)

    if diff_coeff is not None:
        fout.write('diff_coeff = %10.5e\n'%diff_coeff)
    if source_s is not None:
        fout.write('source_s = %10.5e\n'%source_s)
    if source_transfer_coeff is not None:
        fout.write('source_transfer_coeff = %10.5e\n'%source_transfer_coeff)
    if source_saturation is not None:
        fout.write('source_saturation = %10.5e\n'%source_saturation)

    if u_left is not None:
        fout.write('u_left = %10.5e\n'%u_left)
    if u_right is not None:
        fout.write('u_right = %10.5e\n'%u_right)

    if u2_left is not None:
        fout.write('u2_left = %10.5e\n'%u2_left)
    if u2_right is not None:
        fout.write('u2_right = %10.5e\n'%u2_right)

    if u2_left is not None or u2_right is not None:
        if diff_coeff_2 is not None:
            fout.write('diff_coeff_2 = %10.5e\n'%diff_coeff_2)
        if source_s_2 is not None:
            fout.write('source_s_2 = %10.5e\n'%source_s_2)
        if source_transfer_coeff_2 is not None:
            fout.write('source_transfer_coeff_2 = %10.5e\n'%source_transfer_coeff_2)
        if source_saturation_2 is not None:
            fout.write('source_saturation_2 = %10.5e\n'%source_saturation_2)

    if velocity is not None:
        fout.write("velocity = '%10.5e %10.5e %10.5e'\n"%(velocity[0],velocity[1],velocity[2]))

    if qn_bias_left is not None and use_moose_neumann_bc is False:
        fout.write('qn_bias_left = %10.5e\n'%qn_bias_left)
    if transfer_coeff_left is not None:
        fout.write('transfer_coeff_left = %10.5e\n'%transfer_coeff_left)
    if u_reference_left is not None:
        fout.write('u_reference_left = %10.5e\n'%u_reference_left)

    if qn_bias_right is not None and use_moose_neumann_bc is False:
        fout.write('qn_bias_right = %10.5e\n'%qn_bias_right)
    if transfer_coeff_right is not None:
        fout.write('transfer_coeff_right = %10.5e\n'%transfer_coeff_right)
    if u_reference_right is not None:
        fout.write('u_reference_right = %10.5e\n'%u_reference_right)

    fout.write('\n')

    fout.write('[Problem]\n')
    fout.write('  type = FEProblem\n')
    fout.write('  coord_type = XYZ\n')
    fout.write('[]\n')

    fout.write('\n')

    fout.write('[Mesh]\n')

    fout.write('  [omega-1d]\n')
    fout.write('    type = GeneratedMeshGenerator\n')
    fout.write('    dim = 1\n')
    fout.write('    xmin = ${replace xmin}\n')
    fout.write('    xmax = ${replace xmax}\n')
    fout.write('    nx = %i\n'%n_felem)
    if order == 'second':
        fout.write('    elem_type = edge3\n')
    if x_bias is not None:
        fout.write('    bias_x = %3.3e\n'%x_bias)
    fout.write('  []\n')

    fout.write('[]\n')

    fout.write('\n')

    fout.write('[Variables]\n')

    fout.write('  [u]\n')
    fout.write('    order = %s\n'%order)
    fout.write('    family = lagrange\n')
    if u_left is not None and u_right is None:
        fout.write('    initial_condition = ${replace u_left}\n')
    elif u_left is None and u_right is not None:
        fout.write('    initial_condition = ${replace u_right}\n')
    elif u_left is not None and u_right is not None:
        fout.write('    initial_condition = ${fparse (u_right+u_left)/2}\n')
    fout.write('  []\n')

    if u2_left is not None or u2_right is not None:
        fout.write('  [u2]\n')
        fout.write('    order = %s\n'%order)
        fout.write('    family = lagrange\n')
        if u2_left is not None and u2_right is None:
            fout.write('    initial_condition = ${replace u2_left}\n')
        elif u2_left is None and u2_right is not None:
            fout.write('    initial_condition = ${replace u2_right}\n')
        elif u2_left is not None and u2_right is not None:
            fout.write('    initial_condition = ${fparse (u2_right+u2_left)/2}\n')
        fout.write('  []\n')

    fout.write('[]\n')

    fout.write('\n')

    if compute_diffusion_flux:
        fout.write('[AuxVariables]\n')

        fout.write('  [diffFluxU]\n')
        if order in ['second', 'SECOND']:
            fout.write('    order = FIRST\n')
        elif order in ['first', 'FIRST']:
            fout.write('    order = CONSTANT\n')
        else:
            assert False
        fout.write('    family = MONOMIAL_VEC\n')
        fout.write('  []\n')

        if u2_left is not None or u2_right is not None:
            fout.write('  [diffFluxU2]\n')
            if order in ['second', 'SECOND']:
                fout.write('    order = FIRST\n')
            elif order in ['first', 'FIRST']:
                fout.write('    order = CONSTANT\n')
            else:
                assert False
            fout.write('    family = MONOMIAL_VEC\n')
            fout.write('  []\n')

        fout.write('  [diffFluxU_x]\n')
        if order in ['second', 'SECOND']:
            fout.write('    order = FIRST\n')
        elif order in ['first', 'FIRST']:
            fout.write('    order = CONSTANT\n')
        else:
            assert False
        fout.write('    family = MONOMIAL\n')
        fout.write('  []\n')

        if u2_left is not None or u2_right is not None:
            fout.write('  [diffFluxU2_x]\n')
            if order in ['second', 'SECOND']:
                fout.write('    order = FIRST\n')
            elif order in ['first', 'FIRST']:
                fout.write('    order = CONSTANT\n')
            else:
                assert False
            fout.write('    family = MONOMIAL\n')
            fout.write('  []\n')

        fout.write('[]\n')

        fout.write('\n')

    fout.write('[Kernels]\n')

    fout.write('  [diffusion-term]\n')
    fout.write('    type = DiffusionTerm\n')
    fout.write('    variable = u     # produced quantity\n')
    fout.write('    diffCoeff = ${replace diff_coeff}\n')
    fout.write('  []\n')

    if source_s is not None:
        fout.write('  [source-term]\n')
        fout.write('    type = SourceTerm\n')
        fout.write('    variable = u     # add to produced quantity\n')
        if source_s is not None:
            fout.write('    sourceS = ${replace source_s}\n')
        if source_transfer_coeff is not None:
            fout.write('    transferCoeff = ${replace source_transfer_coeff}\n')
        if source_saturation is not None:
            fout.write('    saturation = ${replace source_saturation}\n')
        if u2_left is None and u2_right is None:
            fout.write('    coupledVariable = u\n')
        else:
            fout.write('    coupledVariable = u2\n')
        if source_s_2 is not None:
            fout.write('    sourceSCoupled = ${replace source_s_2}\n')
        if source_transfer_coeff_2 is not None:
            fout.write('    transferCoeffCoupled = ${replace source_transfer_coeff_2}\n')
        if source_saturation_2 is not None:
            fout.write('    saturationCoupled = ${replace source_saturation_2}\n')
        fout.write('  []\n')

    if velocity is not None:
        fout.write('  [convection-term]\n')
        fout.write('    type = ConvectionTerm\n')
        fout.write('    variable = u     # produced quantity\n')
        fout.write('    velocity = ${replace velocity}\n')
        fout.write('  []\n')

    if u2_left is not None or u2_right is not None:

        fout.write('  [diffusion-term-2]\n')
        fout.write('    type = DiffusionTerm\n')
        fout.write('    variable = u2     # produced quantity\n')
        fout.write('    diffCoeff = ${replace diff_coeff_2}\n')
        fout.write('  []\n')

        fout.write('  [source-term-2]\n')
        fout.write('    type = SourceTerm\n')
        fout.write('    variable = u2     # add to produced quantity\n')
        if source_s_2 is not None:
            fout.write('    sourceS = ${replace source_s_2}\n')
        if source_transfer_coeff_2 is not None:
            fout.write('    transferCoeff = ${replace source_transfer_coeff_2}\n')
        if source_saturation_2 is not None:
            fout.write('    saturation = ${replace source_saturation_2}\n')
        if u_left is None and u_right is None:
            fout.write('    coupledVariable = u2\n')
        else:
            fout.write('    coupledVariable = u\n')
        if source_s is not None:
            fout.write('    sourceSCoupled = ${replace source_s}\n')
        if source_transfer_coeff is not None:
            fout.write('    transferCoeffCoupled = ${replace source_transfer_coeff}\n')
        if source_saturation is not None:
            fout.write('    saturationCoupled = ${replace source_saturation}\n')
        fout.write('  []\n')

        if velocity is not None:
            fout.write('  [convection-term-2]\n')
            fout.write('    type = ConvectionTerm\n')
            fout.write('    variable = u2     # produced quantity\n')
            fout.write('    velocity = ${replace velocity}\n')
            fout.write('  []\n')

    fout.write('[]\n')

    fout.write('\n')

    if compute_diffusion_flux:
        fout.write('[AuxKernels]\n')

        fout.write('  [diffusion-flux]\n')
        fout.write('    execute_on = timestep_end\n')
        if use_moose_diffusion_flux:
            fout.write('    type = DiffusionFluxAux     # provided by MOOSE\n')
            fout.write('    diffusion_variable = u\n')
            fout.write('    diffusivity = ${replace diff_coeff}\n')
            fout.write('    component = x\n')
            fout.write('    variable = diffFluxU_x     # produced quantity\n')
            fout.write('  []\n')
        else:
            fout.write('    type = DiffusionFlux\n')
            fout.write('    field = u\n')
            fout.write('    diffCoeff = ${replace diff_coeff}\n')
            fout.write('    variable = diffFluxU     # produced quantity\n')
            fout.write('  []\n')
            fout.write('  [diffusion-flux-x]\n')
            fout.write('    execute_on = timestep_end\n')
            fout.write('    type = VectorVariableComponentAux\n')
            fout.write('    variable = diffFluxU_x    # produced quantity\n')
            fout.write('    component = x\n')
            fout.write('    vector_variable = diffFluxU   \n')
            fout.write('  []\n')

        if u2_left is not None or u2_right is not None:
            fout.write('  [diffusion-flux-2]\n')
            fout.write('    execute_on = timestep_end\n')
            if use_moose_diffusion_flux:
                fout.write('    type = DiffusionFluxAux     # provided by MOOSE\n')
                fout.write('    diffusion_variable = u2\n')
                fout.write('    diffusivity = ${replace diff_coeff_2}\n')
                fout.write('    component = x\n')
                fout.write('    variable = diffFluxU2_x     # produced quantity\n')
                fout.write('  []\n')
            else:
                fout.write('    type = DiffusionFlux\n')
                fout.write('    field = u2\n')
                fout.write('    diffCoeff = ${replace diff_coeff_2}\n')
                fout.write('    variable = diffFluxU2     # produced quantity\n')
                fout.write('  []\n')
                fout.write('  [diffusion-flux-x-2]\n')
                fout.write('    execute_on = timestep_end\n')
                fout.write('    type = VectorVariableComponentAux\n')
                fout.write('    variable = diffFluxU2_x    # produced quantity\n')
                fout.write('    component = x\n')
                fout.write('    vector_variable = diffFluxU2   \n')
                fout.write('  []\n')

        fout.write('[]\n')
        fout.write('\n')

    fout.write('[BCs]\n')

    if u_left is not None:
        fout.write('  [entry-u]\n')
        fout.write('    type = DirichletBC\n')
        fout.write('    variable = u\n')
        fout.write('    boundary = left\n')
        fout.write('    value = ${replace u_left}\n')
        fout.write('  []\n')

    if u2_left is not None:
        fout.write('  [entry-u2]\n')
        fout.write('    type = DirichletBC\n')
        fout.write('    variable = u2\n')
        fout.write('    boundary = left\n')
        fout.write('    value = ${replace u2_left}\n')
        fout.write('  []\n')

    if u_right is not None:
        fout.write('  [exit-u]\n')
        fout.write('    type = DirichletBC\n')
        fout.write('    variable = u\n')
        fout.write('    boundary = right\n')
        fout.write('    value = ${replace u_right}\n')
        fout.write('  []\n')

    if u_right is None and use_moose_neumann_bc is False:
        fout.write('  [exit-u]\n')
        fout.write('    type = NormalFluxBC\n')
        fout.write('    variable = u\n')
        fout.write('    boundary = right\n')
        if qn_bias_right is not None:
           fout.write('    bias = ${replace qn_bias_right}\n')
        if transfer_coeff_right is not None:
            fout.write('    transferCoeff = ${replace transfer_coeff_right}\n')
        if u_reference_right is not None:
            fout.write('    reference = ${replace u_reference_right}\n')
        fout.write('  []\n')
    if u_right is None and use_moose_neumann_bc is True:
        fout.write('  [exit-u]\n')
        fout.write('    type = NeumannBC\n')
        fout.write('    variable = u\n')
        fout.write('    boundary = right\n')
        if qn_bias_right is not None:
            fout.write('    value = %10.5e\n'%(qn_bias_right))
        else:
            fout.write('      value = 0.0\n')
        fout.write('  []\n')

    if u2_right is not None:
        fout.write('  [exit-u2]\n')
        fout.write('    type = DirichletBC\n')
        fout.write('    variable = u2\n')
        fout.write('    boundary = right\n')
        fout.write('    value = ${replace u2_right}\n')
        fout.write('  []\n')

    fout.write('[]\n')

    fout.write('\n')

    if solver=='fdp-newt-full':
        fout.write('[Preconditioning]\n')
        fout.write("  active = 'fdp-newt-full'\n")
        fout.write('  [fdp-newt-full]\n')
        fout.write('    type = FDP\n')
        fout.write('    full = true\n')
        fout.write("    solve_type = 'NEWTON'\n")
        fout.write("    petsc_options_iname = '-pc_type -mat_fd_coloring_err -mat_fd_type'\n")
        #fout.write("    petsc_options_value = 'lu       1e-6                 ds'\n")
        fout.write("    petsc_options_value = 'lu  %10.3e               ds'\n"%fdp_dh)
        fout.write('  []\n')
        fout.write('[]\n')

        fout.write('\n')

        fout.write('[Executioner]\n')
        fout.write('  type = Steady\n')
        fout.write('[]\n')
    elif solver=='pjfnk-hypre-boom':
        fout.write('[Executioner]\n')
        fout.write('  type = Steady\n')
        fout.write("  solve_type = 'PJFNK'\n")
        fout.write("  petsc_options_iname = '-pc_type -pc_hypre_type'\n")
        fout.write("  petsc_options_value = 'hypre boomeramg'\n")
        fout.write('[]\n')
    else:
        assert False,'Invalid solver option = %r'%solver


    if compute_energy:

        fout.write('\n')

        fout.write('[Postprocessors]\n')
        fout.write('  [bulk-energy]\n')
        fout.write('    type = BulkEnergy\n')
        fout.write("    execute_on = 'timestep_end final'\n")
        fout.write("    variable = 'u'     # bulk energy unknown variable\n")
        fout.write('    diffCoeff = ${replace diff_coeff}\n')
        fout.write('    sourceS = ${replace source_s}\n')
        fout.write('  []\n')

        boundary_energy_right = False
        if qn_bias_right is not None or transfer_coeff_right is not None or u_reference_right is not None:
            fout.write('  [boundary-energy-right]\n')
            fout.write('    type = BoundaryEnergy\n')
            fout.write("    execute_on = 'timestep_end final'\n")
            fout.write("    variable = 'u'     # bulk energy unknown variable\n")
            fout.write('    diffCoeff = ${replace diff_coeff}\n')
            if qn_bias_right is not None:
                fout.write("    diffNormalFluxBias = ${replace qn_bias_right}\n")
            if transfer_coeff_right is not None:
                fout.write('    transferCoeff = ${replace transfer_coeff_right}\n')
            if u_reference_right is not None:
                fout.write('    reference = ${replace u_reference_right}\n')
            fout.write('    boundary = right\n')
            fout.write('  []\n')

            boundary_energy_right = True

        boundary_energy_left = False
        if qn_bias_left is not None or transfer_coeff_left is not None or u_reference_left is not None:
            fout.write('  [boundary-energy-left]\n')
            fout.write('    type = BoundaryEnergy\n')
            fout.write("    execute_on = 'timestep_end final'\n")
            fout.write("    variable = 'u'     # bulk energy unknown variable\n")
            fout.write('    diffCoeff = ${replace diff_coeff}\n')
            if qn_bias_left is not None:
                fout.write("    diffNormalFluxBias = ${replace qn_bias_left}\n")
            if transfer_coeff_left is not None:
                fout.write('    transferCoeff = ${replace transfer_coeff_left}\n')
            if u_reference_left is not None:
                fout.write('    reference = ${replace u_reference_left}\n')
            fout.write('    boundary = left\n')
            fout.write('  []\n')

            boundary_energy_left = True

        if boundary_energy_left and boundary_energy_right:

            fout.write('  [total-energy]\n')
            fout.write('    type = LinearCombinationPostprocessor\n')
            fout.write("    pp_names = 'bulk-energy boundary-energy-right boundary-energy-left'\n")
            fout.write("    pp_coefs = '1  1  1'\n")
            fout.write('    b = 0\n')
            fout.write('  []\n')

        elif boundary_energy_right:

            fout.write('  [total-energy]\n')
            fout.write('    type = LinearCombinationPostprocessor\n')
            fout.write("    pp_names = 'bulk-energy boundary-energy-right'\n")
            fout.write("    pp_coefs = '1  1'\n")
            fout.write('    b = 0\n')
            fout.write('  []\n')

        elif boundary_energy_left:

            fout.write('  [total-energy]\n')
            fout.write('    type = LinearCombinationPostprocessor\n')
            fout.write("    pp_names = 'bulk-energy boundary-energy-left'\n")
            fout.write("    pp_coefs = '1  1'\n")
            fout.write('    b = 0\n')
            fout.write('  []\n')

        elif not boundary_energy_left and not boundary_energy_right:

            pass

        else:

            assert False, 'Houston help... '

        fout.write('[]\n')

    fout.write('\n')

    fout.write('[VectorPostprocessors]\n')

    fout.write('  [omega-data]\n')
    fout.write('    type = LineValueSampler\n')
    fout.write("    execute_on = 'timestep_end final'\n")
    if compute_diffusion_flux:
        if u2_left is None and u2_right is None:
            fout.write("    variable = 'u diffFluxU_x'    # output data\n")
        else:
            fout.write("    variable = 'u u2 diffFluxU_x diffFluxU2_x'    # output data\n")
    else:
        if u2_left is None and u2_right is None:
            fout.write("    variable = 'u'    # output data\n")
        else:
            fout.write("    variable = 'u u2'    # output data\n")
    fout.write("    start_point = '${replace xmin} 0 0'\n")
    fout.write("    end_point = '${replace xmax} 0 0'\n")
    fout.write('    num_points = %i\n'%n_plot_pts)
    fout.write('    sort_by = id\n')
    fout.write('  []\n')

    fout.write('[]\n')

    fout.write('\n')

    fout.write('[Outputs]\n')
    fout.write('  console = true\n')
    fout.write('  [csv]\n')
    fout.write('    type = CSV\n')
    fout.write("    file_base = 'output'\n")
    fout.write("    execute_on = 'final'\n")
    #fout.write("    show = 'x-data'\n")
    fout.write('  []\n')
    if compute_energy:
        fout.write('  [console-energy]\n')
        fout.write('    type = Console\n')
        fout.write("    execute_on = 'final linear nonlinear'\n")
        if not boundary_energy_right and not boundary_energy_left:
          fout.write("    show = 'bulk-energy'\n")
        elif boundary_energy_right and not boundary_energy_left:
            fout.write("    show = 'bulk-energy boundary-energy-right total-energy'\n")
        elif not boundary_energy_right and boundary_energy_left:
            fout.write("    show = 'bulk-energy boundary-energy-left total-energy'\n")
        elif boundary_energy_right and boundary_energy_left:
            fout.write("    show = 'bulk-energy boundary-energy-left boundary-energy-right total-energy'\n")
        else:
            assert False, 'Houston help.'
        fout.write('  []\n')
        fout.write('  [file-energy]\n')
        fout.write('    type = CSV\n')
        fout.write("    execute_on = 'final'\n")
        fout.write("    file_base = 'output_energy'\n")
        if qn_bias_left is None and transfer_coeff_left is None and u_reference_left is None:
            fout.write("    show = 'bulk-energy'\n")
        else:
            fout.write("    show = 'bulk-energy boundary-energy total-energy'\n")
        fout.write('  []\n')
    fout.write('[]\n')

    fout.close()
#*********************************************************************************
def write_engy5310_1d_interfacial_coupling_input_file(
                                    x_left, x_right,
                                    x_interface=None,
                                    u1_left=None,
                                    qn_bias_left=None, transfer_coeff_left=None,
                                    u_reference_left=None,
                                    qn_bias_right=None, transfer_coeff_right=None,
                                    u_reference_right=None,
                                    diff_coeff_1=1.0,
                                    source_s_1=None,
                                    source_transfer_coeff=None,
                                    source_saturation=None,
                                    u2_right=None,
                                    diff_coeff_2=None,
                                    source_s_2=None,
                                    source_transfer_coeff_2=None,
                                    source_saturation_2=None,
                                    jump_coeff=1,
                                    is_jump_condition=True,
                                    interface_thickness_percent=0.1,
                                    is_flux_continuity=True,
                                    velocity=None,
                                    n_felem_1=1, x_bias_1=None,
                                    n_felem_2=1, x_bias_2=None,
                                    order='second',
                                    n_plot_pts_1=10,
                                    n_plot_pts_2=10,
                                    compute_diffusion_flux=False,
                                    use_moose_diffusion_flux=False,
                                    compute_energy=False,
                                    use_moose_neumann_bc=False,
                                    solver='pjfnk-hypre-boom',
                                    fdp_dh=1e-6,
                                    file_name=None):

    assert n_felem_1 > 1, 'Must use more than 1 FE.'
    assert n_felem_2 > 1, 'Must use more than 1 FE.'

    if x_interface is None:
        x_interface = (x_left + x_right)/2.0

    interf_thickness = (x_right - x_left)*interface_thickness_percent/100

    if qn_bias_left is not None:
        assert u1_left is None
    if transfer_coeff_left is not None:
        assert transfer_coeff_left >= 0.0
        assert u1_left is None
    if u_reference_left is not None:
        assert u1_left is None

    if velocity is not None:
        assert isinstance(velocity, tuple)
        assert len(velocity) == 3

    if x_bias_1 is not None:
        assert 0.5 <= x_bias_1 <= 2
    if x_bias_2 is not None:
        assert 0.5 <= x_bias_2 <= 2

    if file_name is None:
        fout = open('engy5310p1/input.hit', 'wt')
    else:
        fout = open(file_name, 'wt')

    fout.write('# Engy-5310 Problem 1: Poisson 1D FEM\n')
    fout.write('# UMass Lowell Nuclear Chemical Engineering\n')
    fout.write('# Prof. Valmor F. de Almeida\n')
    time_date_stamp = datetime.datetime.today().strftime('%d%b%y %H:%M:%S')
    fout.write('# %s\n'%time_date_stamp)

    fout.write('\n')

    fout.write('# Parameters\n')
    fout.write('xmin = %10.5e\n'%x_left)
    fout.write('xmax = %10.5e\n'%x_right)
    fout.write('x_interface = %10.5e\n'%x_interface)

    fout.write('\n')

    if u1_left is not None:
        fout.write('u1_left = %10.5e\n'%u1_left)

    if diff_coeff_1 is not None:
        fout.write('diff_coeff_1 = %10.5e\n'%diff_coeff_1)
    if source_s_1 is not None:
        fout.write('source_s_1 = %10.5e\n'%source_s_1)
    if source_transfer_coeff is not None:
        fout.write('source_transfer_coeff = %10.5e\n'%source_transfer_coeff)
    if source_saturation is not None:
        fout.write('source_saturation = %10.5e\n'%source_saturation)

    fout.write('\n')

    if u2_right is not None:
        fout.write('u2_right = %10.5e\n'%u2_right)

    if diff_coeff_2 is not None:
        fout.write('diff_coeff_2 = %10.5e\n'%diff_coeff_2)
    if source_s_2 is not None:
        fout.write('source_s_2 = %10.5e\n'%source_s_2)
    if source_transfer_coeff_2 is not None:
        fout.write('source_transfer_coeff_2 = %10.5e\n'%source_transfer_coeff_2)
    if source_saturation_2 is not None:
        fout.write('source_saturation_2 = %10.5e\n'%source_saturation_2)

    fout.write('\n')

    if velocity is not None:
        fout.write("velocity = '%10.5e %10.5e %10.5e'\n"%(velocity[0],velocity[1],velocity[2]))

    if qn_bias_left is not None and use_moose_neumann_bc is False:
        fout.write('qn_bias_left = %10.5e\n'%qn_bias_left)
    if transfer_coeff_left is not None:
        fout.write('transfer_coeff_left = %10.5e\n'%transfer_coeff_left)
    if u_reference_left is not None:
        fout.write('u_reference_left = %10.5e\n'%u_reference_left)

    if qn_bias_right is not None and use_moose_neumann_bc is False:
        fout.write('qn_bias_right = %10.5e\n'%qn_bias_right)
    if transfer_coeff_right is not None:
        fout.write('transfer_coeff_right = %10.5e\n'%transfer_coeff_right)
    if u_reference_right is not None:
        fout.write('u_reference_right = %10.5e\n'%u_reference_right)

    fout.write('jump_coeff = %10.5e\n'%jump_coeff)

    fout.write('\n')

    fout.write('[Problem]\n')
    fout.write('  type = FEProblem\n')
    fout.write('  coord_type = XYZ\n')
    fout.write('[]\n')

    fout.write('\n')

    fout.write('[Mesh]\n')

    fout.write('  [omega1]\n')
    fout.write('    type = GeneratedMeshGenerator\n')
    fout.write('    dim = 1\n')
    fout.write('    xmin = ${replace xmin}\n')
    fout.write('    xmax = ${replace x_interface}\n')
    fout.write('    nx = %i\n'%n_felem_1)
    if order == 'second':
        fout.write('    elem_type = edge3\n')
    if x_bias_2 is not None:
        fout.write('    bias_x = %3.3e\n'%x_bias_2)
    fout.write('  []\n')

    fout.write('  [omega2]\n')
    fout.write('    type = GeneratedMeshGenerator\n')
    fout.write('    dim = 1\n')
    fout.write('    xmin = ${replace x_interface}\n')
    fout.write('    xmax = ${replace xmax}\n')
    fout.write('    nx = %i\n'%n_felem_2)
    if order == 'second':
        fout.write('    elem_type = edge3\n')
    if x_bias_1 is not None:
        fout.write('    bias_x = %3.3e\n'%x_bias_1)
    fout.write('  []\n')

    fout.write('  [omega]\n')
    fout.write('    type = StitchedMeshGenerator\n')
    fout.write("    inputs = 'omega1 omega2'\n")
    fout.write("    stitch_boundaries_pairs = 'right left'\n")
    fout.write("    clear_stitched_boundary_ids = 'true'\n")
    fout.write('  []\n')

    fout.write('# Create subdomains: Omega_1 and Omega_2\n')
    fout.write('  [mod1]\n')
    fout.write('    type = SubdomainBoundingBoxGenerator\n')
    fout.write('    input = omega\n')
    fout.write('    block_id = 1\n')
    fout.write('    block_name = omega_1\n')
    fout.write("    bottom_left = '${replace xmin} 0 0'\n")
    fout.write("    top_right = '${replace x_interface} 1 0'\n")
    fout.write('  []\n')

    fout.write('  [mod2]\n')
    fout.write('    type = SubdomainBoundingBoxGenerator\n')
    fout.write('    input = mod1\n')
    fout.write('    block_id = 2\n')
    fout.write('    block_name = omega_2\n')
    fout.write("    bottom_left = '${replace x_interface} 0 0'\n")
    fout.write("    top_right = '${replace xmax} 1 0'\n")
    fout.write('  []\n')

    fout.write('# Create interface of subdomains: Omega_1 and Omega_2\n')
    fout.write('  [mod3]\n')
    fout.write('    type = SideSetsBetweenSubdomainsGenerator\n')
    fout.write('    input = mod2\n')
    fout.write('    primary_block = omega_1\n')
    fout.write('    paired_block = omega_2 \n')
    fout.write('    new_boundary = interface_12\n')
    fout.write('  []\n')

    fout.write('# Create boundaries of subdomains: Omega_1 and Omega_2\n')
    fout.write('  [mod4]\n')
    fout.write('    type = SideSetsAroundSubdomainGenerator\n')
    fout.write('    input = mod3\n')
    fout.write('    block = omega_1\n')
    fout.write("    normal = '-1 0 0'\n")
    fout.write('    new_boundary = omega_1_left\n')
    fout.write('  []\n')

    fout.write('  [mod5]\n')
    fout.write('    type = SideSetsAroundSubdomainGenerator\n')
    fout.write('    input = mod4\n')
    fout.write('    block = omega_2\n')
    fout.write("    normal = '1 0 0'\n")
    fout.write('    new_boundary = omega_2_right\n')
    fout.write('  []\n')

    fout.write('[]\n')

    fout.write('\n')

    fout.write('[Variables]\n')

    fout.write('  [u1]\n')
    fout.write('    block = omega_1\n')
    fout.write('    order = %s\n'%order)
    fout.write('    family = lagrange\n')
    if u1_left is not None:
        fout.write('    initial_condition = ${replace u1_left}\n')
    fout.write('  []\n')

    fout.write('  [u2]\n')
    fout.write('    block = omega_2\n')
    fout.write('    order = %s\n'%order)
    fout.write('    family = lagrange\n')
    if u2_right is not None:
        fout.write('    initial_condition = ${replace u2_right}\n')
    fout.write('  []\n')

    fout.write('[]\n')

    fout.write('\n')

    if compute_diffusion_flux:
        fout.write('[AuxVariables]\n')

        fout.write('  [diffFluxU1]\n')
        if order in ['second', 'SECOND']:
            fout.write('    order = FIRST\n')
        elif order in ['first', 'FIRST']:
            fout.write('    order = CONSTANT\n')
        else:
            assert False
        fout.write('    family = MONOMIAL_VEC\n')
        fout.write('  []\n')

        fout.write('  [diffFluxU2]\n')
        if order in ['second', 'SECOND']:
            fout.write('    order = FIRST\n')
        elif order in ['first', 'FIRST']:
            fout.write('    order = CONSTANT\n')
        else:
            assert False
        fout.write('    family = MONOMIAL_VEC\n')
        fout.write('  []\n')

        fout.write('  [diffFluxU1_x]\n')
        if order in ['second', 'SECOND']:
            fout.write('    order = FIRST\n')
        elif order in ['first', 'FIRST']:
            fout.write('    order = CONSTANT\n')
        else:
            assert False
        fout.write('    family = MONOMIAL\n')
        fout.write('  []\n')

        fout.write('  [diffFluxU2_x]\n')
        if order in ['second', 'SECOND']:
            fout.write('    order = FIRST\n')
        elif order in ['first', 'FIRST']:
            fout.write('    order = CONSTANT\n')
        else:
            assert False
        fout.write('    family = MONOMIAL\n')
        fout.write('  []\n')

        fout.write('[]\n')

        fout.write('\n')

    fout.write('[Kernels]\n')

    fout.write('  [diffusion-term-1]\n')
    fout.write('    block = omega_1\n')
    fout.write('    type = DiffusionTerm\n')
    fout.write('    variable = u1     # produced quantity\n')
    fout.write('    diffCoeff = ${replace diff_coeff_1}\n')
    fout.write('  []\n')

    if source_s_1 is not None:
        fout.write('  [source-term-1]\n')
        fout.write('    type = SourceTerm\n')
        fout.write('    block = omega_1\n')
        fout.write('    variable = u1     # add to produced quantity\n')
        fout.write('    sourceS = ${replace source_s_1}\n')
        if source_transfer_coeff is not None:
            fout.write('    transferCoeff = ${replace source_transfer_coeff}\n')
        if source_saturation is not None:
            fout.write('    saturation = ${replace source_saturation}\n')
        fout.write('    coupledVariable = u1\n')
        fout.write('  []\n')

    if velocity is not None:
        fout.write('  [convection-term-1]\n')
        fout.write('    type = ConvectionTerm\n')
        fout.write('    block = omega_1\n')
        fout.write('    variable = u1     # produced quantity\n')
        fout.write('    velocity = ${replace velocity}\n')
        fout.write('  []\n')


    fout.write('  [diffusion-term-2]\n')
    fout.write('    type = DiffusionTerm\n')
    fout.write('    block = omega_2\n')
    fout.write('    variable = u2     # produced quantity\n')
    fout.write('    diffCoeff = ${replace diff_coeff_2}\n')
    fout.write('  []\n')

    if source_s_2 is not None:
        fout.write('  [source-term-2]\n')
        fout.write('    type = SourceTerm\n')
        fout.write('    block = omega_2\n')
        fout.write('    variable = u2     # add to produced quantity\n')
        fout.write('    sourceS = ${replace source_s_2}\n')
        if source_transfer_coeff_2 is not None:
            fout.write('    transferCoeff = ${replace source_transfer_coeff_2}\n')
        if source_saturation_2 is not None:
            fout.write('    saturation = ${replace source_saturation_2}\n')
        fout.write('    coupledVariable = u2\n')
        fout.write('  []\n')

    if velocity is not None:
        fout.write('  [convection-term-2]\n')
        fout.write('    block = omega_2\n')
        fout.write('    type = ConvectionTerm\n')
        fout.write('    variable = u2     # produced quantity\n')
        fout.write('    velocity = ${replace velocity}\n')
        fout.write('  []\n')

    fout.write('[]\n')

    fout.write('\n')

    if is_jump_condition==True or is_flux_continuity==True:
        fout.write('[InterfaceKernels]\n')
        if is_jump_condition==True:
            fout.write('  [jump]\n')
            fout.write('    type = InterfaceJump\n')
            fout.write('    variable = u1\n')
            fout.write('    neighbor_var = u2\n')
            fout.write('    boundary = interface_12\n')
            fout.write('    jumpCoeff = ${replace jump_coeff} # jump coefficient u1 = k * u2\n')
            fout.write('  []\n')

        if is_flux_continuity==True:
            fout.write('  [normal-flux-continuity]\n')
            fout.write('    type = InterfaceNormalFluxContinuity\n')
            fout.write('    variable = u1\n')
            fout.write('    neighbor_var = u2\n')
            fout.write('    boundary = interface_12\n')
            fout.write('    diffCoeff = ${replace diff_coeff_1}\n')
            fout.write('    diffCoeffNeighbor = ${replace diff_coeff_2}\n')
            fout.write('  []\n')

        fout.write('[]\n')
        fout.write('\n')

    if compute_diffusion_flux:
        fout.write('[AuxKernels]\n')

        fout.write('  [diffusion-flux-1]\n')
        fout.write('    execute_on = timestep_end\n')
        if use_moose_diffusion_flux:
            fout.write('    type = DiffusionFluxAux     # provided by MOOSE\n')
            fout.write('    diffusion_variable = u1\n')
            fout.write('    block = omega_1\n')
            fout.write('    diffusivity = ${replace diff_coeff_1}\n')
            fout.write('    component = x\n')
            fout.write('    variable = diffFluxU1_x     # produced quantity\n')
            fout.write('  []\n')
        else:
            fout.write('    type = DiffusionFlux\n')
            fout.write('    field = u1\n')
            fout.write('    block = omega_1\n')
            fout.write('    diffCoeff = ${replace diff_coeff_1}\n')
            fout.write('    variable = diffFluxU1     # produced quantity\n')
            fout.write('  []\n')
            fout.write('  [diffusion-flux-x-1]\n')
            fout.write('    execute_on = timestep_end\n')
            fout.write('    type = VectorVariableComponentAux\n')
            fout.write('    variable = diffFluxU1_x    # produced quantity\n')
            fout.write('    block = omega_1\n')
            fout.write('    component = x\n')
            fout.write('    vector_variable = diffFluxU1   \n')
            fout.write('  []\n')

        fout.write('  [diffusion-flux-2]\n')
        fout.write('    execute_on = timestep_end\n')
        if use_moose_diffusion_flux:
            fout.write('    type = DiffusionFluxAux     # provided by MOOSE\n')
            fout.write('    diffusion_variable = u2\n')
            fout.write('    block = omega_2\n')
            fout.write('    diffusivity = ${replace diff_coeff_2}\n')
            fout.write('    component = x\n')
            fout.write('    variable = diffFluxU2_x     # produced quantity\n')
            fout.write('  []\n')
        else:
            fout.write('    type = DiffusionFlux\n')
            fout.write('    field = u2\n')
            fout.write('    block = omega_2\n')
            fout.write('    diffCoeff = ${replace diff_coeff_2}\n')
            fout.write('    variable = diffFluxU2     # produced quantity\n')
            fout.write('  []\n')
            fout.write('  [diffusion-flux-x-2]\n')
            fout.write('    execute_on = timestep_end\n')
            fout.write('    type = VectorVariableComponentAux\n')
            fout.write('    variable = diffFluxU2_x    # produced quantity\n')
            fout.write('    block = omega_2\n')
            fout.write('    component = x\n')
            fout.write('    vector_variable = diffFluxU2   \n')
            fout.write('  []\n')

        fout.write('[]\n')
        fout.write('\n')

    fout.write('[BCs]\n')

    if u1_left is not None:
        fout.write('  [left]\n')
        fout.write('    type = DirichletBC\n')
        fout.write('    variable = u1\n')
        fout.write('    boundary = omega_1_left\n')
        fout.write('    value = ${replace u1_left}\n')
        fout.write('  []\n')

    if u2_right is not None:
        fout.write('  [right]\n')
        fout.write('    type = DirichletBC\n')
        fout.write('    variable = u2\n')
        fout.write('    boundary = omega_2_right\n')
        fout.write('    value = ${replace u2_right}\n')
        fout.write('  []\n')

    if u1_left is None and use_moose_neumann_bc is False:
        fout.write('  [left]\n')
        fout.write('    type = NormalFluxBC\n')
        fout.write('    variable = u1\n')
        fout.write('    boundary = omega_1_left\n')
        if qn_bias_left is not None:
           fout.write('    bias = ${replace qn_bias_left}\n')
        if transfer_coeff_left is not None:
            fout.write('    transferCoeff = ${replace transfer_coeff_left}\n')
        if u_reference_left is not None:
            fout.write('    reference = ${replace u_reference_left}\n')
        fout.write('  []\n')

    if u1_left is None and use_moose_neumann_bc is True:
        fout.write('  [left]\n')
        fout.write('    type = NeumannBC\n')
        fout.write('    variable = u1\n')
        fout.write('    boundary = omega_1_left\n')
        if qn_bias_left is not None:
            fout.write('    value = %10.5e\n'%(qn_bias_left))
        else:
            fout.write('      value = 0.0\n')
        fout.write('  []\n')

    fout.write('[]\n')

    fout.write('\n')

    if solver=='fdp-newt-full':
        fout.write('[Preconditioning]\n')
        fout.write("  active = 'fdp-newt-full'\n")
        fout.write('  [fdp-newt-full]\n')
        fout.write('    type = FDP\n')
        fout.write('    full = true\n')
        fout.write("    solve_type = 'NEWTON'\n")
        fout.write("    petsc_options_iname = '-pc_type -mat_fd_coloring_err -mat_fd_type'\n")
        #fout.write("    petsc_options_value = 'lu       1e-6                 ds'\n")
        fout.write("    petsc_options_value = 'lu  %10.3e               ds'\n"%fdp_dh)
        fout.write('  []\n')
        fout.write('[]\n')

        fout.write('\n')

        fout.write('[Executioner]\n')
        fout.write('  type = Steady\n')
        fout.write('[]\n')
    elif solver=='pjfnk-hypre-boom':
        fout.write('[Executioner]\n')
        fout.write('  type = Steady\n')
        fout.write("  solve_type = 'PJFNK'\n")
        fout.write("  petsc_options_iname = '-pc_type -pc_hypre_type'\n")
        fout.write("  petsc_options_value = 'hypre boomeramg'\n")
        fout.write('[]\n')
    else:
        assert False,'Invalid solver option = %r'%solver


    if compute_energy:

        fout.write('\n')

        fout.write('[Postprocessors]\n')
        fout.write('  [bulk-energy]\n')
        fout.write('    type = BulkEnergy\n')
        fout.write("    execute_on = 'timestep_end final'\n")
        fout.write("    variable = 'u1'     # bulk energy unknown variable\n")
        fout.write('    diffCoeff = ${replace diff_coeff_1}\n')
        fout.write('    sourceS = ${replace source_s}\n')
        fout.write('  []\n')

        boundary_energy_right = False
        if qn_bias_right is not None or transfer_coeff_right is not None or u_reference_right is not None:
            fout.write('  [boundary-energy-right]\n')
            fout.write('    type = BoundaryEnergy\n')
            fout.write("    execute_on = 'timestep_end final'\n")
            fout.write("    variable = 'u1'     # bulk energy unknown variable\n")
            fout.write('    diffCoeff = ${replace diff_coeff_1}\n')
            if qn_bias_right is not None:
                fout.write("    diffNormalFluxBias = ${replace qn_bias_right}\n")
            if transfer_coeff_right is not None:
                fout.write('    transferCoeff = ${replace transfer_coeff_right}\n')
            if u_reference_right is not None:
                fout.write('    reference = ${replace u_reference_right}\n')
            fout.write('    boundary = right\n')
            fout.write('  []\n')

            boundary_energy_right = True

        boundary_energy_left = False
        if qn_bias_left is not None or transfer_coeff_left is not None or u_reference_left is not None:
            fout.write('  [boundary-energy-left]\n')
            fout.write('    type = BoundaryEnergy\n')
            fout.write("    execute_on = 'timestep_end final'\n")
            fout.write("    variable = 'u1'     # bulk energy unknown variable\n")
            fout.write('    diffCoeff = ${replace diff_coeff_1}\n')
            if qn_bias_left is not None:
                fout.write("    diffNormalFluxBias = ${replace qn_bias_left}\n")
            if transfer_coeff_left is not None:
                fout.write('    transferCoeff = ${replace transfer_coeff_left}\n')
            if u_reference_left is not None:
                fout.write('    reference = ${replace u_reference_left}\n')
            fout.write('    boundary = left\n')
            fout.write('  []\n')

            boundary_energy_left = True

        if boundary_energy_left and boundary_energy_right:

            fout.write('  [total-energy]\n')
            fout.write('    type = LinearCombinationPostprocessor\n')
            fout.write("    pp_names = 'bulk-energy boundary-energy-right boundary-energy-left'\n")
            fout.write("    pp_coefs = '1  1  1'\n")
            fout.write('    b = 0\n')
            fout.write('  []\n')

        elif boundary_energy_right:

            fout.write('  [total-energy]\n')
            fout.write('    type = LinearCombinationPostprocessor\n')
            fout.write("    pp_names = 'bulk-energy boundary-energy-right'\n")
            fout.write("    pp_coefs = '1  1'\n")
            fout.write('    b = 0\n')
            fout.write('  []\n')

        elif boundary_energy_left:

            fout.write('  [total-energy]\n')
            fout.write('    type = LinearCombinationPostprocessor\n')
            fout.write("    pp_names = 'bulk-energy boundary-energy-left'\n")
            fout.write("    pp_coefs = '1  1'\n")
            fout.write('    b = 0\n')
            fout.write('  []\n')

        elif not boundary_energy_left and not boundary_energy_right:

            pass

        else:

            assert False, 'Houston help... '

        fout.write('[]\n')

    fout.write('\n')

    fout.write('[VectorPostprocessors]\n')

    fout.write('  [omega_1]\n')
    fout.write('    type = LineValueSampler\n')
    fout.write("    execute_on = 'timestep_end final'\n")
    if compute_diffusion_flux:
        fout.write("    variable = 'u1 diffFluxU1_x'    # output data\n")
    else:
        fout.write("    variable = 'u1'    # output data\n")
    fout.write("    start_point = '${replace xmin} 0 0'\n")
    fout.write("    end_point = '${fparse x_interface-%5.2e} 0 0'\n"%interf_thickness)
    fout.write('    num_points = %i\n'%n_plot_pts_1)
    fout.write('    sort_by = id\n')
    fout.write('  []\n')

    fout.write('  [omega_2]\n')
    fout.write('    type = LineValueSampler\n')
    fout.write("    execute_on = 'timestep_end final'\n")
    if compute_diffusion_flux:
        fout.write("    variable = 'u2 diffFluxU2_x'    # output data\n")
    else:
        fout.write("    variable = 'u2'    # output data\n")
    fout.write("    start_point = '${fparse x_interface+%5.2e} 0 0'\n"%interf_thickness)
    fout.write("    end_point = '${replace xmax} 0 0'\n")
    fout.write('    num_points = %i\n'%n_plot_pts_2)
    fout.write('    sort_by = id\n')
    fout.write('  []\n')

    fout.write('[]\n')

    fout.write('\n')

    fout.write('[Outputs]\n')
    fout.write('  console = true\n')
    fout.write('  [csv]\n')
    fout.write('    type = CSV\n')
    fout.write("    file_base = 'output'\n")
    fout.write("    execute_on = 'final'\n")
    #fout.write("    show = 'x-data'\n")
    fout.write('  []\n')
    if compute_energy:
        fout.write('  [console-energy]\n')
        fout.write('    type = Console\n')
        fout.write("    execute_on = 'final linear nonlinear'\n")
        if not boundary_energy_right and not boundary_energy_left:
          fout.write("    show = 'bulk-energy'\n")
        elif boundary_energy_right and not boundary_energy_left:
            fout.write("    show = 'bulk-energy boundary-energy-right total-energy'\n")
        elif not boundary_energy_right and boundary_energy_left:
            fout.write("    show = 'bulk-energy boundary-energy-left total-energy'\n")
        elif boundary_energy_right and boundary_energy_left:
            fout.write("    show = 'bulk-energy boundary-energy-left boundary-energy-right total-energy'\n")
        else:
            assert False, 'Houston help.'
        fout.write('  []\n')
        fout.write('  [file-energy]\n')
        fout.write('    type = CSV\n')
        fout.write("    execute_on = 'final'\n")
        fout.write("    file_base = 'output_energy'\n")
        if qn_bias_left is None and transfer_coeff_left is None and u_reference_left is None:
            fout.write("    show = 'bulk-energy'\n")
        else:
            fout.write("    show = 'bulk-energy boundary-energy total-energy'\n")
        fout.write('  []\n')
    fout.write('[]\n')

    fout.close()
#*********************************************************************************
def write_engy5310_p1_2d_input_file(rz_axis=None,
                                 x_left=None, x_right=None, y_bottom=None, y_top=None,
                                 u_left=None, u_right=None, u_bottom=None, u_top=None,
                                 qn_bias_bottom=None, transfer_coeff_bottom=None,
                                 u_reference_bottom=None,
                                 qn_bias_top=None, transfer_coeff_top=None,
                                 u_reference_top=None,
                                 qn_bias_left=None, transfer_coeff_left=None,
                                 u_reference_left=None,
                                 qn_bias_right=None, transfer_coeff_right=None,
                                 u_reference_right=None,
                                 diff_coeff=1.0,
                                 source_s=1.0,
                                 source_transfer_coeff=None,
                                 source_saturation=None,
                                 diff_coeff2=1.0,
                                 u2_left=None, u2_right=None,
                                 diff_coeff_2=None,
                                 source_s_2=None,
                                 source_transfer_coeff_2=None,
                                 source_saturation_2=None,
                                 velocity=None,
                                 n_felem_x=1, n_felem_y=1,
                                 x_bias=None, y_bias=None,
                                 order='second',
                                 x_line_ypos=None, y_line_xpos=None,
                                 n_plot_pts_x=10, n_plot_pts_y=10,
                                 compute_diffusion_flux=False,
                                 use_moose_diffusion_flux=False,
                                 compute_energy=False,
                                 use_moose_neumann_bc=False,
                                 solver='pjfnk-hypre-boom',
                                 fdp_dh=1e-6,
                                 file_name=None):

    if u_right is not None:
        assert qn_bias_right is None
        assert transfer_coeff_right is None
        assert u_reference_right is None
    if qn_bias_right is not None:
        assert u_right is None
    if transfer_coeff_right is not None:
        assert transfer_coeff_right >= 0.0
        assert u_right is None
    if u_reference_right is not None:
        assert u_right is None

    if velocity is not None:
        assert isinstance(velocity, tuple)
        assert len(velocity) == 3

    if x_bias is not None:
        assert 0.5 <= x_bias <= 2
    if y_bias is not None:
        assert 0.5 <= y_bias <= 2

    if x_line_ypos is not None:
        assert y_bottom <= x_line_ypos <= y_top
    if y_line_xpos is not None:
        assert x_left <= y_line_xpos <= x_right

    if file_name is None:
        fout = open('engy5310p1/input.hit', 'wt')
    else:
        fout = open(file_name, 'wt')

    fout.write('# Engy-5310 Problem: Poisson 2D FEM\n')
    fout.write('# UMass Lowell Nuclear Chemical Engineering\n')
    fout.write('# Prof. Valmor F. de Almeida\n')
    time_date_stamp = datetime.datetime.today().strftime('%d%b%y %H:%M:%S')
    fout.write('# %s\n'%time_date_stamp)

    fout.write('\n')

    fout.write('# Parameters\n')
    fout.write('xmin = %10.5e\n'%x_left)
    fout.write('xmax = %10.5e\n'%x_right)
    fout.write('ymin = %10.5e\n'%y_bottom)
    fout.write('ymax = %10.5e\n'%y_top)
    fout.write('diff_coeff = %10.5e\n'%diff_coeff)
    if source_s is not None:
        fout.write('source_s = %10.5e\n'%source_s)
    if source_transfer_coeff is not None:
        fout.write('source_transfer_coeff = %10.5e\n'%source_transfer_coeff)
    if source_saturation is not None:
        fout.write('source_saturation = %10.5e\n'%source_saturation)

    if u_left is not None:
        fout.write('u_left = %10.5e\n'%u_left)
    if u_right is not None:
        fout.write('u_right = %10.5e\n'%u_right)
    if u_bottom is not None:
        fout.write('u_bottom = %10.5e\n'%u_bottom)
    if u_top is not None:
        fout.write('u_top = %10.5e\n'%u_top)

    if velocity is not None:
        fout.write("velocity = '%10.5e %10.5e %10.5e'\n"%(velocity[0],velocity[1],velocity[2]))


    if qn_bias_left is not None:
        fout.write('qn_bias_left = %10.5e\n'%qn_bias_left)
    if transfer_coeff_left is not None:
        fout.write('transfer_coeff_left = %10.5e\n'%transfer_coeff_left)
    if u_reference_left is not None:
        fout.write('u_reference_left = %10.5e\n'%u_reference_left)

    if qn_bias_right is not None:
        fout.write('qn_bias_right = %10.5e\n'%qn_bias_right)
    if transfer_coeff_right is not None:
        fout.write('transfer_coeff_right = %10.5e\n'%transfer_coeff_right)
    if u_reference_right is not None:
        fout.write('u_reference_right = %10.5e\n'%u_reference_right)

    if qn_bias_bottom is not None:
        fout.write('qn_bias_bottom = %10.5e\n'%qn_bias_bottom)
    if transfer_coeff_bottom is not None:
        fout.write('transfer_coeff_bottom_= %10.5e\n'%transfer_coeff_bottom)
    if u_reference_bottom is not None:
        fout.write('u_reference_bottom = %10.5e\n'%u_reference_bottom)

    if qn_bias_top is not None:
        fout.write('qn_bias_top = %10.5e\n'%qn_bias_top)
    if transfer_coeff_top is not None:
        fout.write('transfer_coeff_top_= %10.5e\n'%transfer_coeff_top)
    if u_reference_top is not None:
        fout.write('u_reference_top = %10.5e\n'%u_reference_top)

    fout.write('\n')

    fout.write('[Problem]\n')
    fout.write('  type = FEProblem\n')
    if rz_axis is None:
        fout.write('  coord_type = XYZ\n')
    elif rz_axis == 'y':
        fout.write('  coord_type = RZ\n')
        fout.write('  rz_coord_axis = Y\n')
    else:
        assert False
    fout.write('[]\n')

    fout.write('\n')

    fout.write('[Mesh]\n')
    fout.write('  [2d]\n')
    fout.write('    type = GeneratedMeshGenerator\n')
    fout.write('    dim = 2\n')
    fout.write('    xmin = ${replace xmin}\n')
    fout.write('    xmax = ${replace xmax}\n')
    fout.write('    ymin = ${replace ymin}\n')
    fout.write('    ymax = ${replace ymax}\n')
    fout.write('    nx = %i\n'%n_felem_x)
    fout.write('    ny = %i\n'%n_felem_y)
    if order == 'second':
        fout.write('    elem_type = QUAD9\n')
    if x_bias is not None:
        fout.write('    bias_x = %3.3e\n'%x_bias)
    if y_bias is not None:
        fout.write('    bias_y = %3.3e\n'%y_bias)
    fout.write('  []\n')
    fout.write('[]\n')

    fout.write('\n')

    fout.write('[Variables]\n')
    fout.write('  [u]\n')
    fout.write('    order = %s\n'%order)
    fout.write('    family = lagrange\n')
    values = list()
    values = [u_left, u_right, u_bottom, u_top]
    names = ['u_left', 'u_right', 'u_bottom', 'u_top']
    ids = [i for (i,v) in enumerate(values) if v is not None]
    if len(ids):
        average = str()
        for i in ids:
            average += names[i]+'+'
        average = '('+average[:-1]+')/'+str(len(ids))
        fout.write('    initial_condition = ${fparse '+average+'}\n')
    fout.write('  []\n')
    fout.write('[]\n')

    fout.write('\n')

    if compute_diffusion_flux:
        fout.write('[AuxVariables]\n')

        fout.write('  [diffFluxU]\n')
        if order in ['second', 'SECOND']:
            fout.write('    order = FIRST\n')
        elif order in ['first', 'FIRST']:
            fout.write('    order = CONSTANT\n')
        else:
            assert False
        fout.write('    family = MONOMIAL_VEC\n')
        fout.write('  []\n')

        fout.write('  [diffFluxU_x]\n')
        if order in ['second', 'SECOND']:
            fout.write('    order = FIRST\n')
        elif order in ['first', 'FIRST']:
            fout.write('    order = CONSTANT\n')
        else:
            assert False
        fout.write('    family = MONOMIAL\n')
        fout.write('  []\n')

        fout.write('  [diffFluxU_y]\n')
        if order in ['second', 'SECOND']:
            fout.write('    order = FIRST\n')
        elif order in ['first', 'FIRST']:
            fout.write('    order = CONSTANT\n')
        else:
            assert False
        fout.write('    family = MONOMIAL\n')
        fout.write('  []\n')

        fout.write('[]\n')

        fout.write('\n')

    fout.write('[Kernels]\n')
    fout.write('  [diffusion-term]\n')
    fout.write('    type = DiffusionTerm\n')
    fout.write('    variable = u     # produced quantity\n')
    fout.write('    diffCoeff = ${replace diff_coeff}\n')
    fout.write('  []\n')

    if source_s is not None:
        fout.write('  [source-term]\n')
        fout.write('    type = SourceTerm\n')
        fout.write('    variable = u     # add to produced quantity\n')
        if source_s is not None:
            fout.write('    sourceS = ${replace source_s}\n')
        if source_transfer_coeff is not None:
            fout.write('    transferCoeff = ${replace source_transfer_coeff}\n')
        if source_saturation is not None:
            fout.write('    saturation = ${replace source_saturation}\n')
        if u2_left is None and u2_right is None:
            fout.write('    coupledVariable = u\n')
        else:
            fout.write('    coupledVariable = u2\n')
        if source_s_2 is not None:
            fout.write('    sourceSCoupled = ${replace source_s_2}\n')
        if source_transfer_coeff_2 is not None:
            fout.write('    transferCoeffCoupled = ${replace source_transfer_coeff_2}\n')
        if source_saturation_2 is not None:
            fout.write('    saturationCoupled = ${replace source_saturation_2}\n')
        fout.write('  []\n')

    if velocity is not None:
        fout.write('  [convection-term]\n')
        fout.write('    type = ConvectionTerm\n')
        fout.write('    variable = u     # produced quantity\n')
        fout.write('    velocity = ${replace velocity}\n')
        fout.write('  []\n')
    fout.write('[]\n')

    fout.write('\n')

    if compute_diffusion_flux:
        fout.write('[AuxKernels]\n')
        if use_moose_diffusion_flux:
            fout.write('  [diffusion-flux-x]\n')
            fout.write('    execute_on = timestep_end\n')
            fout.write('    type = DiffusionFluxAux     # provided by MOOSE\n')
            fout.write('    diffusion_variable = u\n')
            fout.write('    diffusivity = ${replace diff_coeff}\n')
            fout.write('    component = x\n')
            fout.write('    variable = diffFluxU_x     # produced quantity\n')
            fout.write('  []\n')
            fout.write('  [diffusion-flux-y]\n')
            fout.write('    execute_on = timestep_end\n')
            fout.write('    type = DiffusionFluxAux     # provided by MOOSE\n')
            fout.write('    diffusion_variable = u\n')
            fout.write('    diffusivity = ${replace diff_coeff}\n')
            fout.write('    component = y\n')
            fout.write('    variable = diffFluxU_y     # produced quantity\n')
            fout.write('  []\n')
        else:
            fout.write('  [diffusion-flux]\n')
            fout.write('    execute_on = timestep_end\n')
            fout.write('    type = DiffusionFlux   # new kernel\n')
            fout.write('    field = u\n')
            fout.write('    diffCoeff = ${replace diff_coeff}\n')
            fout.write('    variable = diffFluxU     # produced quantity\n')
            fout.write('  []\n')
            fout.write('  [diffusion-flux-x]\n')
            fout.write('    execute_on = timestep_end\n')
            fout.write('    type = VectorVariableComponentAux  # provided by MOOSE\n')
            fout.write('    variable = diffFluxU_x    # produced quantity\n')
            fout.write('    component = x\n')
            fout.write('    vector_variable = diffFluxU   \n')
            fout.write('  []\n')
            fout.write('  [diffusion-flux-y]\n')
            fout.write('    execute_on = timestep_end\n')
            fout.write('    type = VectorVariableComponentAux  # provided by MOOSE\n')
            fout.write('    variable = diffFluxU_y    # produced quantity\n')
            fout.write('    component = y\n')
            fout.write('    vector_variable = diffFluxU   \n')
            fout.write('  []\n')
        fout.write('[]\n')
        fout.write('\n')

    fout.write('[BCs]\n')

    fout.write('  [east]\n')
    if u_left is None and use_moose_neumann_bc:
        fout.write('    type = NeumannBC\n')
        fout.write('    variable = u\n')
        fout.write('    boundary = left\n')
        if qn_bias_left is not None:
            fout.write('    value = ${replace qn_bias_left}\n')
        else:
           fout.write('    value = %10.5e\n'%0.0)
    elif u_left is None and use_moose_neumann_bc is False:
        fout.write('    type = NormalFluxBC\n')
        fout.write('    variable = u\n')
        fout.write('    boundary = left\n')
        if qn_bias_left is not None:
            fout.write('    bias = ${replace qn_bias_left}\n')
        if transfer_coeff_left is not None:
            fout.write('    transferCoeff = ${replace transfer_coeff_left}\n')
        if u_reference_left is not None:
            fout.write('    reference = ${replace u_reference_left}\n')
    else:
        fout.write('    type = DirichletBC\n')
        fout.write('    variable = u\n')
        fout.write('    boundary = left\n')
        fout.write('    value = ${replace u_left}\n')
    fout.write('  []\n')

    fout.write('  [west]\n')
    if u_right is None and use_moose_neumann_bc:
        fout.write('    type = NeumannBC\n')
        fout.write('    variable = u\n')
        fout.write('    boundary = right\n')
        if qn_bias_right is not None:
            fout.write('    value = ${replace qn_bias_right}\n')
        else:
           fout.write('    value = %10.5e\n'%0.0)
    elif u_right is None and use_moose_neumann_bc is False:
        fout.write('    type = NormalFluxBC\n')
        fout.write('    variable = u\n')
        fout.write('    boundary = right\n')
        if qn_bias_right is not None:
            fout.write('    bias = ${replace qn_bias_right}\n')
        if transfer_coeff_right is not None:
            fout.write('    transferCoeff = ${replace transfer_coeff_right}\n')
        if u_reference_right is not None:
            fout.write('    reference = ${replace u_reference_right}\n')
    else:
        fout.write('    type = DirichletBC\n')
        fout.write('    variable = u\n')
        fout.write('    boundary = right\n')
        fout.write('    value = ${replace u_right}\n')
    fout.write('  []\n')

    fout.write('  [south]\n')
    if u_bottom is None and use_moose_neumann_bc:
        fout.write('    type = NeumannBC\n')
        fout.write('    variable = u\n')
        fout.write('    boundary = bottom\n')
        if qn_bias_bottom is not None:
            fout.write('    value = ${replace qn_bias_bottom}\n')
        else:
           fout.write('    value = %10.5e\n'%0.0)
    elif u_bottom is None and use_moose_neumann_bc is False:
        fout.write('    type = NormalFluxBC\n')
        fout.write('    variable = u\n')
        fout.write('    boundary = bottom\n')
        if qn_bias_bottom is not None:
            fout.write('    bias = ${replace qn_bias_bottom}\n')
        if transfer_coeff_bottom is not None:
            fout.write('    transferCoeff = ${replace transfer_coeff_bottom}\n')
        if u_reference_bottom is not None:
            fout.write('    reference = ${replace u_reference_bottom}\n')
    else:
        fout.write('    type = DirichletBC\n')
        fout.write('    variable = u\n')
        fout.write('    boundary = bottom\n')
        fout.write('    value = ${replace u_bottom}\n')
    fout.write('  []\n')

    fout.write('  [north]\n')
    if u_top is None and use_moose_neumann_bc:
        fout.write('    type = NeumannBC\n')
        fout.write('    variable = u\n')
        fout.write('    boundary = top\n')
        if qn_bias_top is not None:
            fout.write('    value = ${replace qn_bias_top}\n')
        else:
           fout.write('    value = %10.5e\n'%0.0)
    elif u_top is None and use_moose_neumann_bc is False:
        fout.write('    type = NormalFluxBC\n')
        fout.write('    variable = u\n')
        fout.write('    boundary = top\n')
        if qn_bias_top is not None:
            fout.write('    bias = ${replace qn_bias_top}\n')
        if transfer_coeff_top is not None:
            fout.write('    transferCoeff = ${replace transfer_coeff_top}\n')
        if u_reference_top is not None:
            fout.write('    reference = ${replace u_reference_top}\n')
    else:
        fout.write('    type = DirichletBC\n')
        fout.write('    variable = u\n')
        fout.write('    boundary = top\n')
        fout.write('    value = ${replace u_top}\n')
    fout.write('  []\n')

    fout.write('[]\n')

    fout.write('\n')

    if solver=='fdp-newt-full':
        fout.write('[Preconditioning]\n')
        fout.write("  active = 'fdp-newt-full'\n")
        fout.write('  [fdp-newt-full]\n')
        fout.write('    type = FDP\n')
        fout.write('    full = true\n')
        fout.write("    solve_type = 'NEWTON'\n")
        fout.write("    petsc_options_iname = '-pc_type -mat_fd_coloring_err -mat_fd_type'\n")
        fout.write("    petsc_options_value = 'lu  %10.3e               ds'\n"%fdp_dh)
        fout.write('  []\n')
        fout.write('[]\n')

        fout.write('\n')

        fout.write('[Executioner]\n')
        fout.write('  type = Steady\n')
        fout.write('[]\n')
    elif solver=='pjfnk-hypre-boom':
        fout.write('[Executioner]\n')
        fout.write('  type = Steady\n')
        fout.write("  solve_type = 'PJFNK'\n")
        fout.write("  petsc_options_iname = '-pc_type -pc_hypre_type'\n")
        fout.write("  petsc_options_value = 'hypre boomeramg'\n")
        fout.write('[]\n')
    else:
        assert False,'Invalid solver option = %r'%solver

    if compute_energy:

        fout.write('\n')

        fout.write('[Postprocessors]\n')
        fout.write('  [bulk-energy]\n')
        fout.write('    type = BulkEnergy\n')
        fout.write("    execute_on = 'timestep_end final'\n")
        fout.write("    variable = 'u'     # bulk energy unknown variable\n")
        fout.write('    diffCoeff = ${replace diff_coeff}\n')
        fout.write('    sourceS = ${replace source_s}\n')
        fout.write('  []\n')
        if qn_bias_right is not None or transfer_coeff_right is not None or u_reference_right is not None:
            fout.write('  [boundary-energy]\n')
            fout.write('    type = BoundaryEnergy\n')
            fout.write("    execute_on = 'timestep_end final'\n")
            fout.write("    variable = 'u'     # bulk energy unknown variable\n")
            if qn_bias_right is not None:
                fout.write("    diffNormalFluxBias = ${replace qn_bias_right}\n")
            fout.write('    diffCoeff = ${replace diff_coeff}\n')
            if transfer_coeff_right is not None:
                fout.write('    transferCoeff = ${replace transfer_coeff_right}\n')
            if u_reference_right is not None:
                fout.write('    varReference = ${replace u_reference}\n')
            fout.write('    boundary = right\n')
            fout.write('  []\n')
            fout.write('  [total-energy]\n')
            fout.write('    type = LinearCombinationPostprocessor\n')
            fout.write("    pp_names = 'bulk-energy boundary-energy'\n")
            fout.write("    pp_coefs = '1  1'\n")
            fout.write('    b = 0\n')
            fout.write('  []\n')

        fout.write('[]\n')

    fout.write('\n')

    fout.write('[VectorPostprocessors]\n')
    fout.write('  [x-line]\n')
    fout.write('    type = LineValueSampler\n')
    fout.write("    execute_on = 'timestep_end final'\n")
    if compute_diffusion_flux:
        fout.write("    variable = 'u diffFluxU_x diffFluxU_y'    # output data\n")
    else:
        fout.write("    variable = 'u'    # output data\n")
    if x_line_ypos is not None:
        fout.write("    start_point = '${replace xmin} %10.5e 0'\n"%x_line_ypos)
        fout.write("    end_point = '${replace xmax}   %10.5e 0'\n"%x_line_ypos)
    else:
        fout.write("    start_point = '${replace xmin} ${fparse (ymax+ymin)/2} 0'\n")
        fout.write("    end_point = '${replace xmax} ${fparse (ymax+ymin)/2} 0'\n")
    fout.write('    num_points = %i\n'%n_plot_pts_x)
    fout.write('    sort_by = id\n')
    fout.write('  []\n')
    fout.write('  [y-line]\n')
    fout.write('    type = LineValueSampler\n')
    fout.write("    execute_on = 'timestep_end final'\n")
    if compute_diffusion_flux:
        fout.write("    variable = 'u diffFluxU_x diffFluxU_y'    # output data\n")
    else:
        fout.write("    variable = 'u'    # output data\n")
    if y_line_xpos is not None:
        fout.write("    start_point = '%10.5e ${replace ymin} 0'\n"%y_line_xpos)
        fout.write("    end_point = '%10.5e   ${replace ymax} 0'\n"%y_line_xpos)
    else:
        fout.write("    start_point = '${fparse (xmax+xmin)/2} ${replace ymin} 0'\n")
        fout.write("    end_point = '${fparse (xmax+xmin)/2} ${replace ymax} 0'\n")
    fout.write('    num_points = %i\n'%n_plot_pts_y)
    fout.write('    sort_by = id\n')
    fout.write('  []\n')
    fout.write('[]\n')

    fout.write('\n')

    fout.write('[Outputs]\n')
    fout.write('  console = true\n')
    fout.write('  [vtk]\n')
    fout.write('    type = VTK\n')
    fout.write('    execute_on = final\n')
    fout.write('    file_base = out\n')
    fout.write('  []\n')
    fout.write('  [x]\n')
    fout.write('    type = CSV\n')
    fout.write("    execute_on = 'final'\n")
    fout.write("    show = 'x-line'\n")
    fout.write('    file_base = out-x\n')
    fout.write('  []\n')
    fout.write('  [y]\n')
    fout.write('    type = CSV\n')
    fout.write("    execute_on = 'final'\n")
    fout.write("    show = 'y-line'\n")
    fout.write('    file_base = out-y\n')
    fout.write('  []\n')
    if compute_energy:
        fout.write('  [console-energy]\n')
        fout.write('    type = Console\n')
        fout.write("    execute_on = 'final linear nonlinear'\n")
        if qn_bias_right is None and transfer_coeff_right is None and u_reference_right is None:
            fout.write("    show = 'bulk-energy'\n")
        else:
            fout.write("    show = 'bulk-energy boundary-energy total-energy'\n")
        fout.write('  []\n')
        fout.write('  [file-energy]\n')
        fout.write('    type = CSV\n')
        fout.write("    execute_on = 'final'\n")
        fout.write("    file_base = 'out_energy'\n")
        if qn_bias_right is None and transfer_coeff_right is None and u_reference_right is None:
            fout.write("    show = 'bulk-energy'\n")
        else:
            fout.write("    show = 'bulk-energy boundary-energy total-energy'\n")
        fout.write('  []\n')
    fout.write('[]\n')

    fout.close()
#*********************************************************************************
def write_engy5310_p1_3d_input_file(
                                 x_west=None, x_east=None, y_south=None, y_north=None,
                                 z_bottom=None, z_top=None,
                                 u_west=None, u_east=None, u_south=None, u_north=None,
                                 u_bottom=None, u_top=None,
                                 qn_west_bias=None, transfer_coeff_west=None,
                                 u_reference_west=None,
                                 qn_east_bias=None, transfer_coeff_east=None,
                                 u_reference_east=None,
                                 qn_south_bias=None, transfer_coeff_south=None,
                                 u_reference_south=None,
                                 qn_north_bias=None, transfer_coeff_north=None,
                                 u_reference_north=None, qn_bias_bottom=None,
                                 transfer_coeff_bottom=None, u_reference_bottom=None,
                                 qn_bias_top=None, transfer_coeff_top=None,
                                 u_reference_top=None,
                                 diff_coeff=1.0, source_s=1.0,
                                 n_felem_x=1, n_felem_y=1, n_felem_z=1,
                                 order='second',
                                 n_plot_pts_x=10, n_plot_pts_y=10, n_plot_pts_z=10,
                                 compute_diffusion_flux=False,
                                 use_moose_diffusion_flux=False,
                                 compute_energy=False,
                                 use_moose_neumann_bc=False,
                                 solver='pjfnk-hypre-boom',
                                 fdp_dh=1e-6):

    if u_east is not None:
        assert qn_east_bias is None
        assert transfer_coeff_east is None
        assert u_reference_east is None
    if qn_east_bias is not None:
        assert u_east is None
    if transfer_coeff_east is not None:
        assert transfer_coeff_east >= 0.0
        assert u_east is None
    if u_reference_east is not None:
        assert u_east is None

    fout = open('engy5310p1/input.hit', 'wt')

    fout.write('# Engy-5310 Problem: Poisson 3D FEM\n')
    fout.write('# UMass Lowell Nuclear Chemical Engineering\n')
    fout.write('# Prof. Valmor F. de Almeida\n')
    time_date_stamp = datetime.datetime.today().strftime('%d%b%y %H:%M:%S')
    fout.write('# %s\n'%time_date_stamp)

    fout.write('\n')

    fout.write('# Parameters\n')
    fout.write('xmin = %10.5e\n'%x_west)
    fout.write('xmax = %10.5e\n'%x_east)
    fout.write('ymin = %10.5e\n'%y_south)
    fout.write('ymax = %10.5e\n'%y_north)
    fout.write('zmin = %10.5e\n'%z_bottom)
    fout.write('zmax = %10.5e\n'%z_top)
    fout.write('diff_coeff = %10.5e\n'%diff_coeff)
    fout.write('source_s = %10.5e\n'%source_s)

    if u_west is not None:
        fout.write('u_west = %10.5e\n'%u_west)
    if u_east is not None:
        fout.write('u_east = %10.5e\n'%u_east)
    if u_south is not None:
        fout.write('u_south = %10.5e\n'%u_south)
    if u_north is not None:
        fout.write('u_north = %10.5e\n'%u_north)
    if u_bottom is not None:
        fout.write('u_bottom = %10.5e\n'%u_bottom)
    if u_top is not None:
        fout.write('u_top = %10.5e\n'%u_top)

    if qn_west_bias:
        fout.write('qn_west_bias = %10.5e\n'%qn_west_bias)
    if transfer_coeff_west:
        fout.write('transfer_coeff_west = %10.5e\n'%transfer_coeff_west)
    if u_reference_west:
        fout.write('u_reference_west = %10.5e\n'%u_reference_west)

    if qn_east_bias:
        fout.write('qn_east_bias = %10.5e\n'%qn_east_bias)
    if transfer_coeff_east:
        fout.write('transfer_coeff_east = %10.5e\n'%transfer_coeff_east)
    if u_reference_east:
        fout.write('u_reference_east = %10.5e\n'%u_reference_east)

    if qn_south_bias:
        fout.write('qn_south_bias = %10.5e\n'%qn_south_bias)
    if transfer_coeff_south:
        fout.write('transfer_coeff_south_= %10.5e\n'%transfer_coeff_south)
    if u_reference_south:
        fout.write('u_reference_south = %10.5e\n'%u_reference_south)

    if qn_north_bias:
        fout.write('qn_north_bias = %10.5e\n'%qn_north_bias)
    if transfer_coeff_north:
        fout.write('transfer_coeff_north_= %10.5e\n'%transfer_coeff_north)
    if u_reference_north:
        fout.write('u_reference_north = %10.5e\n'%u_reference_north)

    if qn_bias_bottom:
        fout.write('qn_bias_bottom = %10.5e\n'%qn_bias_bottom)
    if transfer_coeff_bottom:
        fout.write('transfer_coeff_bottom_= %10.5e\n'%transfer_coeff_bottom)
    if u_reference_bottom:
        fout.write('u_reference_bottom = %10.5e\n'%u_reference_bottom)

    if qn_bias_top:
        fout.write('qn_bias_top = %10.5e\n'%qn_bias_top)
    if transfer_coeff_top:
        fout.write('transfer_coeff_top_= %10.5e\n'%transfer_coeff_top)
    if u_reference_top:
        fout.write('u_reference_top = %10.5e\n'%u_reference_top)

    fout.write('\n')

    fout.write('[Problem]\n')
    fout.write('  type = FEProblem\n')
    fout.write('  coord_type = XYZ\n')
    fout.write('[]\n')

    fout.write('\n')

    fout.write('[Mesh]\n')
    fout.write('  [3d]\n')
    fout.write('    type = GeneratedMeshGenerator\n')
    fout.write('    dim = 3\n')
    fout.write('    xmin = ${replace xmin}\n')
    fout.write('    xmax = ${replace xmax}\n')
    fout.write('    ymin = ${replace ymin}\n')
    fout.write('    ymax = ${replace ymax}\n')
    fout.write('    zmin = ${replace zmin}\n')
    fout.write('    zmax = ${replace zmax}\n')
    fout.write('    nx = %i\n'%n_felem_x)
    fout.write('    ny = %i\n'%n_felem_y)
    fout.write('    nz = %i\n'%n_felem_z)
    if order == 'second':
        fout.write('    elem_type = HEX27\n')
    fout.write('  []\n')
    fout.write('[]\n')

    fout.write('\n')

    fout.write('[Variables]\n')
    fout.write('  [u]\n')
    fout.write('    order = %s\n'%order)
    fout.write('    family = lagrange\n')
    values = list()
    values = [u_west, u_east, u_south, u_north, u_bottom, u_top]
    names = ['u_west', 'u_east', 'u_south', 'u_north', 'u_bottom', 'u_top']
    ids = [i for (i,v) in enumerate(values) if v is not None]
    if len(ids):
        average = str()
        for i in ids:
            average += names[i]+'+'
        average = '('+average[:-1]+')/'+str(len(ids))
        fout.write('    initial_condition = ${fparse '+average+'}\n')
    fout.write('  []\n')
    fout.write('[]\n')

    fout.write('\n')

    if compute_diffusion_flux:
        fout.write('[AuxVariables]\n')

        fout.write('  [diffFluxU]\n')
        if order in ['second', 'SECOND']:
            fout.write('    order = FIRST\n')
        elif order in ['first', 'FIRST']:
            fout.write('    order = CONSTANT\n')
        else:
            assert False
        fout.write('    family = MONOMIAL_VEC\n')
        fout.write('  []\n')

        fout.write('  [diffFluxU_x]\n')
        if order in ['second', 'SECOND']:
            fout.write('    order = FIRST\n')
        elif order in ['first', 'FIRST']:
            fout.write('    order = CONSTANT\n')
        else:
            assert False
        fout.write('    family = MONOMIAL\n')
        fout.write('  []\n')

        fout.write('  [diffFluxU_y]\n')
        if order in ['second', 'SECOND']:
            fout.write('    order = FIRST\n')
        elif order in ['first', 'FIRST']:
            fout.write('    order = CONSTANT\n')
        else:
            assert False
        fout.write('    family = MONOMIAL\n')
        fout.write('  []\n')

        fout.write('  [diffFluxU_z]\n')
        if order in ['second', 'SECOND']:
            fout.write('    order = FIRST\n')
        elif order in ['first', 'FIRST']:
            fout.write('    order = CONSTANT\n')
        else:
            assert False
        fout.write('    family = MONOMIAL\n')
        fout.write('  []\n')

        fout.write('[]\n')

        fout.write('\n')

    fout.write('[Kernels]\n')
    fout.write('  [diffusion-term]\n')
    fout.write('    type = DiffusionTerm\n')
    fout.write('    variable = u     # produced quantity\n')
    fout.write('    diffCoeff = ${replace diff_coeff}\n')
    fout.write('  []\n')
    fout.write('  [source-term]\n')
    fout.write('    type = SourceTerm\n')
    fout.write('    variable = u     \n')
    fout.write('    coupledVariable = u   \n')
    fout.write('    sourceS = ${replace source_s}\n')
    fout.write('  []\n')
    fout.write('[]\n')

    fout.write('\n')

    if compute_diffusion_flux:
        fout.write('[AuxKernels]\n')
        if use_moose_diffusion_flux:
            fout.write('  [diffusion-flux-x]\n')
            fout.write('    execute_on = timestep_end\n')
            fout.write('    type = DiffusionFluxAux     # provided by MOOSE\n')
            fout.write('    diffusion_variable = u\n')
            fout.write('    diffusivity = ${replace diff_coeff}\n')
            fout.write('    component = x\n')
            fout.write('    variable = diffFluxU_x     # produced quantity\n')
            fout.write('  []\n')
            fout.write('  [diffusion-flux-y]\n')
            fout.write('    execute_on = timestep_end\n')
            fout.write('    type = DiffusionFluxAux     # provided by MOOSE\n')
            fout.write('    diffusion_variable = u\n')
            fout.write('    diffusivity = ${replace diff_coeff}\n')
            fout.write('    component = y\n')
            fout.write('    variable = diffFluxU_y     # produced quantity\n')
            fout.write('  []\n')
            fout.write('  [diffusion-flux-z]\n')
            fout.write('    execute_on = timestep_end\n')
            fout.write('    type = DiffusionFluxAux     # provided by MOOSE\n')
            fout.write('    diffusion_variable = u\n')
            fout.write('    diffusivity = ${replace diff_coeff}\n')
            fout.write('    component = z\n')
            fout.write('    variable = diffFluxU_z     # produced quantity\n')
            fout.write('  []\n')
        else:
            fout.write('  [diffusion-flux]\n')
            fout.write('    execute_on = timestep_end\n')
            fout.write('    type = DiffusionFlux   # new kernel\n')
            fout.write('    field = u\n')
            fout.write('    diffCoeff = ${replace diff_coeff}\n')
            fout.write('    variable = diffFluxU     # produced quantity\n')
            fout.write('  []\n')
            fout.write('  [diffusion-flux-x]\n')
            fout.write('    execute_on = timestep_end\n')
            fout.write('    type = VectorVariableComponentAux  # provided by MOOSE\n')
            fout.write('    variable = diffFluxU_x    # produced quantity\n')
            fout.write('    component = x\n')
            fout.write('    vector_variable = diffFluxU   \n')
            fout.write('  []\n')
            fout.write('  [diffusion-flux-y]\n')
            fout.write('    execute_on = timestep_end\n')
            fout.write('    type = VectorVariableComponentAux  # provided by MOOSE\n')
            fout.write('    variable = diffFluxU_y    # produced quantity\n')
            fout.write('    component = y\n')
            fout.write('    vector_variable = diffFluxU   \n')
            fout.write('  []\n')
            fout.write('  [diffusion-flux-z]\n')
            fout.write('    execute_on = timestep_end\n')
            fout.write('    type = VectorVariableComponentAux  # provided by MOOSE\n')
            fout.write('    variable = diffFluxU_z    # produced quantity\n')
            fout.write('    component = z\n')
            fout.write('    vector_variable = diffFluxU   \n')
            fout.write('  []\n')
        fout.write('[]\n')
        fout.write('\n')

    fout.write('[BCs]\n')

    fout.write('  [west]\n')
    if u_west is None and use_moose_neumann_bc:
        fout.write('    type = NeumannBC\n')
        fout.write('    variable = u\n')
        fout.write('    boundary = left\n')
        if qn_west_bias:
            fout.write('    value = ${replace qn_west_bias}\n')
        else:
           fout.write('    value = %10.5e\n'%0.0)
    elif u_west is None and use_moose_neumann_bc is False:
        fout.write('    type = NormalFluxBC\n')
        fout.write('    variable = u\n')
        fout.write('    boundary = left\n')
        if qn_west_bias:
            fout.write('    bias = ${replace qn_west_bias}\n')
        if transfer_coeff_west:
            fout.write('    transferCoeff = ${replace transfer_coeff_west}\n')
        if u_reference_east:
            fout.write('    reference = ${replace u_reference_west}\n')
    else:
        fout.write('    type = DirichletBC\n')
        fout.write('    variable = u\n')
        fout.write('    boundary = left\n')
        fout.write('    value = ${replace u_west}\n')
    fout.write('  []\n')

    fout.write('  [east]\n')
    if u_east is None and use_moose_neumann_bc:
        fout.write('    type = NeumannBC\n')
        fout.write('    variable = u\n')
        fout.write('    boundary = right\n')
        if qn_east_bias:
            fout.write('    value = ${replace qn_east_bias}\n')
        else:
           fout.write('    value = %10.5e\n'%0.0)
    elif u_east is None and use_moose_neumann_bc is False:
        fout.write('    type = NormalFluxBC\n')
        fout.write('    variable = u\n')
        fout.write('    boundary = right\n')
        if qn_east_bias:
            fout.write('    bias = ${replace qn_east_bias}\n')
        if transfer_coeff_east:
            fout.write('    transferCoeff = ${replace transfer_coeff_east}\n')
        if u_reference_east:
            fout.write('    reference = ${replace u_reference_east}\n')
    else:
        fout.write('    type = DirichletBC\n')
        fout.write('    variable = u\n')
        fout.write('    boundary = right\n')
        fout.write('    value = ${replace u_east}\n')
    fout.write('  []\n')

    fout.write('  [south]\n')
    if u_south is None and use_moose_neumann_bc:
        fout.write('    type = NeumannBC\n')
        fout.write('    variable = u\n')
        fout.write('    boundary = bottom\n')
        if qn_south_bias:
            fout.write('    value = ${replace qn_south_bias}\n')
        else:
           fout.write('    value = %10.5e\n'%0.0)
    elif u_south is None and use_moose_neumann_bc is False:
        fout.write('    type = NormalFluxBC\n')
        fout.write('    variable = u\n')
        fout.write('    boundary = bottom\n')
        if qn_south_bias:
            fout.write('    bias = ${replace qn_south_bias}\n')
        if transfer_coeff_south:
            fout.write('    transferCoeff = ${replace transfer_coeff_south}\n')
        if u_reference_south:
            fout.write('    reference = ${replace u_reference_south}\n')
    else:
        fout.write('    type = DirichletBC\n')
        fout.write('    variable = u\n')
        fout.write('    boundary = bottom\n')
        fout.write('    value = ${replace u_south}\n')
    fout.write('  []\n')

    fout.write('  [north]\n')
    if u_north is None and use_moose_neumann_bc:
        fout.write('    type = NeumannBC\n')
        fout.write('    variable = u\n')
        fout.write('    boundary = top\n')
        if qn_north_bias:
            fout.write('    value = ${replace qn_north_bias}\n')
        else:
           fout.write('    value = %10.5e\n'%0.0)
    elif u_north is None and use_moose_neumann_bc is False:
        fout.write('    type = NormalFluxBC\n')
        fout.write('    variable = u\n')
        fout.write('    boundary = top\n')
        if qn_north_bias:
            fout.write('    bias = ${replace qn_north_bias}\n')
        if transfer_coeff_north:
            fout.write('    transferCoeff = ${replace transfer_coeff_north}\n')
        if u_reference_north:
            fout.write('    reference = ${replace u_reference_north\n')
    else:
        fout.write('    type = DirichletBC\n')
        fout.write('    variable = u\n')
        fout.write('    boundary = top\n')
        fout.write('    value = ${replace u_north}\n')
    fout.write('  []\n')

    fout.write('  [bottom]\n')
    if u_bottom is None and use_moose_neumann_bc:
        fout.write('    type = NeumannBC\n')
        fout.write('    variable = u\n')
        fout.write('    boundary = back\n')
        if qn_bias_bottom:
            fout.write('    value = ${replace qn_bias_bottom}\n')
        else:
           fout.write('    value = %10.5e\n'%0.0)
    elif u_bottom is None and use_moose_neumann_bc is False:
        fout.write('    type = NormalFluxBC\n')
        fout.write('    variable = u\n')
        fout.write('    boundary = back\n')
        if qn_bias_bottom:
            fout.write('    bias = ${replace qn_bias_bottom}\n')
        if transfer_coeff_bottom:
            fout.write('    transferCoeff = ${replace transfer_coeff_bottom}\n')
        if u_reference_bottom:
            fout.write('    reference = ${replace u_reference_bottom}\n')
    else:
        fout.write('    type = DirichletBC\n')
        fout.write('    variable = u\n')
        fout.write('    boundary = back\n')
        fout.write('    value = ${replace u_bottom}\n')
    fout.write('  []\n')

    fout.write('  [top]\n')
    if u_top is None and use_moose_neumann_bc:
        fout.write('    type = NeumannBC\n')
        fout.write('    variable = u\n')
        fout.write('    boundary = front\n')
        if qn_bias_top:
            fout.write('    value = ${replace qn_bias_top}\n')
        else:
           fout.write('    value = %10.5e\n'%0.0)
    elif u_top is None and use_moose_neumann_bc is False:
        fout.write('    type = NormalFluxBC\n')
        fout.write('    variable = u\n')
        fout.write('    boundary = front\n')
        if qn_bias_top:
            fout.write('    bias = ${replace qn_bias_top}\n')
        if transfer_coeff_top:
            fout.write('    transferCoeff = ${replace transfer_coeff_top}\n')
        if u_reference_top:
            fout.write('    reference = ${replace u_reference_top\n')
    else:
        fout.write('    type = DirichletBC\n')
        fout.write('    variable = u\n')
        fout.write('    boundary = front\n')
        fout.write('    value = ${replace u_top}\n')
    fout.write('  []\n')

    fout.write('[]\n')

    fout.write('\n')

    if solver=='fdp-newt-full':
        fout.write('[Preconditioning]\n')
        fout.write("  active = 'fdp-newt-full'\n")
        fout.write('  [fdp-newt-full]\n')
        fout.write('    type = FDP\n')
        fout.write('    full = true\n')
        fout.write("    solve_type = 'NEWTON'\n")
        fout.write("    petsc_options_iname = '-pc_type -mat_fd_coloring_err -mat_fd_type'\n")
        fout.write("    petsc_options_value = 'lu  %10.3e               ds'\n"%fdp_dh)
        fout.write('  []\n')
        fout.write('[]\n')

        fout.write('\n')

        fout.write('[Executioner]\n')
        fout.write('  type = Steady\n')
        fout.write('[]\n')
    elif solver=='pjfnk-hypre-boom':
        fout.write('[Executioner]\n')
        fout.write('  type = Steady\n')
        fout.write("  solve_type = 'PJFNK'\n")
        fout.write("  petsc_options_iname = '-pc_type -pc_hypre_type'\n")
        fout.write("  petsc_options_value = 'hypre boomeramg'\n")
        fout.write('[]\n')
    else:
        assert False,'Invalid solver option = %r'%solver

    if compute_energy:

        fout.write('\n')

        fout.write('[Postprocessors]\n')
        fout.write('  [bulk-energy]\n')
        fout.write('    type = BulkEnergy\n')
        fout.write("    execute_on = 'timestep_end final'\n")
        fout.write("    variable = 'u'     # bulk energy unknown variable\n")
        fout.write('    diffCoeff = ${replace diff_coeff}\n')
        fout.write('    sourceS = ${replace source_s}\n')
        fout.write('  []\n')
        if qn_east_bias is not None or transfer_coeff_east is not None or u_reference_east is not None:
            fout.write('  [boundary-energy]\n')
            fout.write('    type = BoundaryEnergy\n')
            fout.write("    execute_on = 'timestep_end final'\n")
            fout.write("    variable = 'u'     # bulk energy unknown variable\n")
            if qn_east_bias is not None:
                fout.write("    diffNormalFluxBias = ${replace qn_east_bias}\n")
            fout.write('    diffCoeff = ${replace diff_coeff}\n')
            if transfer_coeff_east is not None:
                fout.write('    transferCoeff = ${replace transfer_coeff_east}\n')
            if u_reference_east is not None:
                fout.write('    varReference = ${replace u_reference}\n')
            fout.write('    boundary = east\n')
            fout.write('  []\n')
            fout.write('  [total-energy]\n')
            fout.write('    type = LinearCombinationPostprocessor\n')
            fout.write("    pp_names = 'bulk-energy boundary-energy'\n")
            fout.write("    pp_coefs = '1  1'\n")
            fout.write('    b = 0\n')
            fout.write('  []\n')

        fout.write('[]\n')

    fout.write('\n')

    fout.write('[VectorPostprocessors]\n')
    fout.write('  [x-line]\n')
    fout.write('    type = LineValueSampler\n')
    fout.write("    execute_on = 'timestep_end final'\n")
    if compute_diffusion_flux:
        fout.write("    variable = 'u diffFluxU_x diffFluxU_y diffFluxU_z'    # output data\n")
    else:
        fout.write("    variable = 'u'    # output data\n")
    fout.write("    start_point = '${replace xmin} ${fparse (ymax+ymin)/2} ${fparse (zmax+zmin)/2}'\n")
    fout.write("    end_point = '${replace xmax} ${fparse (ymax+ymin)/2} ${fparse (zmax+zmin)/2}'\n")
    fout.write('    num_points = %i\n'%n_plot_pts_x)
    fout.write('    sort_by = id\n')
    fout.write('  []\n')

    fout.write('  [y-line]\n')
    fout.write('    type = LineValueSampler\n')
    fout.write("    execute_on = 'timestep_end final'\n")
    if compute_diffusion_flux:
        fout.write("    variable = 'u diffFluxU_x diffFluxU_y diffFluxU_z'    # output data\n")
    else:
        fout.write("    variable = 'u'    # output data\n")
    fout.write("    start_point = '${fparse (xmax+xmin)/2} ${replace ymin} ${fparse (zmax+zmin)/2}'\n")
    fout.write("    end_point = '${fparse (xmax+xmin)/2} ${replace ymax} ${fparse (zmax+zmin)/2}'\n")
    fout.write('    num_points = %i\n'%n_plot_pts_y)
    fout.write('    sort_by = id\n')
    fout.write('  []\n')

    fout.write('  [z-line]\n')
    fout.write('    type = LineValueSampler\n')
    fout.write("    execute_on = 'timestep_end final'\n")
    if compute_diffusion_flux:
        fout.write("    variable = 'u diffFluxU_x diffFluxU_y diffFluxU_z'  # output data\n")
    else:
        fout.write("    variable = 'u'    # output data\n")
    fout.write("    start_point = '${fparse (xmax+xmin)/2} ${fparse (ymax+ymin)/2} ${replace zmin}'\n")
    fout.write("    end_point = '${fparse (xmax+xmin)/2} ${fparse (ymax+ymin)/2} ${replace zmax}'\n")
    fout.write('    num_points = %i\n'%n_plot_pts_z)
    fout.write('    sort_by = id\n')
    fout.write('  []\n')

    fout.write('[]\n')

    fout.write('\n')

    fout.write('[Outputs]\n')
    fout.write('  console = true\n')
    fout.write('  [vtk]\n')
    fout.write('    type = VTK\n')
    fout.write('    execute_on = final\n')
    fout.write('    file_base = out\n')
    fout.write('  []\n')
    fout.write('  [x]\n')
    fout.write('    type = CSV\n')
    fout.write("    execute_on = 'final'\n")
    fout.write("    show = 'x-line'\n")
    fout.write('    file_base = out-x\n')
    fout.write('  []\n')
    fout.write('  [y]\n')
    fout.write('    type = CSV\n')
    fout.write("    execute_on = 'final'\n")
    fout.write("    show = 'y-line'\n")
    fout.write('    file_base = out-y\n')
    fout.write('  []\n')
    fout.write('  [z]\n')
    fout.write('    type = CSV\n')
    fout.write("    execute_on = 'final'\n")
    fout.write("    show = 'z-line'\n")
    fout.write('    file_base = out-z\n')
    fout.write('  []\n')
    if compute_energy:
        fout.write('  [console-energy]\n')
        fout.write('    type = Console\n')
        fout.write("    execute_on = 'final linear nonlinear'\n")
        if qn_east_bias is None and transfer_coeff_east is None and u_reference_east is None:
            fout.write("    show = 'bulk-energy'\n")
        else:
            fout.write("    show = 'bulk-energy boundary-energy total-energy'\n")
        fout.write('  []\n')
        fout.write('  [file-energy]\n')
        fout.write('    type = CSV\n')
        fout.write("    execute_on = 'final'\n")
        fout.write("    file_base = 'out_energy'\n")
        if qn_east_bias is None and transfer_coeff_east is None and u_reference_east is None:
            fout.write("    show = 'bulk-energy'\n")
        else:
            fout.write("    show = 'bulk-energy boundary-energy total-energy'\n")
        fout.write('  []\n')
    fout.write('[]\n')

    fout.close()
#*********************************************************************************
def run_engy5310_p1_moose_app():

    #cmd = 'cat /home/dealmeida/OneDrive/uml-courses/engy-5310/2021-01-05-spring/jupynb-repo/notebooks/engy5310p1/new_input.i'
    #cmd = 'ls /home/dealmeida/OneDrive/uml-courses/engy-5310/2021-01-05-spring/jupynb-repo/notebooks/engy5310p1'

    app = r'/home/dealmeida/OneDrive/uml-courses/engy-5310/2021-01-05-spring/jupynb-repo/notebooks/engy5310p1/engy5310p1-opt'
    input_file = r'/home/dealmeida/OneDrive/uml-courses/engy-5310/2021-01-05-spring/jupynb-repo/notebooks/engy5310p1/new_input.i'

    cmd = app+' -i '+input_file

    process = subprocess.run(cmd,
                             shell=True,
                             #check=True,
                             text=True,
                             capture_output=True)

    #print(process)
    #print(process.stdout)
    #print(process.stderr)
#*********************************************************************************
