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

from scipy.interpolate import interp1d
from scipy.integrate import quad
from matplotlib import pyplot as plt # import the pyplot function of the matplotlib package
#*********************************************************************************
def get_domain_partition(degree, n_elem, x_min, x_max, bc_x_min='', bc_x_max='essential'):
    #assert degree == 1
    # Local node numbering on parent domain
    # --0--------------1---->
    #  -1      0      +1    zeta
    gnodes_x = np.linspace(x_min, x_max, n_elem+1, dtype=np.float64)
    patches = list()
    local_to_global_node_id_map = []
    for e in range(n_elem):
        gnode_id_2 = e + 1            # right
        gnode_id_1 = gnode_id_2 - 1   # left
        x1 = gnodes_x[gnode_id_1]
        x2 = gnodes_x[gnode_id_2]
        # Local node id:  0   1
        patches.append((x1, x2))
        # Local node id:                        0           1
        local_to_global_node_id_map.append([gnode_id_1, gnode_id_2])
    if bc_x_min == 'essential':
        local_to_global_node_id_map[0][0] = -1
    if bc_x_max == 'essential':
        local_to_global_node_id_map[-1][-1] = -1
    return (patches, gnodes_x, local_to_global_node_id_map)
#*********************************************************************************
def get_parent_mapping():
    # zeta in [-1,1]
    parent_mapping = lambda zeta, x_e_bar, h_e: x_e_bar + h_e/2*zeta # compute x
    parent_mapping_prime = lambda h_e: h_e/2# compute mapping derivative wrt zeta
    # x in Omega_e
    inverse_parent_mapping = lambda x, x_e_bar, h_e: (x - x_e_bar)*2/h_e # compute zeta
    inverse_parent_mapping_prime = lambda h_e: 2/h_e # compute zeta
    #inverse_parent_mapping_der =
    return (parent_mapping, parent_mapping_prime, inverse_parent_mapping, inverse_parent_mapping_prime)
#*********************************************************************************
def get_parent_basis_functions():
    parent_basis_func_list = list()
    parent_basis_func_prime_list = []
    parent_basis_func_list.append(lambda zeta: -(zeta-1)/2)  # left
    parent_basis_func_list.append(lambda zeta:  (zeta+1)/2)  # right
    parent_basis_func_prime_list.append(lambda zeta: -1/2) # left
    parent_basis_func_prime_list.append(lambda zeta:  1/2) # right
    return (parent_basis_func_list, parent_basis_func_prime_list)
#*********************************************************************************
def global_basis_function(i, x, domain_partition, parent_mapping, parent_basis_functions):
    """Evaluate the ith global FE basis function and its derivative on x points.

    This is never needed in practice. It is here for demonstrating the theory.
    """
    try:
        len(x)
    except TypeError:
        x = np.array([x])

    if not isinstance(x, np.ndarray):
       assert isinstance(x, list) or isinstance(x, tuple)
       x = np.array(x)

    phi_i_x = np.copy(x) * 0.0 # initialization
    phi_prime_i_x = np.copy(x) * 0.0 # initialization

    patches = domain_partition[0]
    local_to_global_node_id_map = domain_partition[2]
    inverse_parent_mapping = parent_mapping_prime[2]
    inverse_parent_mapping_prime = parent_mapping[3]
    parent_basis_func_list = parent_basis_functions[0]
    parent_basis_func_prime_list = parent_basis_functions[1]

    # expensive reverse lookup
    for j, x_j in enumerate(x):
        for e, nodes_x in enumerate(patches):
            if nodes_x[0] <= x_j <= nodes_x[1]:
                n_lnodes = len(nodes_x)
                for I in range(n_lnodes):
                    if local_to_global_node_id_map[e][I] == i:
                        x_e_bar = (nodes_x[0] + nodes_x[1])/2
                        h_e = nodes_x[1] - nodes_x[0]
                        zetta = inverse_parent_mapping(x_j, x_e_bar, h_e)
                        d_zetta = inverse_parent_mapping_prime(h_e)

                        phi_i_x[j] = parent_basis_func_list[I](zetta)
                        phi_prime_i_x[j] = parent_basis_func_prime_list[I](d_zetta)
                break
    return [phi_i_x, phi_prime_i_x]
#*********************************************************************************
def get_global_basis_functions(domain_partition, parent_mapping, parent_basis_functions, global_basis_function):

    basis_func_list = list()
    basis_func_prime_list = []
    n_gnodes = domain_partition[1].size

    local_to_global_node_id_map = domain_partition[2]

    phi_i = lambda i, x: global_basis_function(i,x, domain_partition, parent_mapping, parent_basis_functions)[0]
    phi_prime_i = lambda i, x: global_basis_function(i,x, domain_partition, parent_mapping, parent_basis_functions)[1]

    visited = [False]*n_gnodes
    for e in range(n_elem):
        for I in range(len(local_to_global_node_id_map[e])):
            gnode_id = local_to_global_node_id_map[e][I]
            if gnode_id >= 0 and not visited[gnode_id]:
                      basis_func_list.append(lambda x, i=gnode_id: phi_i(i,x))
                      basis_func_prime_list.append(lambda x, i=gnode_id: phi_prime_i(i,x))
                      visited[gnode_id] = True

    assert len(basis_func_list) >= 1, 'There are no basis functions to build.'

    return [basis_func_list, basis_func_prime_list]
#*********************************************************************************
def FourierBasis(x,N,shift_fourier,Kappa):
    function=np.zeros((len(x),2*N+1))
    for i in range (len(x)):
        for j in range((2*N+1)):
            if (j==0):
                function[i,j]=1
            elif (j%2)==1:
                function[i,j]=np.cos(((j//2)+1)*Kappa*(x[i]-shift_fourier))
            elif (j%2)==0:
                function[i,j]=np.sin((j//2)*Kappa*(x[i]-shift_fourier))
    return np.array(function)
#*********************************************************************************
