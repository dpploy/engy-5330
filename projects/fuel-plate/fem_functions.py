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

    parent_mapping, parent_mapping_prime, inverse_parent_mapping, inverse_parent_mapping_prime = get_parent_mapping()
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
def get_global_basis_functions(domain_partition, parent_mapping, parent_basis_functions, global_basis_function,n_elem):

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









#*********************************************************************************
def convertToNumpy(inputData, twoOutputs=False):
    numpyArray = np.zeros((len(inputData[:]), len(inputData[0])))
    for tupCount, tup in enumerate(inputData[:]):
        for eleCount, ele in enumerate(tup):
            numpyArray[tupCount][eleCount] = ele
    if twoOutputs == True:
        return (numpyArray[:, 0], numpyArray[:, 1])
    else:
        return (numpyArray)
#*********************************************************************************
# TODO Rename this here and in `get_parent_basis_functions`
def _extracted_from_get_parent_basis_functions_10(parent_basis_func_list, parent_basis_func_prime_list):
    parent_basis_func_list.append(lambda zetta: -((zetta - 1) / 2) * (-zetta))
    parent_basis_func_list.append(lambda zetta: -(zetta - 1) * (zetta + 1))
    parent_basis_func_list.append(lambda zetta: (zetta + 1) / 2 * zetta)

    parent_basis_func_prime_list.append(lambda zetta: (2 * zetta - 1) / 2)
    parent_basis_func_prime_list.append(lambda zetta: -2 * zetta)
    parent_basis_func_prime_list.append(lambda zetta: (2 * zetta + 1) / 2)
#*********************************************************************************
def plot_func(domain_partition, phi_list, x_min,x_max,title='Lagrange Basis Functions'):
    import matplotlib.pyplot as plt

    plt.style.use('classic')
    plt.figure(1, figsize=(14, 5))

    npts = 500
    x_pts = np.linspace(x_min, x_max, npts)
    # x_pts = domain_partition[1]
    for (i, phi_i) in enumerate(phi_list):
        plt.plot(x_pts, phi_i(x_pts), '-', label=r'$\phi_{%i}$' % i)

    gnodes_x = domain_partition[1]
    plt.scatter(gnodes_x, np.zeros(gnodes_x.size), color='red', marker='x', s=80, label='nodes')

    plt.title(title, fontsize=20)
    plt.ylabel(r'$\phi_i(x)$', fontsize=18)
    plt.xlabel(r'$x$', fontsize=18)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.legend(loc='best', fontsize=12)
    plt.grid(True)
    plt.show()
#*********************************************************************************    
def global_basis_deriv(i, x, domain_partition, parent_mapping, parent_basis_functions):
    """Evaluate the ith global FE basis function and its derivative on x points.

    This is never needed in practice. It is here for demonstrating the theory.
    """

    try:
        len(x)
    except TypeError:
        x = np.array([x])

    if not isinstance(x, np.ndarray):
        assert isinstance(x, (list, tuple))
        x = np.array(x)

    phi_i_x = np.copy(x) * 0.0  # initialization
    phi_prime_i_x = np.copy(x) * 0.0  # initialization

    patches = domain_partition[0]
    local_to_global_node_id_map = domain_partition[2]
    inverse_parent_mapping = parent_mapping[2]
    parent_basis_func_list = parent_basis_functions[0]
    parent_basis_deriv_list = parent_basis_functions[1]
    n_elem = len(patches)
    # expensive reverse lookup
    for j, x_j in enumerate(x):
        for e, nodes_x in enumerate(patches):
            if nodes_x[0] <= x_j <= nodes_x[1]:
                # n_lnodes = len(nodes_x)
                n_lnodes = len(parent_basis_func_list)
                for I in range(n_lnodes):
                    if local_to_global_node_id_map[e][I] == i:
                        x_e_bar = (nodes_x[0] + nodes_x[1]) / 2
                        h_e = nodes_x[1] - nodes_x[0]
                        zetta = inverse_parent_mapping(x_j, x_e_bar, h_e)

                        phi_prime_i_x[j] = parent_basis_deriv_list[I](zetta)
                break
    return phi_prime_i_x
#*********************************************************************************
def assemble_gram_matrix(a_mtrx, domain_partition, parent_mapping, parent_basis_func_list):
    n = len(parent_basis_func_list[0])

    local_to_global_node_id_map = domain_partition[2]

    local_a_mtrx = np.zeros((n, n), dtype=np.float64)

    for I in range(n):
        parent_basis_func_I = parent_basis_func_list[0][I]
        for J in range(n):
            integrand = lambda zeta: parent_basis_func_I(zeta) * parent_basis_func_list[0][J](zeta)
            (local_a_mtrx[I, J], _) = scipy.integrate.quad(integrand, -1, 1)
    for (e, gnode_ids) in enumerate(local_to_global_node_id_map):
        (x_0, x_1) = domain_partition[0][e]
        h_e = x_1 - x_0
        parent_mapping_jacobian = parent_mapping[1](h_e)
        for (I, i) in enumerate(gnode_ids):
            for (J, j) in enumerate(gnode_ids):
                 a_mtrx[i, j] += local_a_mtrx[I, J] * parent_mapping_jacobian
#*********************************************************************************
def assemble_load_vector(b_vec, domain_partition, parent_mapping, parent_basis_func_list, f,constrain=False):
    n = len(parent_basis_func_list[0])
    local_to_global_node_id_map = domain_partition[2]

    local_b_vec = np.zeros(n, dtype=np.float64)

    for (e, gnode_ids) in enumerate(local_to_global_node_id_map):
        (x_0, x_1) = domain_partition[0][e]
        h_e = x_1 - x_0
        parent_mapping_jacobian = parent_mapping[1](h_e)
        x_e_bar = (x_0 + x_1) / 2.0

        for I in range(n):
            parent_basis_func_I = parent_basis_func_list[0][I]
            f_parent = lambda zeta, x_e_bar=x_e_bar, h_e=h_e: f(parent_mapping[0](zeta, x_e_bar, h_e))
            integrand = lambda zeta: f_parent(zeta) * parent_basis_func_I(zeta)
            (local_b_vec[I], _) = scipy.integrate.quad(integrand, -1, 1)


        for (I, i) in enumerate(gnode_ids):
            b_vec[i] += local_b_vec[I] * parent_mapping_jacobian
#*********************************************************************************
def genCollocationPts(x, y, x_tilde_pts):
    slopeList = []
    interceptList = []
    for i in range(1, len(x)):
        slope = (y[i] - y[i - 1]) / (x[i] - x[i - 1])
        intercept = y[i] - (slope * x[i])
        slopeList.append(slope)
        interceptList.append(intercept)
    y_tilde_pts = []
    for i in range(len(x_tilde_pts)):
        result = np.where(x <= x_tilde_pts[i])
        x_tildeLoc = result[-1][-1]
        if x_tildeLoc >= len(slopeList):
            x_tildeLoc = len(slopeList) - 1
        y_tilde = x_tilde_pts[i] * slopeList[x_tildeLoc] + interceptList[x_tildeLoc]
        y_tilde_pts.append(y_tilde)
    return (np.asarray(y_tilde_pts))
#*********************************************************************************

'''Functions for Gram Matrix'''

def index_hand_fourier(p):
    t0=0
    if p>=1:
        for _ in range(p):
            t0 += 2*N+1
    return t0

def FourierBasis_single(x_pts,j,K,shift_fourier):
    if (j==0):
        mat=1
    elif (j%2)==1:
        mat=np.cos(((j//2)+1)*K*(x_pts-shift_fourier))
    elif (j%2)==0:
        mat=np.sin((j//2)*K*(x_pts-shift_fourier))
    return mat

def derivative_FourierBasis_single(x_pts,j,K,shift_fourier):
    if (j==0):
        mat = 0
    elif (j%2)==1:
        mat = -(((j//2)+1)*K*np.sin(((j//2)+1)*K*(x_pts-shift_fourier)))
    elif (j%2)==0:
        mat = (((j//2)+1)*K*np.cos((j//2)*K*(x_pts-shift_fourier)))
    return mat

def FourierBasis_matrix(x_pts,N,K,shift_fourier):
    matrix = np.zeros((len(x_pts),2*N+1))
    for i in range (len(x_pts)):
        for j in range((2*N+1)):
            matrix = FourierBasis_single(x_pts[i],j,K,shift_fourier)
    return np.array(matrix)

def derivative_g_c(x,i,j,K,shift_fourier):
    return derivative_FourierBasis_single(
        x, i, K, shift_fourier
    ) * derivative_FourierBasis_single(x, j, K, shift_fourier)

def g_c(x,i,j,K,shift_fourier):
    return FourierBasis_single(x, i, K, shift_fourier) * FourierBasis_single(
        x, j, K, shift_fourier
    )

def f_das(x):
    dx = np.array(10**(-6))
    return ((u(x+dx)-u(x))/dx)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
def Fourierbasis_gram_matrix(N,K,x_min,x_max,shift_fourier):
    gram_matrx_con = np.zeros(((np.sum(2*N)+1),(np.sum(2*N)+1)))
    for p in range(1):
        for i in range(2*N+1):
            for j in range(2*N+1):
                gram_matrx_con[int(index_hand_fourier(p)+i),int(index_hand_fourier(p)+j)] = (quad(g_c,x_min,x_max,args = (i,j,K,shift_fourier), limit = 10000)[0])+g_c(x_min,i,j,K,shift_fourier)+g_c(x_max,i,j,K,shift_fourier)+derivative_g_c(x_min,i,j,K,shift_fourier)+derivative_g_c(x_max,i,j,K,shift_fourier)
    return gram_matrx_con

def Load_vec_con(x,i,K,shift_fourier):
    return (FourierBasis_single(x,i,K,shift_fourier)*u(x))

def derivative_Load_vec_con(x,i,K,shift_fourier):
    return (derivative_FourierBasis_single(x,i,K,shift_fourier)*u(x))

def Fourier_load_vector(N,K,x_min,x_max,shift_fourier):
    b_vec_c =np.zeros((np.sum(2*N)+1))
    for i in range(2*N + 1):
        b_vec_c[i] = (quad(Load_vec_con, x_min, x_max, args = (i,K,shift_fourier), limit = 10000)[0])+Load_vec_con(x_min,i,K,shift_fourier)+Load_vec_con(x_max,i,K,shift_fourier)+derivative_Load_vec_con(x_min,i,K,shift_fourier)+derivative_Load_vec_con(x_max,i,K,shift_fourier)
    return b_vec_c

def bestg_vec_func_con(x,N,c_star_vec_con,K,shift_fourier):
    a_mtrx = FourierBasis(x,N,shift_fourier,K)
    return a_mtrx@c_star_vec_con
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#



def FourierBasis_Mat(x,N,shift_fourier,Kappa):
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

def Derivative_FourierBasis_Mat(x,N,shift_fourier,Kappa):
    function = np.zeros((len(x),2*N+1))
    for i in range (len(x)):
        for j in range((2*N+1)):
            if (j==0):
                function[i,j]=0
            elif (j%2)==1:
                function[i,j]=-(((j//2)+1)*Kappa*np.sin(((j//2)*Kappa*(x[i]-shift_fourier))))
            elif (j%2)==0:
                function[i,j]=(((j//2)+1)*Kappa*np.cos((j//2)*Kappa*(x[i]-shift_fourier)))
    return np.array(function)

def Derivative2_FourierBasis_Mat(x,N,shift_fourier,Kappa):
    function = np.zeros((len(x),2*N+1))
    for i in range (len(x)):
        for j in range((2*N+1)):
            if (j==0):
                function[i,j]=0
            elif (j%2)==1:
                function[i,j]=-(((((j//2)+1)*Kappa)**2)*np.cos(((j//2)*Kappa*(x[i]-shift_fourier))))
            elif (j%2)==0:
                function[i,j]=-(((((j//2)+1)*Kappa)**2)*np.sin((j//2)*Kappa*(x[i]-shift_fourier)))
    return np.array(function)

def evaluation_matrix(x,N,shift_fourier,Kappa):
    return FourierBasis_Mat(x,N,shift_fourier,Kappa)

def evaluation_matrix_derivative1(x,N,shift_fourier,Kappa):
    return Derivative_FourierBasis_Mat(x,N,shift_fourier,Kappa)

def evaluation_matrix_derivative2(x,N,shift_fourier,Kappa):
    return Derivative2_FourierBasis_Mat(x,N,shift_fourier,Kappa)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
