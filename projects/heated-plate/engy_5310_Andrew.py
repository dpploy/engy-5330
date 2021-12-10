import numpy as np
from matplotlib import pyplot as plt
import scipy
import scipy.integrate
import scipy.interpolate

def convertToNumpy(inputData, twoOutputs=False):

    numpyArray = np.zeros((len(inputData[:]), len(inputData[0])))
    tupCount = 0
    for tup in inputData[:]:
        eleCount = 0
        for ele in tup:
            numpyArray[tupCount][eleCount] = ele
            eleCount += 1
        tupCount += 1
    if twoOutputs == True:
        return (numpyArray[:, 0], numpyArray[:, 1])
    else:
        return (numpyArray)

def genFuncDerivAtPoint(x,y, xVal):
    slopeList = []
    interceptList = []
    for i in range(1, len(x)):
        slope = (y[i] - y[i - 1]) / (x[i] - x[i - 1])
        intercept = y[i] - (slope * x[i])
        slopeList.append(slope)
        interceptList.append(intercept)



    result = np.where(x <= xVal)
    x_tildeLoc = result[-1][-1]

    if x_tildeLoc >= len(slopeList):
        x_tildeLoc = len(slopeList) - 1
    slope = slopeList[x_tildeLoc]


    return(slope)

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

def genCollocationPtsDeriv(x, y, x_tilde_pts):
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
        y_tilde = slopeList[x_tildeLoc]
        y_tilde_pts.append(y_tilde)
    return (np.asarray(y_tilde_pts))


def genCollocationPts_Step(x, y, x_colloc_pts):
    xpts = [x[0]]
    ypts = [y[0]]
    for i in range(1, len(x)):
        xloc = (x[i] - x[i - 1]) / 2 + x[i - 1]
        xpts.append(xloc)
        xpts.append(xloc)
        ypts.append(y[i - 1])
        ypts.append(y[i])

        xpts.append(x[i])
        ypts.append(y[i])
    x = np.asarray(xpts)
    y = np.asarray(ypts)
    y_colloc_pts = []
    for i in range(len(x_colloc_pts)):
        result = np.where(x <= x_colloc_pts[i])
        x_collocLoc = result[-1][-1]
        y_colloc_pts.append(y[x_collocLoc])
    return (np.asarray(y_colloc_pts))


def genStepCollocPts(x, y, x_colloc_pts):
    xpts = [x[0]]
    ypts = [y[0]]
    for i in range(1, len(x)):
        xloc = (x[i] - x[i - 1]) / 2 + x[i - 1]
        xpts.append(xloc)

    x = np.asarray(xpts)
    y = np.asarray(ypts)
    y_colloc_pts = []
    for i in range(len(x_colloc_pts)):
        result = np.where(x < x_colloc_pts[i])
        #print(result)
        x_collocLoc = result[-1][-1]
        y_colloc_pts.append(y[x_collocLoc])
    return (np.asarray(y_colloc_pts))


def build_fourier_linear_system(mu, n_pairs, time_vec, no_mass_cc_vec):
    a_mtrx = np.zeros((len(time_vec), 2 * n_pairs + 1))
    b_vec = np.copy(no_mass_cc_vec)

    for j in range(len(time_vec)):
        for k in range(n_pairs):
            a_mtrx[:, 0] = 1
            a_mtrx[j, (k * 2) + 1] = np.cos((k + 1) * mu * time_vec[j])
            a_mtrx[j, (k + 1) * 2] = np.sin((k + 1) * mu * time_vec[j])

    return (a_mtrx, b_vec)


def genStepFunction(x, y, offset=0.0):
    xpts = [x[0]]
    ypts = [y[0]]
    for i in range(1, len(x)):
        xloc = (x[i] - x[i - 1]) / 2 + x[i - 1]
        xpts.append(xloc)
        xpts.append(xloc + offset)
        ypts.append(y[i - 1])
        ypts.append(y[i])

        xpts.append(x[i])
        ypts.append(y[i])
    return (np.asarray(xpts), np.asarray(ypts))


def gen_A_mtrx(kappa, m, n, colloc_x_vec,centered=False):
    if centered:
        x_bar = (colloc_x_vec[-1]-colloc_x_vec[0])/2
    else:
        x_bar=0
    A_mtrx = np.zeros((m, n))
    A_mtrx[:, 0] = 1
    for colNum in range(1, int((n - 1) / 2) + 1):
        A_mtrx[:, 2 * colNum - 1] = np.cos((colloc_x_vec[:]-x_bar) * kappa * (colNum))
        A_mtrx[:, 2 * colNum] = np.sin((colloc_x_vec-x_bar) * kappa * (colNum))
    return (A_mtrx)


def fourierSolution(x, kappa, n, c_vec,centered=False):
    solMtrx = gen_A_mtrx(kappa, m=1, n=n, colloc_x_vec=np.array([x]),centered=centered)
    solValue = solMtrx @ c_vec
    return (solValue[0])


def wavSolution(x, NList, bunch_pts, shift, sigmas, kappa, c_vec, centered=False):
    solMtrx = genWavMtrx_Updated(NList, bunch_pts, shift, sigmas, kappa, x_colloc_pts=np.array([x]),centered=centered)
    solValue = solMtrx @ c_vec
    return (solValue[0])


def wavSolutionSquared(x, NList, bunch_pts, shift, sigmas, kappa, c_vec,centered=False):
    result = wavSolution(x, NList, bunch_pts, shift, sigmas, kappa, c_vec,centered=centered)
    return (result ** 2)


def gen_wav_residual(x, NList, bunch_pts, shift, sigmas, kappa, c_vec, xList, yList):
    f = genCollocationPts(xList, yList, np.array([x]))[0]
    g = wavSolution(x, NList, bunch_pts, shift, sigmas, kappa, c_vec)
    return (f - g)


def gen_f_pt_wav_squared(x, xList, yList):
    result = genCollocationPts(xList, yList, np.array([x]))[0]
    return (result ** 2)


def fourierSolutionSquared(x, kappa, n, c_vec,centered=False):
    solMtrx = gen_A_mtrx(kappa, m=1, n=n, colloc_x_vec=np.array([x]),centered=centered)
    solValue = solMtrx @ c_vec
    return (solValue[0] ** 2)


def genStepChangeList(x, y):
    xpts = [x[0]]
    ypts = [y[0]]
    for i in range(1, len(x)):
        xloc = (x[i] - x[i - 1]) / 2 + x[i - 1]
        xpts.append(xloc)
        xpts.append(xloc)
        ypts.append(y[i - 1])
        ypts.append(y[i])

        xpts.append(x[i])
        ypts.append(y[i])
    x = np.asarray(xpts)
    y = np.asarray(ypts)
    return (x, y)


def gen_f_pt_squared(x, xStepChangeList, yStepChangeList):
    result = np.where(x < xStepChangeList)
    return ((yStepChangeList[result[-1][0]] ** 2))


def genWavMtrx(N, bunch_pts, shift, sigma, kappa, x_colloc_pts):
    # Old function which sums the bounding gaussians for each packet
    # Can only accept 1 sigma, # modes for all packets
    P = len(bunch_pts)
    m = len(x_colloc_pts)
    wav_mtrx = np.zeros((m, int(2 * N * P) + 1))
    k = 0
    x = x_colloc_pts
    for p in bunch_pts:
        wav_mtrx[:, 0] += (1 / (np.sqrt(2 * np.pi) * sigma)) * np.exp(-((x - p) ** 2) / (2 * (sigma ** 2)))

    for p_loc in range(1, len(bunch_pts) + 1):
        p = bunch_pts[p_loc - 1]
        for k in range(1, N + 1):
            loc = (p_loc - 1) * N * 2 + 2 * (k - 1) + 1
            wav_mtrx[:, loc] = (1 / (np.sqrt(2 * np.pi) * sigma)) * np.exp(-((x - p) ** 2) / (2 * (sigma ** 2))) * (
                np.cos(k * kappa * x + (shift * x ** 2) / 2))
            wav_mtrx[:, loc + 1] = (1 / (np.sqrt(2 * np.pi) * sigma)) * np.exp(-((x - p) ** 2) / (2 * (sigma ** 2))) * (
                np.sin(k * kappa * x + (shift * x ** 2) / 2))
    wav_mtrx = wav_mtrx * np.sqrt(2 * np.pi) * sigma
    return (wav_mtrx)


def genEquallySpacedBunchPts(x_min, x_max, n):
    # Generate equally spaced means for wavelet packets
    wavelength = x_max - x_min
    mean_dist = wavelength / (n + 1)
    bunch_pts = []
    for point_num in range(1, n + 1):
        bunch_pt = x_min + point_num * mean_dist
        bunch_pts.append(bunch_pt)
    return (bunch_pts)


def genWavMtrx_Updated(NList, bunch_pts, shift, sigmas, kappa, x_colloc_pts,centered=False):
    # Initialize wavelet matrix
    P = len(bunch_pts)
    m = len(x_colloc_pts)
    num_eqns = int(sum(NList) * 2 + len(NList))
    wav_mtrx = np.zeros((m, num_eqns))
    k = 0
    x = x_colloc_pts
    for N_index, N in enumerate(NList):
        # Iterates through each packet

        # Location of the column where the packet starts
        baseColNum = int(sum(NList[:N_index])) * 2 + N_index

        # Initialize sigma and mean for gaussian
        sigma = sigmas[N_index]
        p = bunch_pts[N_index]

        if centered:
            x_bar = p
        else:
            x_bar = 0

        # Create bounding gaussian for the target packet
        wav_mtrx[:, baseColNum] = np.exp(-((x - p) ** 2) / (2 * (sigma ** 2)))
        for modeNum in range(N):
            # Iterates through the rest of the packet

            # loc is col num where mode begins
            loc = baseColNum + (modeNum) * 2 + 1

            k = modeNum + 1

            # Fill sin and cos funcs
            wav_mtrx[:, loc] = np.exp(-((x - p) ** 2) / (2 * (sigma ** 2))) * (
                np.cos(k * kappa * (x-x_bar) + (shift * x ** 2) / 2))
            wav_mtrx[:, loc + 1] = np.exp(-((x - p) ** 2) / (2 * (sigma ** 2))) * (
                np.sin(k * kappa * (x-x_bar) + (shift * x ** 2) / 2))

    return (wav_mtrx)


def plotWavelet(NList, bunch_pts, shift, sigmas, kappa, x_min, x_max, x, y, num_plotting_pts=1000,
                continuous=False,centered=False):

    x_colloc_pts = np.linspace(x_min, x_max, num_plotting_pts)

    if continuous:
        phiVec = wavBasisFunctionVec(bunch_pts, NList, shift, sigmas, kappa, x_min, x_max, x, y,centered).func

        fig = plt.figure(figsize=(18, 5))
        ax = plt.subplot(111)
        plt.grid()
        plt.title(
            f"Wavelet Basis Function # Modes = {NList}, Packets = {bunch_pts} | Plotted from basis function vector")
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

        for N_index, N in enumerate(NList):
            baseColNum = int(sum(NList[:N_index])) * 2 + N_index
            plt.plot(x_colloc_pts, phiVec[baseColNum](x_colloc_pts), "k-", label="exp(.)")
            for modeNum in range(N):
                # Iterates through an individual packet
                loc = baseColNum + (modeNum) * 2 + 1

                k = modeNum + 1
                plt.plot(x_colloc_pts, phiVec[loc](x_colloc_pts), label=r"exp(.)cos(%i $\kappa x + \phi$ )" % k)
                plt.plot(x_colloc_pts, phiVec[loc+1](x_colloc_pts), label=r"exp(.)sin(%i $\kappa x + \phi$ )" % k)
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        plt.show()

    else:
        # num_plotting_pts = 1000
        wav_mtrx = genWavMtrx_Updated(NList, bunch_pts, shift, sigmas, kappa, x_colloc_pts)

        fig = plt.figure(figsize=(18, 5))
        ax = plt.subplot(111)

        plt.grid()
        plt.title(f"Wavelet Basis Function # Modes = {NList}, Packets = {bunch_pts}")
        for N_index, N in enumerate(NList):
            baseColNum = int(sum(NList[:N_index])) * 2 + N_index
            plt.plot(x_colloc_pts, wav_mtrx[:, baseColNum], "k-", label="exp(.)")
            for modeNum in range(N):
                # Iterates through an individual packet
                loc = baseColNum + (modeNum) * 2 + 1

                k = modeNum + 1

                plt.plot(x_colloc_pts, wav_mtrx[:, loc], label=r"exp(.)cos(%i $\kappa x + \phi$ )" % k)
                plt.plot(x_colloc_pts, wav_mtrx[:, loc + 1], label=r"exp(.)sin(%i $\kappa x + \phi$ )" % k)

        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        # plt.legend()
        plt.show()


def genExpanded_g_vec(NList, bunch_pts, shift, sigmas, kappa, c_vec, x_min, x_max, num_pts=10000,centered=False):
    x_colloc_pts = np.linspace(x_min, x_max, num_pts)

    wav_mtrx = genWavMtrx_Updated(NList, bunch_pts, shift, sigmas, kappa, x_colloc_pts,centered)
    g_vec = wav_mtrx @ c_vec
    return (g_vec)

def genExpanded_g_vec_fourier(n, kappa, c_vec, x_min, x_max, num_pts=10000,centered=False):
    x_colloc_pts = np.linspace(x_min, x_max, num_pts)

    A_mtrx = gen_A_mtrx(kappa,num_pts,int(2*n+1),x_colloc_pts,centered=centered)
    g_vec = A_mtrx @ c_vec
    return (g_vec)


def genLoadMtrx(NList, bunch_pts, shift, sigmas, kappa, x_min, x_max, xList, yList, constrain=False, centered=False):
    # Initialize load matrix
    num_eqns = int(sum(NList) * 2 + len(NList))
    loadMtrx = np.zeros(num_eqns)

    # Creates vector with all the basis functions [phi_0, phi_1, ..., phi_n]
    phiVec = wavBasisFunctionVec(bunch_pts, NList, shift, sigmas, kappa, x_min, x_max, xList, yList, centered)

    # Iterates through each packet
    for N_index, N in enumerate(NList):

        # Location of the column where the packet starts
        baseColNum = int(sum(NList[:N_index])) * 2 + N_index

        # Create bounding gaussian for the target packet
        phi_i = phiVec.func[baseColNum]
        f = lambda x: genCollocationPts(xList, yList, np.asarray([x]))[0]

        targetFunc = lambda x: f(x) * phi_i(x)

        if constrain:
            loadMtrx[baseColNum] = scipy.integrate.quad(targetFunc, x_min, x_max, limit=1000,
                            epsrel=10e-5, epsabs=10e-5)[0] \
                                   + phi_i(x_min) * f(x_min) + phi_i(x_max) * f(x_max)
        else:
            loadMtrx[baseColNum] = \
            scipy.integrate.quad(targetFunc, x_min, x_max, limit=1000,
                            epsrel=10e-5, epsabs=10e-5)[0]

        # Iterates through the rest of the packet
        for modeNum in range(1,N+1):

            # loc is col num where mode begins
            loc = baseColNum + (modeNum) * 2 - 1

            # Fill Cosine Func
            phi_i = phiVec.func[loc]
            targetFunc = lambda x: f(x) * phi_i(x)
            if constrain:
                loadMtrx[loc] = scipy.integrate.quad(targetFunc, x_min, x_max, limit=1000,
                            epsrel=10e-5, epsabs=10e-5)[0]\
                                + phi_i(x_min) * f(x_min) + phi_i(x_max) * f(x_max)
            else:
                loadMtrx[loc] = scipy.integrate.quad(targetFunc, x_min, x_max, limit=1000,
                            epsrel=10e-5, epsabs=10e-5)[0]


            # Fill Sine Func
            phi_i = phiVec.func[loc + 1]
            targetFunc = lambda x: f(x) * phi_i(x)
            if constrain:
                loadMtrx[loc + 1] = scipy.integrate.quad(targetFunc, x_min, x_max, limit=1000,
                            epsrel=10e-5, epsabs=10e-5)[0]\
                                    + phi_i(x_min) * f(x_min) + phi_i(x_max) * f(x_max)
            else:
                loadMtrx[loc + 1] = scipy.integrate.quad(targetFunc, x_min, x_max, limit=1000,
                            epsrel=10e-5, epsabs=10e-5)[0]

    return (loadMtrx)

"""
def genGramMtrx(NList, bunch_pts, shift, sigmas, kappa, x_min, x_max, xList, yList, constrain=False, centered=False):
    '''
    Creates an nxn matrix of zeros then populates it. It begins iterating through each column,
    filling in all rows above and including the diagonal. Because the Gram matrix is symmetric,
    it then copies all entries across the diagonal.
    '''

    # Creates vector with all the basis functions [phi_0, phi_1, ..., phi_n]
    phiVec = wavBasisFunctionVec(bunch_pts, NList, shift, sigmas, kappa, x_min, x_max, xList, yList, centered)

    # Initilize parameters relating to mode
    num_eqns = int(sum(NList) * 2 + len(NList))

    # Generate nxn Gram matrix
    gramMtrx = np.zeros((num_eqns, num_eqns))

    for N_index, N in enumerate(NList):
        # Iterates through each packet

        # Location of the column where the packet starts
        baseColNum = int(sum(NList[:N_index])) * 2 + N_index

        # Create bounding gaussian for the target packet
        for distToDiag in range(baseColNum + 1):
            phi_j = phiVec.func[baseColNum]
            phi_i = phiVec.func[distToDiag]
            targetFunc = lambda x: phi_j(x) * phi_i(x)

            if constrain:
                gramMtrx[distToDiag, baseColNum] = scipy.integrate.quad(targetFunc, x_min, x_max, limit=1000,
                            epsrel=10e-5, epsabs=10e-5)[0] + phi_i(x_min) * phi_j(x_min) + phi_i(x_max) * phi_j(x_max)
            else:
                gramMtrx[distToDiag, baseColNum] = scipy.integrate.quad(targetFunc, x_min, x_max, limit=1000,
                            epsrel=10e-5, epsabs=10e-5)[0]


        # Iterates through the rest of the packet
        for modeNum in range(N):
            # loc is col num where mode begins
            loc = baseColNum + (modeNum) * 2 + 1

            # Generates all exp(.)cos(.) * phi_i functions down the j'th column
            phi_j = phiVec.func[loc]
            for distToDiag in range(loc + 1):
                phi_i = phiVec.func[distToDiag]
                targetFunc = lambda x: phi_i(x) * phi_j(x)

                if constrain:
                    gramMtrx[distToDiag, loc] = scipy.integrate.quad(targetFunc, x_min, x_max, limit=1000,
                        epsrel=10e-5, epsabs=10e-5)[0] + phi_i(x_min) * phi_j(x_min) + phi_i(x_max) * phi_j(x_max)
                else:
                    gramMtrx[distToDiag, loc] = scipy.integrate.quad(targetFunc, x_min, x_max, limit=1000, epsrel=10e-5,
                        epsabs=10e-5)[0]


            # Generate exp(.)sin(.) * phi_i functions down the j+1'th column
            phi_j = phiVec.func[loc + 1]
            for distToDiag in range(loc + 2):
                phi_i = phiVec.func[distToDiag]
                targetFunc = lambda x: phi_i(x) * phi_j(x)

                if constrain:
                    gramMtrx[distToDiag, loc + 1] = scipy.integrate.quad(targetFunc, x_min, x_max, limit=1000,
                        epsrel=10e-5, epsabs=10e-5)[0] + phi_i(x_min) * phi_j(x_min) + phi_i(x_max) * phi_j(x_max)
                else:
                    gramMtrx[distToDiag, loc + 1] = scipy.integrate.quad(targetFunc, x_min, x_max, limit=1000,
                         epsrel=10e-5, epsabs=10e-5)[0]


        # Copies elements of the gram matrix across the diagonal
        for rowNum in range(num_eqns):
            for i in range(rowNum):
                gramMtrx[rowNum, i] = gramMtrx[i, rowNum]

    return (gramMtrx)
"""
def genGramMtrx(NList, bunch_pts, shift, sigmas, kappa, x_min, x_max, xList, yList, constrain=False, centered=False):
    '''
    Creates an nxn matrix of zeros then populates it. It begins iterating through each column,
    filling in all rows above and including the diagonal. Because the Gram matrix is symmetric,
    it then copies all entries across the diagonal.
    '''

    # Creates vector with all the basis functions [phi_0, phi_1, ..., phi_n]
    phiVec = wavBasisFunctionVec(bunch_pts, NList, shift, sigmas, kappa, x_min, x_max, xList, yList, centered)

    # Initilize parameters relating to mode
    num_eqns = int(sum(NList) * 2 + len(NList))

    # Generate nxn Gram matrix
    gramMtrx = np.zeros((num_eqns, num_eqns))

    for N_index, N in enumerate(NList):
        # Iterates through each packet

        # Location of the column where the packet starts
        baseColNum = int(sum(NList[:N_index])) * 2 + N_index

        # Create bounding gaussian for the target packet
        for distToDiag in range(baseColNum + 1):
            phi_j = phiVec.func[baseColNum]
            phi_i = phiVec.func[distToDiag]
            targetFunc = lambda x: phi_j(x) * phi_i(x)

            if constrain:
                gramMtrx[distToDiag, baseColNum] = scipy.integrate.quad(targetFunc, x_min, x_max, limit=100)[0] + phi_i(x_min) * phi_j(x_min) + phi_i(x_max) * phi_j(x_max)
            else:
                gramMtrx[distToDiag, baseColNum] = scipy.integrate.quad(targetFunc, x_min, x_max, limit=100)[0]


        # Iterates through the rest of the packet
        for modeNum in range(1,N+1):
            # loc is col num where mode begins
            loc = baseColNum + (modeNum) * 2 - 1

            # Generates all exp(.)cos(.) * phi_i functions down the j'th column
            phi_j = phiVec.func[loc]
            for distToDiag in range(loc + 1):
                phi_i = phiVec.func[distToDiag]
                targetFunc = lambda x: phi_i(x) * phi_j(x)

                if constrain:
                    gramMtrx[distToDiag, loc] = scipy.integrate.quad(targetFunc, x_min, x_max, limit=100)[0] + phi_i(x_min) * phi_j(x_min) + phi_i(x_max) * phi_j(x_max)
                else:
                    gramMtrx[distToDiag, loc] = scipy.integrate.quad(targetFunc, x_min, x_max, limit=100)[0]


            # Generate exp(.)sin(.) * phi_i functions down the j+1'th column
            phi_j = phiVec.func[loc + 1]
            for distToDiag in range(loc + 2):
                phi_i = phiVec.func[distToDiag]
                targetFunc = lambda x: phi_i(x) * phi_j(x)

                if constrain:
                    gramMtrx[distToDiag, loc + 1] = scipy.integrate.quad(targetFunc, x_min, x_max, limit=100)[0] + phi_i(x_min) * phi_j(x_min) + phi_i(x_max) * phi_j(x_max)
                else:
                    gramMtrx[distToDiag, loc + 1] = scipy.integrate.quad(targetFunc, x_min, x_max, limit=100)[0]


        # Copies elements of the gram matrix across the diagonal
        for rowNum in range(num_eqns):
            for i in range(rowNum):
                gramMtrx[rowNum, i] = gramMtrx[i, rowNum]

    return (gramMtrx)


class wavBasisFunctionVec:
    """
    Creates a vector of wavelet basis functions that can be called as a list with the .func method.
    Really could just be done with just a list, but I had already built this before I realized that...
    """

    def __init__(self, bunch_pts, NList, shift, sigmas, kappa, x_min, x_max, xList, yList, centered=False):
        self.func = []
        self.initBasisFuncs(bunch_pts, NList, shift, sigmas, kappa, centered)

    def initBasisFuncs(self, bunch_pts, NList, shift, sigmas, kappa, centered):
        # Iterates through each packet
        for N_index, N in enumerate(NList):
            # Location of the column where the packet starts
            baseColNum = int(sum(NList[:N_index])) * 2 + N_index

            # Initialize sigma and mean for gaussian
            sigma = sigmas[N_index]
            p = bunch_pts[N_index]

            if centered:
                x_bar = p
            else:
                x_bar = 0
            # Fill bounding gaussian function
            def targetFunc(x, p=p, sigma=sigma):
                return (np.exp(-((x - p) ** 2) / (2 * (sigma ** 2))))

            self.func.append(targetFunc)

            for modeNum in range(N):
                # Iterates through the rest of the packet

                # Calculates k value
                k = modeNum + 1

                # Fill Cosine Func
                def targetFunc(x, p=p, sigma=sigma, k=k, kappa=kappa, shift=shift,x_bar=x_bar):
                    return (np.exp(-((x - p) ** 2) / (2 * (sigma ** 2))) * (
                        np.cos(k * kappa * (x-x_bar) + (shift * x ** 2) / 2)))

                self.func.append(targetFunc)

                # Fill Sine Func
                def targetFunc(x, p=p, sigma=sigma, k=k, kappa=kappa, shift=shift,x_bar=x_bar):
                    return (np.exp(-((x - p) ** 2) / (2 * (sigma ** 2))) * (
                        np.sin(k * kappa * (x-x_bar) + (shift * x ** 2) / 2)))

                self.func.append(targetFunc)

class fourierBasisFunctionVec:
    """
    Creates a vector of fourier basis functions that can be called as a list with the .func method.
    Really could just be done with just a list, but I had already built this before I realized that...
    """

    def __init__(self, n, kappa, x_min,x_max, centered=False):
        self.func = []
        self.deriv = []
        self.initBasisFuncs(n, kappa, x_min, x_max, centered)
        self.initDerivFuncs(n, kappa, x_min, x_max, centered)

    def initBasisFuncs(self, n, kappa, x_min,x_max, centered):
        if centered:
            x_bar = (x_max + x_min) / 2
        else:
            x_bar = 0

        def targetFunc(x):
            return(x*0 + 1)
        self.func.append(targetFunc)

        for colNum in range(1, n+1):
            k = colNum

            def cosFunc(x, x_bar=x_bar, kappa=kappa,k=k):
                return (np.cos((x-x_bar) * kappa * k))
            self.func.append(cosFunc)

            def sinFunc(x, x_bar=x_bar, kappa=kappa,k=k):
                return (np.sin((x-x_bar) * kappa * k))
            self.func.append(sinFunc)

    def initDerivFuncs(self,n,kappa,x_min,x_max,centered):
        if centered:
            x_bar = (x_max + x_min) / 2
        else:
            x_bar = 0

        def targetFunc(x):
            return(0*x)
        self.deriv.append(targetFunc)

        for colNum in range(1, n+1):
            k = colNum

            def targetFunc(x, x_bar=x_bar, kappa=kappa,k=k):
                return (-k * kappa * np.sin((x - x_bar) * kappa * k))
            self.deriv.append(targetFunc)

            def targetFunc(x, x_bar=x_bar, kappa=kappa,k=k):
                return (k * kappa * np.cos((x - x_bar) * kappa * k))
            self.deriv.append(targetFunc)

def plotFourier(n,kappa,x_min,x_max,num_plotting_pts = 1000,centered=False,plotDeriv=False):

    if plotDeriv:
        basisFuncs = fourierBasisFunctionVec(n, kappa, x_min,x_max, centered).deriv
    else:
        basisFuncs = fourierBasisFunctionVec(n, kappa, x_min, x_max, centered).func

    x_colloc_pts = np.linspace(x_min, x_max, num_plotting_pts)

    fig = plt.figure(figsize=(18, 5))
    ax = plt.subplot(111)
    plt.grid()
    plt.title(f"Fourier Basis Function # Modes = {n}")
    plt.plot(x_colloc_pts,basisFuncs[0](x_colloc_pts),"k-",label="1")


    for k in range(1,n+1):
        colNum = k*2 - 1

        plt.plot(x_colloc_pts,basisFuncs[colNum](x_colloc_pts), label=r"cos(%i $\kappa$(x-$\overline{x}$))" %k)
        plt.plot(x_colloc_pts,basisFuncs[colNum+1](x_colloc_pts), label=r"sin(%i $\kappa$(x-$\overline{x}$))" %k)

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    # plt.legend()
    plt.show()

def genGramMtrxFourier(n, kappa, x_min, x_max, xList, yList, constrain=False, centered=False, constrainDeriv=False):
    '''
    Creates an nxn matrix of zeros then populates it. It begins iterating through each column,
    filling in all rows above and including the diagonal. Because the Gram matrix is symmetric,
    it then copies all entries across the diagonal.
    '''

    # Creates vector with all the basis functions [phi_0, phi_1, ..., phi_n]
    phiObj = fourierBasisFunctionVec(n, kappa, x_min,x_max, centered)
    phiVec = phiObj.func
    derivVec = phiObj.deriv

    # Initilize parameters relating to mode
    num_eqns = int(2*n+1)

    # Generate nxn Gram matrix
    gramMtrx = np.zeros((num_eqns, num_eqns))

    if centered:
        gramMtrx[0,0] = x_max-x_min + 2
    else:
        gramMtrx[0,0] = x_max-x_min

    if constrainDeriv:
        phi_i = phiVec[0]
        phi_i_prime = derivVec[0]
        phi_j = phiVec[0]
        phi_j_prime = phiVec[0]
        targetFunc = lambda x: phi_i(x) * phi_j(x)
        gramMtrx[0,0] = scipy.integrate.quad(targetFunc, x_min, x_max, limit=1000,
                            epsrel=10e-5, epsabs=10e-5)[0] + phi_i(x_min) * phi_j(x_min) + phi_i(x_max) * phi_j(x_max) \
                            + phi_i_prime(x_min) * phi_j_prime(x_min) + phi_i_prime(x_max) * phi_j_prime(x_max)

    for colNum in range(1, num_eqns):
        # Generates all exp(.)cos(.) * phi_i functions down the j'th column
        phi_j = phiVec[colNum]
        phi_j_prime = derivVec[colNum]
        for distToDiag in range(colNum + 1):
            phi_i = phiVec[distToDiag]
            phi_i_prime = derivVec[distToDiag]
            targetFunc = lambda x: phi_i(x) * phi_j(x)

            if constrain:
                gramMtrx[distToDiag, colNum] = scipy.integrate.quad(targetFunc, x_min, x_max, limit=1000,
                                                                 epsrel=10e-5, epsabs=10e-5)[0] + phi_i(x_min) * phi_j(
                    x_min) + phi_i(x_max) * phi_j(x_max)
            elif constrainDeriv:
                gramMtrx[distToDiag, colNum] = scipy.integrate.quad(targetFunc, x_min, x_max, limit=1000,
                            epsrel=10e-5, epsabs=10e-5)[0] + phi_i(x_min) * phi_j(x_min) + phi_i(x_max) * phi_j(x_max) \
                            + phi_i_prime(x_min) * phi_j_prime(x_min) + phi_i_prime(x_max) * phi_j_prime(x_max)
            else:
                gramMtrx[distToDiag, colNum] = scipy.integrate.quad(targetFunc, x_min, x_max, limit=1000, epsrel=10e-5,
                                                                 epsabs=10e-5)[0]

        """# Generates all exp(.)cos(.) * phi_i functions down the j'th column
        phi_j = phiVec[colNum+1]
        phi_j_prime = derivVec[colNum+1]
        for distToDiag in range(colNum + 2):
            phi_i = phiVec[distToDiag]
            phi_i_prime = derivVec[distToDiag]
            targetFunc = lambda x: phi_i(x) * phi_j(x)

            if constrain:
                gramMtrx[distToDiag, colNum+1] = scipy.integrate.quad(targetFunc, x_min, x_max, limit=1000,
                                                    epsrel=10e-5, epsabs=10e-5)[0] + phi_i(x_min) * phi_j(x_min)\
                                                 + phi_i(x_max) * phi_j(x_max)
            if constrainDeriv:
                gramMtrx[distToDiag, colNum] = scipy.integrate.quad(targetFunc, x_min, x_max, limit=1000,
                            epsrel=10e-5, epsabs=10e-5)[0] + phi_i(x_min) * phi_j(x_min) + phi_i(x_max) * phi_j(x_max) \
                            + phi_i_prime(x_min) * phi_j_prime(x_min) + phi_i_prime(x_max) * phi_j_prime(x_max)
            else:
                gramMtrx[distToDiag, colNum+1] = \
                scipy.integrate.quad(targetFunc, x_min, x_max, limit=1000, epsrel=10e-5,
                                     epsabs=10e-5)[0]"""

    # Copies elements of the gram matrix across the diagonal
    for rowNum in range(num_eqns):
        for i in range(rowNum):
            gramMtrx[rowNum, i] = gramMtrx[i, rowNum]

    return (gramMtrx)

def genLoadMtrxFourier(n, kappa, x_min, x_max, xList, yList, constrain=False, centered=False,constrainDeriv=False):
    # Creates vector with all the basis functions [phi_0, phi_1, ..., phi_n]
    phiObj = fourierBasisFunctionVec(n, kappa, x_min, x_max, centered)
    phiVec = phiObj.func
    derivVec = phiObj.deriv

    # Initilize parameters relating to mode
    num_eqns = int(2*n + 1)

    # Initialize load matrix
    loadMtrx = np.zeros(num_eqns)

    for colNum in range(num_eqns):
        phi_i = phiVec[colNum]
        phi_i_prime = derivVec[colNum]
        f = lambda x: genCollocationPts(xList, yList, np.asarray([x]))[0]
        f_prime = lambda x: genFuncDerivAtPoint(xList,yList,x)
        targetFunc = lambda x: f(x) * phi_i(x)

        if constrain:
            loadMtrx[colNum] = scipy.integrate.quad(targetFunc, x_min, x_max, limit=1000, epsrel=10e-5, epsabs=10e-5)[0] + phi_i(x_min) * f(x_min) + phi_i(x_max) * f(x_max)
        elif constrainDeriv:
            loadMtrx[colNum] = scipy.integrate.quad(targetFunc, x_min, x_max, limit=1000, epsrel=10e-5, epsabs=10e-5)[0] \
                        + phi_i(x_min) * f(x_min) + phi_i(x_max) * f(x_max) + phi_i_prime(x_min) * f_prime(x_min) \
                        + phi_i_prime(x_max) * f_prime(x_max)
        else:
            loadMtrx[colNum] = scipy.integrate.quad(targetFunc, x_min, x_max,  limit=1000, epsrel=10e-5, epsabs=10e-5)[0]

    return (loadMtrx)

def genThetaHatVec(k):
    """
    Generates a legendre parent function valid on the domain x_hat = [-1,1]
    Returns a function with arguments (x,degree)
    Uses recursion to generate legendre functions of all previous degrees, may be slow for high degree cases
    """
    thetaHatVec = []
    #k = degree
    if k == 0:
        def targetFunc(x,k=k):
            return(x*k*0 + 1)

    elif k == 1:
        def thetaHat1(zetta):
            return(- (zetta-1)/2)

        def thetaHat2(zetta):
            return((zetta+1)/2)
        thetaHatVec.append(thetaHat1)
        thetaHatVec.append(thetaHat2)

    elif k==2:
        def thetaHat1(zetta):
            return(-((zetta-1)/2)* (-zetta))

        def thetaHat2(zetta):
            return(zetta * (zetta+1)/2)

        def thetaHat3(zetta):
            return(-(zetta-1)*(zetta+1))

        thetaHatVec.append(thetaHat1)
        thetaHatVec.append(thetaHat2)
        thetaHatVec.append(thetaHat3)

    return (thetaHatVec)

def legendreFunc(k):
    """
    Generates a legendre function valid on the domain x_hat = [-1,1]
    Returns a function with arguments (x,degree)
    Uses recursion to generate legendre functions of all previous degrees, may be slow for high degree cases
    """
    #k = degree
    if k == 0:
        def targetFunc(x,k=k):
            return(x*k*0 + 1)

    elif k == 1:
        def targetFunc(x,k=k):
            return(x + 0*k)

    else:
        def targetFunc(x,k=k):
            return(((2*k - 1)/k) * x * legendreFunc(k-1)(x,k-1) - ((k-1)/k) * legendreFunc(k-2)(x,k-2))

    return (targetFunc)

def genPartitions(x_min, x_max, locList, numPartsList):
    """
    Splits a domain into partitions and those partitions into sub-partitions.
    Whole domain is given as between x_min and x_max.

    locList = list of x locations where primary positions are located
    numPartsList = number of sub-partitions to be generated for each primary partition
    """

    partList = []

    partList.append(np.linspace(x_min, locList[0], numPartsList[0] + 1))
    for i, loc in np.ndenumerate(locList[:-1]):
        partList.append(np.linspace(loc, locList[i[0] + 1], numPartsList[i[0] + 1] + 1))

    if locList[-1] != x_max:

        partList.append(np.linspace(locList[-1], x_max, numPartsList[-1] + 1))
    return (partList)

def genGlobalNodes(x_min, x_max, locList, numPartsList):
    partList = genPartitions(x_min, x_max, locList, numPartsList)
    partList = [item for elem in partList for item in elem]
    partList = set(partList)
    partList = list(partList)
    partList.sort()
    return(partList)




#######################################################################################################################
#######################################################################################################################
"""
The code within this commented block were written by Prof. Valmor de Almeida, Fall 2021, ENGY-5310 at UMass Lowell
Expanded by Andrew Hamel Fall 2021 
"""
'''Domain partition'''


"""def get_domain_partition(degree, locList, numPartsList, x_min, x_max, bc_x_min='flux', bc_x_max='flux'):
    partList = genGlobalNodes(x_min, x_max, locList, numPartsList)
    n_elem = len(partList) - 1

    if degree == 2:
        # Local node numbering on parent domain
        # --0--------------1---->
        #  -1      0      +1    zeta
        gnodes_x = []
        for i, part in enumerate(partList[1:]):
            h_e = abs(part - partList[i])
            gnodes_x.append(partList[i])
            gnodes_x.append(partList[i] + h_e * (1 / 2))
        gnodes_x.append(x_max)

        gnodes_x = np.asarray(gnodes_x)

        patches = list()
        local_to_global_node_id_map = list()
        for e in range(n_elem):
            gnode_id_2 = 2 * e + 2  # right
            gnode_id_3 = 2 * e + 1
            gnode_id_1 = 2 * e  # left

            x1 = gnodes_x[gnode_id_1]
            x2 = gnodes_x[gnode_id_2]

            # Local node id:  0   1
            patches.append((x1, x2))
            # Local node id:                        0           1
            local_to_global_node_id_map.append([gnode_id_1, gnode_id_3, gnode_id_2])
        if bc_x_min == 'essential':
            local_to_global_node_id_map[0][0] = -1
        if bc_x_max == 'essential':
            local_to_global_node_id_map[-1][-1] = -1
    if degree == 1:
        gnodes_x = np.linspace(x_min, x_max, n_elem + 1, dtype=np.float64)
        patches = list()
        local_to_global_node_id_map = list()
        for e in range(n_elem):
            gnode_id_2 = e + 1  # right
            gnode_id_1 = gnode_id_2 - 1  # left
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

    return (patches, gnodes_x, local_to_global_node_id_map)"""


def get_domain_partition(degree, locList, numPartsList, x_min, x_max, bc_x_min='flux', bc_x_max='flux'):
    partList = genGlobalNodes(x_min, x_max, locList, numPartsList)
    n_elem = len(partList) - 1
    numLocalNodes = degree + 1
    gnodes_x = []
    for i, part in enumerate(partList[1:]):
        h_e = abs(part - partList[i])
        for localNode in range(degree):
            gnodes_x.append(partList[i] + h_e * (localNode / degree))

    gnodes_x.append(x_max)

    gnodes_x = np.asarray(gnodes_x)

    patches = list()
    local_to_global_node_id_map = list()
    for e in range(n_elem):
        gnodeList = [degree * e, degree * e + degree]

        for localNode in range(1, degree):
            gnodeList.append(degree * e + localNode)

        x1 = gnodes_x[gnodeList[0]]
        x2 = gnodes_x[gnodeList[1]]

        # Local node id:  0   1
        patches.append((x1, x2))

        # Local node id:                        0           1
        local_to_global_node_id_map.append([x for x in range(degree * e, degree * e + degree + 1)])
    if bc_x_min == 'dirichlet':
        local_to_global_node_id_map[0][0] = -1
    if bc_x_max == 'dirichlet':
        local_to_global_node_id_map[-1][-1] = -1


    return (patches, gnodes_x, local_to_global_node_id_map)

'''Parent mapping'''
def get_parent_mapping():
    # zeta in [-1,1]
    parent_mapping = lambda zeta, x_e_bar, h_e: x_e_bar + h_e/2 * zeta # compute x
    parent_mapping_prime = lambda h_e: h_e/2                           # compute mapping derivative wrt zeta
    # x in Omega_e
    inverse_parent_mapping = lambda x, x_e_bar, h_e: (x - x_e_bar)*2/h_e # compute zeta
    inverse_parent_mapping_der = lambda h_e: 2/h_e
    return (parent_mapping, parent_mapping_prime, inverse_parent_mapping,inverse_parent_mapping_der)


'''Parent basis functions'''

"""def get_parent_basis_functions(degree):
    parent_basis_func_list = list()
    parent_basis_func_prime_list = list()
    if degree == 1:
        parent_basis_func_list.append(lambda zetta: -(zetta - 1) / 2)  # left
        parent_basis_func_list.append(lambda zetta: (zetta + 1) / 2)  # right
        parent_basis_func_prime_list.append(lambda zetta: -1 / 2)  # left
        parent_basis_func_prime_list.append(lambda zetta: 1 / 2)  # right
    elif degree == 2:
        parent_basis_func_list.append(lambda zetta: -((zetta - 1) / 2) * (-zetta))
        parent_basis_func_list.append(lambda zetta: -(zetta - 1) * (zetta + 1))
        parent_basis_func_list.append(lambda zetta: (zetta + 1) / 2 * zetta)

        parent_basis_func_prime_list.append(lambda zetta: (2 * zetta - 1) / 2)
        parent_basis_func_prime_list.append(lambda zetta: -2 * zetta)
        parent_basis_func_prime_list.append(lambda zetta: (2 * zetta + 1) / 2)

    return (parent_basis_func_list, parent_basis_func_prime_list)



    return (parent_basis_func_list, parent_basis_func_prime_list)"""
def get_parent_basis_functions(degree):
    num_funcs = degree + 1
    local_nodes = np.linspace(-1, 1, num_funcs)
    phi_vec = []
    phi_prime_vec = []
    for i in range(num_funcs):
        selectedNodeList = np.zeros(num_funcs)
        selectedNodeList[i] = 1

        phi = scipy.interpolate.lagrange(local_nodes, selectedNodeList)
        phi_prime = np.polyder(phi)

        phi_vec.append(phi)
        phi_prime_vec.append(phi_prime)



    sorted_phi_vec = [phi_vec[0], phi_vec[-1]]
    sorted_phi_prime_vec = [phi_prime_vec[0], phi_prime_vec[-1]]

    for i,phi in enumerate(phi_vec[1:-1]):
        sorted_phi_vec.append(phi)
        sorted_phi_prime_vec.append(phi_prime_vec[i+1])
    
    parent_basis_func_list = phi_vec
    parent_basis_func_prime_list = phi_prime_vec

    return (parent_basis_func_list, parent_basis_func_prime_list)



'''Any global basis function'''


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
                        phi_i_x[j] = parent_basis_func_list[I](zetta)
                        #phi_prime_i_x[j] = parent_basis_deriv_list[I](zetta)
                break
    return phi_i_x

def global_basis_deriv(i, x, domain_partition, parent_mapping, parent_basis_functions):
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

    phi_i_x = np.copy(x) * 0.0  # initialization
    phi_prime_i_x = np.copy(x) * 0.0  # initialization

    patches = domain_partition[0]
    local_to_global_node_id_map = domain_partition[2]
    parent_mapping_g2p = parent_mapping[0]
    parent_mapping_prime = parent_mapping[1]
    inverse_parent_mapping = parent_mapping[2]
    inverse_parent_mapping_der = parent_mapping[3]

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
                        #zetta = inverse_parent_mapping(x_j,x_e_bar,h_e)/parent_mapping_prime(inverse_parent_mapping(x_j, x_e_bar, h_e))
                        #zetta = inverse_parent_mapping_der(h_e)
                        zetta = inverse_parent_mapping(x_j,x_e_bar, h_e)
                        #print(x_j,x_e_bar,h_e)

                        #print(zetta)
                        #print("*******")
                        #zetta = x_j
                        #phi_prime_i_x[j] = parent_basis_deriv_list[I](x)/ parent_mapping_prime(inverse_parent_mapping(x_j,x_e_bar, h_e))
                        phi_prime_i_x[j] = parent_basis_deriv_list[I](zetta)
                break
    return phi_prime_i_x

def global_basis_deriv_forIntegrating(i, x, domain_partition, parent_mapping, parent_basis_functions):
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

    phi_i_x = np.copy(x) * 0.0  # initialization
    phi_prime_i_x = np.copy(x) * 0.0  # initialization

    patches = domain_partition[0]
    local_to_global_node_id_map = domain_partition[2]
    parent_mapping_g2p = parent_mapping[0]
    parent_mapping_prime = parent_mapping[1]
    inverse_parent_mapping = parent_mapping[2]
    inverse_parent_mapping_der = parent_mapping[3]

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
                        zetta = inverse_parent_mapping(x_j,x_e_bar, h_e)
                        phi_prime_i_x[j] = parent_basis_deriv_list[I](zetta) * 2/h_e
                break
    return phi_prime_i_x
'''All global basis functions'''


def get_global_basis_functions(domain_partition, parent_mapping, parent_basis_functions, global_basis_function):
    basis_func_list = list()
    n_gnodes = domain_partition[1].size
    local_to_global_node_id_map = domain_partition[2]
    phi_i = lambda i, x: global_basis_function(i, x, domain_partition, parent_mapping, parent_basis_functions)

    visited = [False] * n_gnodes
    n_elem = len(domain_partition[0])
    for e in range(n_elem):
        for I in range(len(local_to_global_node_id_map[e])):
            gnode_id = local_to_global_node_id_map[e][I]
            if gnode_id >= 0 and not visited[gnode_id]:
                basis_func_list.append(lambda x, i=gnode_id: phi_i(i, x))
                visited[gnode_id] = True

    assert len(basis_func_list) >= 1, 'There are no basis functions to build.'

    return basis_func_list


def get_global_basis_deriv(domain_partition, parent_mapping, parent_basis_functions, global_basis_deriv):
    basis_func_list = list()
    n_gnodes = domain_partition[1].size
    local_to_global_node_id_map = domain_partition[2]
    phi_i = lambda i, x: global_basis_deriv(i, x, domain_partition, parent_mapping, parent_basis_functions)

    visited = [False] * n_gnodes
    n_elem = len(domain_partition[0])
    for e in range(n_elem):
        for I in range(len(local_to_global_node_id_map[e])):
            gnode_id = local_to_global_node_id_map[e][I]
            if gnode_id >= 0 and not visited[gnode_id]:
                basis_func_list.append(lambda x, i=gnode_id: phi_i(i, x))
                visited[gnode_id] = True

    assert len(basis_func_list) >= 1, 'There are no basis functions to build.'

    return basis_func_list

def get_global_basis_deriv_forIntegrating(domain_partition, parent_mapping, parent_basis_functions, global_basis_deriv_forIntegrating):
    basis_func_list = list()
    n_gnodes = domain_partition[1].size
    local_to_global_node_id_map = domain_partition[2]
    phi_i = lambda i, x: global_basis_deriv_forIntegrating(i, x, domain_partition, parent_mapping, parent_basis_functions)

    visited = [False] * n_gnodes
    n_elem = len(domain_partition[0])
    for e in range(n_elem):
        for I in range(len(local_to_global_node_id_map[e])):
            gnode_id = local_to_global_node_id_map[e][I]
            if gnode_id >= 0 and not visited[gnode_id]:
                basis_func_list.append(lambda x, i=gnode_id: phi_i(i, x))
                visited[gnode_id] = True

    assert len(basis_func_list) >= 1, 'There are no basis functions to build.'

    return basis_func_list
'''Plot global basis functions'''


def plot_func(domain_partition,phi_list, x_min,x_max,title='Lagrange Basis Functions'):
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

def inner_product(u, v, patches):
    integrand = lambda x: u(x) * v(x)
    inner_product = 0.0
    for nodes_x in patches:
        (inner_product_e, _) = scipy.integrate.quad(integrand, nodes_x[0], nodes_x[1],limit=1000, epsrel=10e-7, epsabs=10e-7)
        inner_product += inner_product_e
    return inner_product

def get_gram_matrix(basis_func_list, domain_partition, inner_product,constrain=False):
    N = len(basis_func_list)
    patches = domain_partition[0]
    x_min = domain_partition[1][0]
    x_max = domain_partition[1][-1]

    a_mtrx = np.zeros((N, N), dtype=np.float64)
    for i,phi_i in enumerate(basis_func_list):
        for j,phi_j in enumerate(basis_func_list):
            if constrain:
                a_mtrx[i, j] = inner_product(phi_i, phi_j, patches) + phi_i(x_min)*phi_j(x_min) + phi_i(x_max)*phi_j(x_max)
            else:
                a_mtrx[i, j] = inner_product(phi_i, phi_j, patches)
    return a_mtrx


def assemble_a_mtrx(femlb):
    n_dof = femlb.n_dof  # total number of degrees of freedom
    a_mtrx = np.zeros((n_dof + femlb.degree, n_dof + femlb.degree), dtype=np.float64)

    patches = femlb.domain_partition[0]
    n = len(femlb.parent_basis_func_list)

    local_a_diff_mtrx = np.zeros((n, n), dtype="object")
    local_a_source_mtrx = np.zeros((n, n), dtype="object")

    for I in range(n):
        parent_basis_func_prime_I = femlb.parent_basis_func_prime_list[I]
        parent_basis_func_I = femlb.parent_basis_func_list[I]
        for J in range(n):

            integrand = lambda zeta, parent_basis_func_prime_I=parent_basis_func_prime_I: parent_basis_func_prime_I(
                zeta) * femlb.parent_basis_func_prime_list[J](zeta)
            local_a_diff_mtrx[I, J] = integrand
            integrand = lambda zeta, parent_basis_func_I=parent_basis_func_I: - parent_basis_func_I(zeta) * \
                                                                              femlb.parent_basis_func_list[J](zeta)
            local_a_source_mtrx[I, J] = integrand

    for (e, dof_ids) in enumerate(femlb.domain_partition[2]):

        patch_nodes_x = patches[e]

        h_e = patch_nodes_x[1] - patch_nodes_x[0]
        x_e_bar = (patch_nodes_x[0] + patch_nodes_x[1]) / 2.0
        parent_mapping_jacobian = femlb.parent_mapping_prime(h_e)

        diff_coeff_parent = lambda zeta, x_e_bar=x_e_bar, h_e=h_e: \
            femlb.diff_coeff(femlb.parent_mapping[0](zeta, x_e_bar, h_e))
        source_slope_parent = lambda zeta, x_e_bar=x_e_bar, h_e=h_e: \
            femlb.source_slope(femlb.parent_mapping[0](zeta, x_e_bar, h_e))

        source_bias_parent = lambda zeta, x_e_bar=x_e_bar, h_e=h_e: \
            femlb.source_bias(femlb.parent_mapping[0](zeta, x_e_bar, h_e))

        w_lift_prime_parent = lambda zeta, x_e_bar=x_e_bar, h_e=h_e: \
            femlb.w_lift_prime(femlb.parent_mapping[0](zeta, x_e_bar, h_e)) * h_e / 2

        w_lift_parent = lambda zeta, x_e_bar=x_e_bar, h_e=h_e: \
            femlb.w_lift(femlb.parent_mapping[0](zeta, x_e_bar, h_e))
        for (I, i) in enumerate(dof_ids):
            if i < 0: continue
            for (J, j) in enumerate(dof_ids):
                if j < 0: continue
                a_mtrx[i, j] += \
                scipy.integrate.quad(lambda zeta: local_a_diff_mtrx[I, J](zeta) * diff_coeff_parent(zeta), -1, 1,
                                     limit=1000, epsrel=10e-7, epsabs=10e-7)[0] / parent_mapping_jacobian + \
                scipy.integrate.quad(lambda zeta: local_a_source_mtrx[I, J](zeta) * source_slope_parent(zeta), -1, 1,
                                     limit=1000, epsrel=10e-7, epsabs=10e-7)[0] * parent_mapping_jacobian

    if femlb.bc_x_min == femlb.bc_x_max == "dirichlet":
        a_mtrx = a_mtrx[1:-1, 1:-1]

    return a_mtrx

def assemble_a_mtrx_neumann_robin(femlb):
    n_dof = femlb.n_dof  # total number of degrees of freedom
    a_mtrx = np.zeros((n_dof , n_dof), dtype=np.float64)

    patches = femlb.domain_partition[0]
    n = len(femlb.parent_basis_func_list)

    local_a_diff_mtrx = np.zeros((n, n), dtype="object")
    local_a_source_mtrx = np.zeros((n, n), dtype="object")

    for I in range(n):
        parent_basis_func_prime_I = femlb.parent_basis_func_prime_list[I]
        parent_basis_func_I = femlb.parent_basis_func_list[I]
        for J in range(n):

            integrand = lambda zeta, parent_basis_func_prime_I=parent_basis_func_prime_I: parent_basis_func_prime_I(
                zeta) * femlb.parent_basis_func_prime_list[J](zeta)
            local_a_diff_mtrx[I, J] = integrand
            integrand = lambda zeta, parent_basis_func_I=parent_basis_func_I: - parent_basis_func_I(zeta) * \
                                                                              femlb.parent_basis_func_list[J](zeta)
            local_a_source_mtrx[I, J] = integrand

    for (e, dof_ids) in enumerate(femlb.domain_partition[2]):

        patch_nodes_x = patches[e]

        h_e = patch_nodes_x[1] - patch_nodes_x[0]
        x_e_bar = (patch_nodes_x[0] + patch_nodes_x[1]) / 2.0
        parent_mapping_jacobian = femlb.parent_mapping_prime(h_e)

        diff_coeff_parent = lambda zeta, x_e_bar=x_e_bar, h_e=h_e: \
            femlb.diff_coeff(femlb.parent_mapping[0](zeta, x_e_bar, h_e))
        source_slope_parent = lambda zeta, x_e_bar=x_e_bar, h_e=h_e: \
            femlb.source_slope(femlb.parent_mapping[0](zeta, x_e_bar, h_e))

        source_bias_parent = lambda zeta, x_e_bar=x_e_bar, h_e=h_e: \
            femlb.source_bias(femlb.parent_mapping[0](zeta, x_e_bar, h_e))

        w_lift_prime_parent = lambda zeta, x_e_bar=x_e_bar, h_e=h_e: \
            femlb.w_lift_prime(femlb.parent_mapping[0](zeta, x_e_bar, h_e)) * h_e / 2

        w_lift_parent = lambda zeta, x_e_bar=x_e_bar, h_e=h_e: \
            femlb.w_lift(femlb.parent_mapping[0](zeta, x_e_bar, h_e))
        for (I, i) in enumerate(dof_ids):
            if i < 0: continue
            for (J, j) in enumerate(dof_ids):
                if j < 0: continue
                a_mtrx[i, j] += \
                scipy.integrate.quad(lambda zeta: local_a_diff_mtrx[I, J](zeta) * diff_coeff_parent(zeta), -1, 1,
                                     limit=1000, epsrel=10e-7, epsabs=10e-7)[0] / parent_mapping_jacobian + \
                scipy.integrate.quad(lambda zeta: local_a_source_mtrx[I, J](zeta) * source_slope_parent(zeta), -1, 1,
                                     limit=1000, epsrel=10e-7, epsabs=10e-7)[0] * parent_mapping_jacobian
                
                if  (j == femlb.n_dof-1 and e==len(femlb.domain_partition[2])-1):
                    a_mtrx[i,j] += (femlb.h * femlb.parent_basis_func_list[J](1)*femlb.parent_basis_func_list[I](1))

    return a_mtrx

def assemble_b_vec(femlb):
    n_dof = femlb.n_dof  # total number of degrees of freedom
    b_vec = np.zeros(n_dof + 2, dtype=np.float64)

    patches = femlb.domain_partition[0]

    n = len(femlb.parent_basis_func_list)

    local_b_bias_vec = np.zeros(n, dtype=np.float64)
    local_b_diff_vec = np.zeros(n, dtype=np.float64)
    local_b_source_vec = np.zeros(n, dtype=np.float64)

    for (e, dof_ids) in enumerate(femlb.domain_partition[2]):
        patch_nodes_x = patches[e]
        h_e = patch_nodes_x[1] - patch_nodes_x[0]
        x_e_bar = (patch_nodes_x[0] + patch_nodes_x[1]) / 2.0
        parent_mapping_jacobian = femlb.parent_mapping_prime(h_e)

        source_bias_parent = lambda zeta, x_e_bar=x_e_bar, h_e=h_e: \
            femlb.source_bias(femlb.parent_mapping[0](zeta, x_e_bar, h_e))

        w_lift_prime_parent = lambda zeta, x_e_bar=x_e_bar, h_e=h_e: \
            femlb.w_lift_prime(femlb.parent_mapping[0](zeta, x_e_bar, h_e)) * h_e / 2

        w_lift_parent = lambda zeta, x_e_bar=x_e_bar, h_e=h_e: \
            femlb.w_lift(femlb.parent_mapping[0](zeta, x_e_bar, h_e))

        for I in range(n):
            parent_basis_func_I = femlb.parent_basis_func_list[I]
            parent_basis_func_prime_I = femlb.parent_basis_func_prime_list[I]

            integrand = lambda zeta: source_bias_parent(zeta) * parent_basis_func_I(zeta)
            (local_b_bias_vec[I], _) =scipy.integrate.quad(integrand, -1, 1)

            integrand = lambda zeta: - femlb.diff_coeff(
                femlb.parent_mapping[0](zeta, x_e_bar, h_e)) * w_lift_prime_parent(zeta) * parent_basis_func_prime_I(
                zeta)
            (local_b_diff_vec[I], _) =scipy.integrate.quad(integrand, -1, 1, limit=1000, epsrel=10e-7, epsabs=10e-7)

            integrand = lambda zeta: femlb.source_slope(femlb.parent_mapping[0](zeta, x_e_bar, h_e)) * w_lift_parent(
                zeta) * parent_basis_func_I(zeta)
            (local_b_source_vec[I], _) =scipy.integrate.quad(integrand, -1, 1, limit=1000, epsrel=10e-7, epsabs=10e-7)

        for (I, i) in enumerate(dof_ids):
            if i < 0: continue
            b_vec[i] += local_b_bias_vec[I] * parent_mapping_jacobian + \
                        local_b_diff_vec[I] / parent_mapping_jacobian + \
                        local_b_source_vec[I] * parent_mapping_jacobian
    if femlb.bc_x_min == femlb.bc_x_max == "dirichlet":
        b_vec = b_vec[1:-1]

    return b_vec

def assemble_b_vec_neumann_robin(femlb):
    n_dof = femlb.n_dof  # total number of degrees of freedom
    b_vec = np.zeros(n_dof, dtype=np.float64)

    patches = femlb.domain_partition[0]

    n = len(femlb.parent_basis_func_list)

    local_b_bias_vec = np.zeros(n, dtype=np.float64)
    local_b_diff_vec = np.zeros(n, dtype=np.float64)
    local_b_source_vec = np.zeros(n, dtype=np.float64)

    for (e, dof_ids) in enumerate(femlb.domain_partition[2]):
        patch_nodes_x = patches[e]
        h_e = patch_nodes_x[1] - patch_nodes_x[0]
        x_e_bar = (patch_nodes_x[0] + patch_nodes_x[1]) / 2.0
        parent_mapping_jacobian = femlb.parent_mapping_prime(h_e)

        source_bias_parent = lambda zeta, x_e_bar=x_e_bar, h_e=h_e: \
            femlb.source_bias(femlb.parent_mapping[0](zeta, x_e_bar, h_e))

        w_lift_prime_parent = lambda zeta, x_e_bar=x_e_bar, h_e=h_e: \
            femlb.w_lift_prime(femlb.parent_mapping[0](zeta, x_e_bar, h_e)) * h_e / 2

        w_lift_parent = lambda zeta, x_e_bar=x_e_bar, h_e=h_e: \
            femlb.w_lift(femlb.parent_mapping[0](zeta, x_e_bar, h_e))

        for I in range(n):
            parent_basis_func_I = femlb.parent_basis_func_list[I]
            parent_basis_func_prime_I = femlb.parent_basis_func_prime_list[I]

            integrand = lambda zeta: source_bias_parent(zeta) * parent_basis_func_I(zeta)
            (local_b_bias_vec[I], _) =scipy.integrate.quad(integrand, -1, 1)

            integrand = lambda zeta: - femlb.diff_coeff(
                femlb.parent_mapping[0](zeta, x_e_bar, h_e)) * w_lift_prime_parent(zeta) * parent_basis_func_prime_I(
                zeta)
            (local_b_diff_vec[I], _) =scipy.integrate.quad(integrand, -1, 1, limit=1000, epsrel=10e-7, epsabs=10e-7)

            integrand = lambda zeta: femlb.source_slope(femlb.parent_mapping[0](zeta, x_e_bar, h_e)) * w_lift_parent(
                zeta) * parent_basis_func_I(zeta)
            (local_b_source_vec[I], _) =scipy.integrate.quad(integrand, -1, 1, limit=1000, epsrel=10e-7, epsabs=10e-7)

        for (I, i) in enumerate(dof_ids):
            if i < 0: continue
            b_vec[i] += local_b_bias_vec[I] * parent_mapping_jacobian
            if i ==0:
                b_vec[i] += femlb.x_min_flux
            if i == len(b_vec)-1:
                b_vec[i] += femlb.h * femlb.x_max_u_ref

    return b_vec

def get_b_vec_inhomogeneous_bc(femlb):
    n = len(femlb.basis_func_list)

    b_vec = np.zeros(n, dtype=np.float64)
    xpts = np.linspace(femlb.x_a, femlb.x_b, 500)
    for i, phi_i in enumerate(femlb.basis_func_list):
        (b_vec[i], _) = femlb.inner_product(femlb.source_bias, phi_i)

        d_x_w_prime = lambda x: femlb.diff_coeff(x) * femlb.w_lift_prime(x)
        (term1, _) = femlb.inner_product(d_x_w_prime, femlb.basis_func_prime_list[i])
        b_vec[i] -= term1

        s_x_w = lambda x: femlb.source_slope(x) * femlb.w_lift(x)
        (term2, _) = femlb.inner_product(s_x_w, phi_i)
        b_vec[i] += term2

    return b_vec

def assemble_gram_matrix(a_mtrx, domain_partition, parent_mapping, parent_basis_func_list):
    n = len(parent_basis_func_list[0])

    local_to_global_node_id_map = domain_partition[2]

    local_a_mtrx = np.zeros((n, n), dtype=np.float64)

    for I in range(n):
        parent_basis_func_I = parent_basis_func_list[0][I]
        for J in range(n):
            integrand = lambda zeta: parent_basis_func_I(zeta) * parent_basis_func_list[0][J](zeta)
            (local_a_mtrx[I, J], _) = scipy.integrate.quad(integrand, -1, 1,limit=1000, epsrel=10e-7, epsabs=10e-7)
    for (e, gnode_ids) in enumerate(local_to_global_node_id_map):
        (x_0, x_1) = domain_partition[0][e]
        h_e = x_1 - x_0
        parent_mapping_jacobian = parent_mapping[1](h_e)
        for (I, i) in enumerate(gnode_ids):
            for (J, j) in enumerate(gnode_ids):
                 a_mtrx[i, j] += local_a_mtrx[I, J] * parent_mapping_jacobian
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
            (local_b_vec[I], _) = scipy.integrate.quad(integrand, -1, 1,limit=1000, epsrel=10e-7, epsabs=10e-7)


        for (I, i) in enumerate(gnode_ids):
            b_vec[i] += local_b_vec[I] * parent_mapping_jacobian


def g_Legendre(x, phi_list, c_star_vec):
    eval_mtrx = np.zeros(len(phi_list), dtype=np.float64)

    for (j, phi_j) in enumerate(phi_list):
        eval_mtrx[j] = phi_j(x)
    g_best = eval_mtrx @ c_star_vec
    return (g_best)


def g_LegendreSquared(x, phi_list, c_star_vec):
    eval_mtrx = np.zeros(len(phi_list), dtype=np.float64)
    for (j, phi_j) in enumerate(phi_list):
        eval_mtrx[j] = phi_j(x)
    g_best = eval_mtrx @ c_star_vec
    return (g_best ** 2)


def f_Squared(x, f):
    return (f(x) ** 2)


def g_minus_f_squared(x, g_Legendre, f, phi_list, c_star_vec):
    gVal = g_Legendre(x, phi_list, c_star_vec)
    fVal = f(x)

    result = (gVal - fVal) ** 2
    return (result)

def get_load_vector(phi_list,domain_partition,f,constrain=False):
    b_vec = np.zeros(len(phi_list), dtype=np.float64)
    patches = domain_partition[0]
    x_min = domain_partition[1][0]
    x_max = domain_partition[1][-1]

    for i, phi_i in enumerate(phi_list):
        if constrain:
            b_vec[i] = inner_product(f, phi_i, patches) + f(x_min)*phi_i(x_min) + f(x_max)*phi_i(x_max)
        else:
            b_vec[i] = inner_product(f, phi_i, patches)
    return(b_vec)

def get_a_mtrx(femlb):
    n = len(femlb.basis_func_list)

    a_mtrx = np.zeros((n, n), dtype=np.float64)

    for i, phi_prime_i in enumerate(femlb.basis_func_prime_list_forIntegrating):
        for j, phi_prime_j in enumerate(femlb.basis_func_prime_list_forIntegrating):
            d_x_phi_prime_j = lambda x: femlb.diff_coeff(x) * phi_prime_i(x)

            (a_ij, _) = femlb.inner_product(d_x_phi_prime_j, phi_prime_j)

            a_mtrx[i, j] = a_ij

    for i, phi_i in enumerate(femlb.basis_func_list):
        for j, phi_j in enumerate(femlb.basis_func_list):
            s_x_phi_j = lambda x: femlb.source_slope(x) * phi_j(x)
            (a_ij, _) = femlb.inner_product(s_x_phi_j, phi_i)

            a_mtrx[i, j] -= a_ij

    return a_mtrx
#######################################################################################################################
#######################################################################################################################

class femLagrange:
    """
    Creates a vector of wavelet basis functions that can be called as a list with the .func method.
    Really could just be done with just a list, but I had already built this before I realized that...
    """

    def __init__(self, x_a,x_b,degree,liftDegree,f_shape_pts,k_shape_pts,s_shape_pts,w_shape_pts,domain_partition,u_a,u_b,u_a_0,u_b_0,bc_x_min="flux",bc_x_max="flux",h=0, x_min_flux = 0, x_max_u_ref=0, s_func=False):
        self.x_a = x_a
        self.x_b = x_b
        self.u_a = u_a
        self.u_b = u_b
        self.u_a_0 = u_a_0
        self.u_b_0 = u_b_0
        self.h = h
        self.x_min_flux = x_min_flux
        self.x_max_u_ref = x_max_u_ref

        self.w_shape_pts = w_shape_pts

        self.bc_x_min = bc_x_min
        self.bc_x_max = bc_x_max
        self.parent_basis_functions = get_parent_basis_functions(degree)

        self.parent_basis_func_list = get_parent_basis_functions(degree)[0]
        self.parent_basis_func_prime_list = get_parent_basis_functions(degree)[1]

        self.domain_partition = domain_partition

        self.parent_mapping = get_parent_mapping()
        self.parent_mapping_prime = self.parent_mapping[1]
        self.inv_parent_mapping_prime = self.parent_mapping[3]
        self.phi_list = get_global_basis_functions(domain_partition, self.parent_mapping, self.parent_basis_functions, global_basis_function)
        self.phi_prime_list = get_global_basis_deriv(domain_partition, self.parent_mapping, self.parent_basis_functions, global_basis_deriv)
        self.phi_prime_list_forIntegrating = get_global_basis_deriv_forIntegrating(domain_partition, self.parent_mapping, self.parent_basis_functions, global_basis_deriv_forIntegrating)

        self.liftDegree = liftDegree
        self.degree = degree
        self.n_pts = 500


        self.degree_partition = [degree]

        self.basis_func_list = self.phi_list
        self.basis_func_prime_list = self.phi_prime_list
        self.basis_func_prime_list_forIntegrating = self.phi_prime_list_forIntegrating


        self.n_dof = len(self.basis_func_list)
        self.liftDomainPartition = get_domain_partition(self.liftDegree, [x_b], [1], x_a, x_b)
        self.lift_parent_basis_functions = get_parent_basis_functions(liftDegree)
        self.lift_basis_func_list = get_global_basis_functions(self.liftDomainPartition, self.parent_mapping,self.lift_parent_basis_functions,global_basis_function)
        self.lift_basis_func_prime_list = get_global_basis_deriv(self.liftDomainPartition, self.parent_mapping, self.lift_parent_basis_functions, global_basis_deriv)


        self.diff_coeff = self.genTargetFunction(k_shape_pts)
        if s_func:
            self.source_slope = s_func
        else:
            self.source_slope = self.genTargetFunction(s_shape_pts)
        self.source_bias = self.genTargetFunction(f_shape_pts)

        self.parent_diff_coeff = self.genTargetFunctionZeta(k_shape_pts)
        self.parent_source_slope = self.genTargetFunctionZeta(s_shape_pts)
        self.parent_source_bias = self.genTargetFunctionZeta(f_shape_pts)


        self.alpha_vec = self.genAlphaVec()

        self.evaluation_matrix = self.genEvalMtrx
        self.lift_evaluation_mtrx = self.genLiftEvalMtrx


        self.w_lift = self.genWLift
        self.w_lift_prime = self.genWLiftPrime


    def genAlphaVec(self):
        alpha_vec = np.zeros(len(self.lift_basis_func_list))

        alpha_vec[0] = self.w_shape_pts[0][1]

        alpha_vec[-1] = self.w_shape_pts[1][1]

        return(alpha_vec)

    def genWLift(self,x):
        return(self.lift_evaluation_mtrx(x)@self.alpha_vec)

    def genWLiftPrime(self,x):
        a_mtrx = self.lift_evaluation_mtrx(x, derivative=True)

        return(a_mtrx @ self.alpha_vec)

    def inner_product(self, u, v):
        integrand = lambda x: u(x) * v(x)
        (inner_product_e) = scipy.integrate.quad(integrand, self.x_a,self.x_b,limit=1000, epsrel=10e-7, epsabs=10e-7)

        return(inner_product_e)


    def genTargetFunction(self,shape_pts,deriv=False):
        (x, y) = convertToNumpy(shape_pts, twoOutputs=True)
        if deriv:
            def targetFunction(xVal, x=x, y=y):
                if isinstance(xVal, np.ndarray):
                    return (genCollocationPtsDeriv(x, y, xVal))
                else:
                    xVal = np.asarray([xVal])
                    return (genCollocationPtsDeriv(x, y, xVal)[0])
            return(targetFunction)

        else:
            def targetFunction(xVal, x=x, y=y):
                if isinstance(xVal, np.ndarray):
                    return (genCollocationPts(x, y, xVal))
                else:
                    xVal = np.asarray([xVal])
                    return (genCollocationPts(x, y, xVal)[0])
            return(targetFunction)



    def genTargetFunctionZeta(self,shape_pts):
        (x, y) = convertToNumpy(shape_pts, twoOutputs=True)
        h_e = x[-1] - x[0]
        x_bar = (x[-1] + x[0]) / 2
        x = (x - x_bar) * 2/h_e
        def targetFunction(xVal, x=x, y=y):
            if isinstance(xVal, np.ndarray):
                return (genCollocationPts(x, y, xVal))
            else:
                xVal = np.asarray([xVal])
                return (genCollocationPts(x, y, xVal)[0])
        return(targetFunction)

    def genEvalMtrx(self,x_pts,derivative=False):
        if not(isinstance(x_pts, np.ndarray)):
            x_pts = np.asarray([x_pts])
        eval_mtrx = np.zeros((len(x_pts), len(self.basis_func_list)), dtype=np.float64)
        if derivative:
            for (j, phi_j) in enumerate(self.basis_func_prime_list_forIntegrating):
                eval_mtrx[:, j] = phi_j(x_pts)
        else:
            for (j, phi_j) in enumerate(self.basis_func_list):
                eval_mtrx[:, j] = phi_j(x_pts)

        return(eval_mtrx)

    def genLiftEvalMtrx(self,x_pts,derivative=False):
        if not(isinstance(x_pts, np.ndarray)):
            x_pts = np.asarray([x_pts])

        eval_mtrx = np.zeros((len(x_pts), len(self.lift_basis_func_prime_list)), dtype=np.float64)
        if derivative:
            for (j, phi_j) in enumerate(self.lift_basis_func_prime_list):
                eval_mtrx[:, j] = phi_j(x_pts)

            return(eval_mtrx)
        else:
            for (j, phi_j) in enumerate(self.lift_basis_func_list):
                eval_mtrx[:, j] = phi_j(x_pts)

            return(eval_mtrx)


def u_star(x, femlb, c_star_vec):
    u_0 = femlb.evaluation_matrix(x)@c_star_vec
    w = femlb.w_lift(x)
    return u_0 + w

def u_prime_star(x, femlb, c_star_vec):
    u_0 = femlb.evaluation_matrix(x, derivative=True)@c_star_vec
    #w = femlb.lift_evaluation_mtrx(x, derivative=True)@femlb.alpha_vec
    w = femlb.w_lift_prime(x)
    return u_0 + w

def diff_flux_x_star(x, femlb, c_star_vec):
    return -femlb.diff_coeff(x)*u_prime_star(x, femlb, c_star_vec)

