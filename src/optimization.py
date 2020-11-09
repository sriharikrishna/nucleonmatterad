import os
import numpy as np
from scipy.optimize import minimize

# a wrapper for the SNM call

def snm_2d_objective(x):
    # convert x into strings
    xstr = ["%.16f" % elem for elem in x]

    # write call string
    callstr = "".join(["./script_nm_snm.sh ",xstr[0]," ",xstr[1]])

    # call the call string
    os.system(callstr)

    # the output was written to out_snm.txt
    # we read in f and g from the .xt:
    file = open("out_snm.txt")
    line = file.read().replace(","," ")
    file.close()
    line = line.split()

    # write f
    f = np.float(line[0])

    # write g
    g = np.zeros((2,1))
    g[0] = np.float(line[3])
    g[1] = np.float(line[4])

    return f, g

def main():
    # initial point
    x0 = [4.0, 0.75]

    # tolerance (definition changes with solver)
    tolerance = 1e-8

    # some options to play with
    options = {'disp': True, 'maxcor': 10, 'ftol': 2.220446049250313e-09, 'gtol': 1e-05, 'eps': 1e-08, 'maxfun': 15000,
               'maxiter': 15000, 'iprint': - 1, 'maxls': 20, 'finite_diff_rel_step': None}
    res = minimize(snm_2d_objective, x0, method='L-BFGS-B', jac = True, tol = tolerance, options=options)

    #

if __name__ == '__main__':
    main()