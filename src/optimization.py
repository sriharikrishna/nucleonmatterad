import os
import os.path
import numpy as np
from scipy.optimize import minimize
import nlopt

# this class wraps functions and gradient functions to write all evaluations to a file
class Funcgradmon(object):
    """ Funcgradmn: wrap f() and grad(), save all x[] f[] grad[] to plot or restart

    Example: minimize, save, restart --

    fg = Funcgradmon( func, gradfunc, verbose=1 )
        # fg(x): f(x), g(x)  for minimize( jac=True )

        # run 100 iter (if linesearch, 200-300 calls of fg()) --
    options = dict( maxiter=100 )  # ...
    min0 = minimize( fg, x0, jac=True, options=options )
    fg.savez( "0.npz", paramstr="..." )  # to plot or restart

        # restart from x[50] --
        # (won't repeat the previous path from 50
        # unless you save and restore the whole state of the optimizer)
    x0 = fg.restart( 50 )
    # change params ...
    min50 = minimize( fg, x0, jac=True, options=options )
    """

    def __init__( self, func, gradfunc, filename, verbose=1 ):
        self.func = func
        self.gradfunc = gradfunc
        self.verbose = verbose
        self.filename = filename
        self.x, self.f, self.g = [], [], []  # growing lists
        self.t = 0

    def __call__( self, x): #( self, x, grad):
        """ f, g = func(x), gradfunc(x); save them; return f, g """
        x = np.asarray_chkfinite( x )  # always
        f = self.func(x)
        g = self.gradfunc(x)
        g = np.asarray_chkfinite( g )
        self.x.append( np.copy(x) )
        self.f.append( _copy( f ))
        self.g.append( np.copy(g) )
        if self.verbose:
            print("%3d:" % self.t)
            fmt = "%-12g" if np.isscalar(f)  else "%s\t"
            print(fmt % f)
            print("x: %s" % x)  # with user's np.set_printoptions
            print("\tgrad: %s" % g)
                # better df dx dg
        # callback: plot
        self.savez(self.filename)
        self.t += 1
        #grad[:] = g
        return f, g

    def restart( self, n ):
        """ x0 = fg.restart( n )  returns x[n] to minimize( fg, x0 )
        """
        x0 = self.x[n]  # minimize from here
        del self.x[:n]
        del self.f[:n]
        del self.g[:n]
        self.t = n
        if self.verbose:
            print("Funcgradmon: restart from x[%d] %s" % (n, x0))
        return x0

    def savez( self, npzfile, **kw ):
        """ np.savez( npzfile, x= f= g= ) """
        x, f, g = map( np.array, [self.x, self.f, self.g] )
        if self.verbose:
            asum = "f: %s \nx: %s \ng: %s" % (
                _asum(f), _asum(x), _asum(g) )
            print("Funcgradmon: saving to %s: \n%s \n" % (npzfile, asum))
        np.savez( npzfile, x=x, f=f, g=g, **kw )

    def load( self, npzfile ):
        load = np.load( npzfile )
        x, f, g = load["x"], load["f"], load["g"]
        if self.verbose:
            asum = "f: %s \nx: %s \ng: %s" % (
                _asum(f), _asum(x), _asum(g) )
            print("Funcgradmon: load %s: \n%s \n" % (npzfile, asum))
        self.x = list( x )
        self.f = list( f )
        self.g = list( g )
        self.loaddict = load
        return self.restart( len(x) - 1 )


def _asum( X ):
    """ one-line array summary: "shape type min av max" """
    if not hasattr( X, "dtype" ):
        return str(X)
    return "%s %s  min av max %.3g %.3g %.3g" % (
            X.shape, X.dtype, X.min(), X.mean(), X.max() )

def _copy( x ):
    return x if x is None  or np.isscalar(x) \
        else np.copy( x )
# a wrapper for the SNM call

def snm_2d_objective(x):
    # convert x into strings
    xstr = ["%.17f" % elem for elem in x]
    absxstr = ["%.17f" % elem for elem in np.abs(x)]

    # write call string
    callstr = "".join(["./script_nm_snm.sh ",xstr[0]," ",xstr[1]])

    # call the call string
    os.system(callstr)

    # the output was written to out_snm_abs(x(i))....txt
    output_file = "out_snm"
    for j in range(len(xstr)):
        #output_file = output_file.join(["_", absxstr[j]])
        output_file = output_file + "_" + absxstr[j]
    output_file = output_file + ".txt" #output_file.join([".txt"])

    # we read in f and g from the .txt:

    file = open(output_file)
    line = file.read().replace(","," ")
    file.close()
    line = line.split()

    # write f
    f = np.float(line[0])

    # SAFETY
    if np.isnan(f) or np.isinf(f):
        f = 1e2

    return f

def snm_2d_objective_der(x):
    xstr = ["%.17f" % elem for elem in x]
    absxstr = ["%.17f" % elem for elem in np.abs(x)]

    # write call string
    callstr = "".join(["./script_nm_snm.sh ", xstr[0], " ", xstr[1]])

    ## the output was written to out_snm_abs(x(i))....txt
    output_file = "out_snm"
    for j in range(len(xstr)):
        # output_file = output_file.join(["_", absxstr[j]])
        output_file = output_file + "_" + absxstr[j]
    output_file = output_file + ".txt"  # output_file.join([".txt"])

    if not os.path.exists(output_file):
        # call the call string
        os.system(callstr)

    # we read in f and g from the .xt:
    file = open(output_file)
    line = file.read().replace(","," ")
    file.close()
    line = line.split()

    # write g
    g = np.zeros(2)
    g[0] = np.float(line[3])
    g[1] = np.float(line[4])

    remove_string = "rm " + output_file
    os.system(remove_string)

    output_file = "out_tap_all_nucmat_snm"
    for j in range(len(xstr)):
        # output_file = output_file.join(["_", absxstr[j]])
        output_file = output_file + "_" + absxstr[j]
    remove_string = "rm " + output_file
    os.system(remove_string)

    output_file = "temp"
    for j in range(len(xstr)):
        # output_file = output_file.join(["_", absxstr[j]])
        output_file = output_file + "_" + absxstr[j]
    output_file = output_file + ".dat"  # output_file.join([".txt"])
    remove_string = "rm " + output_file
    os.system(remove_string)

    return g

def snm_5d_objective(x):
    # convert x into strings
    xstr = ["%.17f" % elem for elem in x]
    absxstr = ["%.17f" % elem for elem in np.abs(x)]

    # write call string
    callstr = "".join(["./script_nm_snm_fullx.sh ",xstr[0]," ",xstr[1]," ",xstr[2]," ",xstr[3]," ",xstr[4]])

    # call the call string
    os.system(callstr)

    ## the output was written to out_snm_abs(x(i))....txt
    output_file = "out_snm"
    for j in range(len(xstr)):
        # output_file = output_file.join(["_", absxstr[j]])
        output_file = output_file + "_" + absxstr[j]
    output_file = output_file + ".txt"  # output_file.join([".txt"])

    # we read in f and g from the .xt:

    file = open(output_file)
    line = file.read().replace(",", " ")
    file.close()
    line = line.split()

    # write f
    f = np.float(line[0])

    # SAFETY
    if np.isnan(f) or np.isinf(f):
        f = 1e2

    return f

def snm_5d_objective_der(x):
    xstr = ["%.17f" % elem for elem in x]
    absxstr = ["%.17f" % elem for elem in np.abs(x)]

    # write call string
    callstr = "".join(
        ["./script_nm_snm_fullx.sh ", xstr[0], " ", xstr[1], " ", xstr[2], " ", xstr[3], " ", xstr[4]])

    ## the output was written to out_snm_abs(x(i))....txt
    output_file = "out_snm"
    for j in range(len(xstr)):
        # output_file = output_file.join(["_", absxstr[j]])
        output_file = output_file + "_" + absxstr[j]
    output_file = output_file + ".txt"  # output_file.join([".txt"])

    # we read in f and g from the .xt:
    if not os.path.exists(output_file):
        # call the call string
        os.system(callstr)

    file = open(output_file)
    line = file.read().replace(",", " ")
    file.close()
    line = line.split()

    # write g
    g = np.zeros(5)
    g[0] = np.float(line[6])
    g[1] = np.float(line[7])
    g[2] = np.float(line[8])
    g[3] = np.float(line[9])
    g[4] = np.float(line[10])

    # SAFETY:
    f = np.float(line[0])
    if np.isnan(f) or np.isinf(f):
        g = np.zeros(5)

    remove_string = "rm " + output_file
    os.system(remove_string)

    output_file = "out_tap_all_nucmat_snm"
    for j in range(len(xstr)):
        # output_file = output_file.join(["_", absxstr[j]])
        output_file = output_file + "_" + absxstr[j]
    remove_string = "rm " + output_file
    os.system(remove_string)

    output_file = "temp"
    for j in range(len(xstr)):
        # output_file = output_file.join(["_", absxstr[j]])
        output_file = output_file + "_" + absxstr[j]
    output_file = output_file + ".dat"  # output_file.join([".txt"])
    remove_string = "rm " + output_file
    os.system(remove_string)

    return g

def pnm_4d_objective(x):
    # convert x into strings
    xstr = ["%.17f" % elem for elem in x]
    absxstr = ["%.17f" % elem for elem in np.abs(x)]

    # write call string
    callstr = "".join(["./script_nm_pnm_f_only.sh ",xstr[0]," ",xstr[1]," ",xstr[2]," ",xstr[3]])

    # call the call string
    os.system(callstr)

    # the output was written to out_snm_abs(x(i))....txt
    output_file = "out_pnm"
    for j in range(len(xstr)):
        #output_file = output_file.join(["_", absxstr[j]])
        output_file = output_file + "_" + absxstr[j]
    output_file = output_file + ".txt" #output_file.join([".txt"])

    # we read in f and g from the .txt:

    file = open(output_file)
    line = file.read().replace(","," ")
    file.close()
    line = line.split()

    # write f
    f = np.float(line[0])

    # SAFETY
    if np.isnan(f) or np.isinf(f):
        f = 1e3

    return f

def pnm_4d_objective_der(x):
    xstr = ["%.17f" % elem for elem in x]
    absxstr = ["%.17f" % elem for elem in np.abs(x)]

    # write call string
    callstr = "".join(["./script_nm_pnm.sh ", xstr[0], " ", xstr[1], " ", xstr[2], " ", xstr[3]])

    ## the output was written to out_snm_abs(x(i))....txt
    output_file = "out_pnm"
    for j in range(len(xstr)):
        # output_file = output_file.join(["_", absxstr[j]])
        output_file = output_file + "_" + absxstr[j]
    output_file = output_file + ".txt"  # output_file.join([".txt"])

    if not os.path.exists(output_file):
        # call the call string
        os.system(callstr)

    # we read in f and g from the .txt:
    file = open(output_file)
    line = file.read().replace(","," ")
    file.close()
    line = line.split()

    # write g
    f = np.float(line[0])
    g = np.zeros(4)
    g[0] = np.float(line[5])
    g[1] = np.float(line[6])
    g[2] = np.float(line[7])
    g[3] = np.float(line[8])

    # SAFETY:
    if np.isnan(f) or np.isinf(f):
        #g = 1000*np.ones(4)
        g = np.zeros(4)
    remove_string = "rm " + output_file
    os.system(remove_string)

    output_file = "out_tap_all_nucmat_pnm"
    for j in range(len(xstr)):
        # output_file = output_file.join(["_", absxstr[j]])
        output_file = output_file + "_" + absxstr[j]
    remove_string = "rm " + output_file
    os.system(remove_string)

    output_file = "temp"
    for j in range(len(xstr)):
        # output_file = output_file.join(["_", absxstr[j]])
        output_file = output_file + "_" + absxstr[j]
    output_file = output_file + ".dat"  # output_file.join([".txt"])
    remove_string = "rm " + output_file
    os.system(remove_string)

    return g

def pnm_5d_objective(x):
    # convert x into strings
    xstr = ["%.17f" % elem for elem in x]
    absxstr = ["%.17f" % elem for elem in np.abs(x)]

    # write call string
    callstr = "".join(["./script_nm_pnm_fullx.sh ",xstr[0]," ",xstr[1]," ",xstr[2]," ",xstr[3]," ",xstr[4]])

    # call the call string
    os.system(callstr)

    ## the output was written to out_snm_abs(x(i))....txt
    output_file = "out_pnm"
    for j in range(len(xstr)):
        # output_file = output_file.join(["_", absxstr[j]])
        output_file = output_file + "_" + absxstr[j]
    output_file = output_file + ".txt"  # output_file.join([".txt"])

    # we read in f and g from the .xt:

    file = open(output_file)
    line = file.read().replace(",", " ")
    file.close()
    line = line.split()

    # write f
    f = np.float(line[0])

    if np.isnan(f) or np.isinf(f):
        f = 1e3

    return f

def pnm_5d_objective_der(x):
    xstr = ["%.17f" % elem for elem in x]
    absxstr = ["%.17f" % elem for elem in np.abs(x)]

    # write call string
    callstr = "".join(
        ["./script_nm_pnm_fullx.sh ", xstr[0], " ", xstr[1], " ", xstr[2], " ", xstr[3], " ", xstr[4]])

    ## the output was written to out_snm_abs(x(i))....txt
    output_file = "out_pnm"
    for j in range(len(xstr)):
        # output_file = output_file.join(["_", absxstr[j]])
        output_file = output_file + "_" + absxstr[j]
    output_file = output_file + ".txt"  # output_file.join([".txt"])

    # we read in f and g from the .xt:
    if not os.path.exists(output_file):
        # call the call string
        os.system(callstr)

    file = open(output_file)
    line = file.read().replace(",", " ")
    file.close()
    line = line.split()

    # write g
    f = np.float(line[0])
    g = np.zeros(5)
    g[0] = np.float(line[6])
    g[1] = np.float(line[7])
    g[2] = np.float(line[8])
    g[3] = np.float(line[9])
    g[4] = np.float(line[10])

    # SAFETY:
    if np.isnan(f) or np.isinf(f):
        g = np.zeros(5)

    remove_string = "rm " + output_file
    os.system(remove_string)

    output_file = "out_tap_all_nucmat_pnm"
    for j in range(len(xstr)):
        # output_file = output_file.join(["_", absxstr[j]])
        output_file = output_file + "_" + absxstr[j]
    remove_string = "rm " + output_file
    os.system(remove_string)

    output_file = "temp"
    for j in range(len(xstr)):
        # output_file = output_file.join(["_", absxstr[j]])
        output_file = output_file + "_" + absxstr[j]
    output_file = output_file + ".dat"  # output_file.join([".txt"])
    remove_string = "rm " + output_file
    os.system(remove_string)



    return g

def main():
    # args:
    # arg1: index of first starting point
    # arg2: index of second starting point
    # arg3: identifier of objective function to use

    import sys
    args = sys.argv
    start_index = int(args[1])
    end_index = int(args[2])
    problem = args[3]
    solver = args[4]

    # make clean and make!
    if problem == "snm2" and (solver == "lbfgs" or solver == "scipy_neldermead"):
        os.system("make -f MakefileTapf clean; make -f Makefile clean; make -f Makefile CASE=snm prep ; make -f MakefileTapf ALL=1 NUCMAT=1 CASE=snm")
    elif problem == "snm2" and solver == "neldermead":
        os.system("make -f MakefileTapf clean; make -f Makefile clean; make -f Makefile CASE=snm prep ; make -f MakefileTapf ALL=1 CASE=snm CUSTOM_INPUTS=1")
        #print("ok")
    elif problem == "snm5":
        os.system("make -f MakefileTapf clean; make -f Makefile clean; make -f Makefile CASE=snm prep ; make -f MakefileTapf ALL=1 FULLX=1 NUCMAT=1 CASE=snm")
    elif problem == "pnm4" and (solver == "lbfgs" or solver == "scipy_neldermead"):
        #os.system("make -f MakefileTapf clean; make -f Makefile clean; make -f Makefile CASE=pnm prep ; make -f MakefileTapf ALL=1 NUCMAT=1 CASE=pnm")
        print("ok")
    elif problem == "pnm4" and solver == "neldermead":
        os.system( "make -f MakefileTapf clean; make -f Makefile clean; make -f Makefile CASE=pnm prep ; make -f MakefileTapf ALL=1 CASE=pnm CUSTOM_INPUTS=1")
        #print("ok")
    elif problem == "pnm5":
        os.system("make -f MakefileTapf clean; make -f Makefile clean; make -f Makefile CASE=pnm prep ; make -f MakefileTapf ALL=1 FULLX=1 NUCMAT=1 CASE=pnm")
    else:
        raise Exception('invalid problem name {} supplied'.format(problem))

    # initial point
    if problem == "snm2":
        x0 = np.array([3.997, 0.796])
        dim = 2
    if problem == "snm5":
        x0 = np.array([3.997, 0.796, 0.796, 1.000, 1.000])
        dim = 5
    if problem == "pnm4":
        x0 = np.array([1.750, 0.531, 1.564, 3.747])
        dim = 4
    if problem == "pnm5":
        x0 = np.array([1.750, 0.531, 0.531, 1.564, 3.747])
        dim = 5
    # tolerance (definition changes with solver)
    tolerance = 1e-8

    # some options to play with
    options = {'disp': True, 'maxcor': dim, 'ftol': 1e-16, 'gtol': 1e-8, 'eps': 1e-08, 'maxfun': 100,
               'maxiter': 100, 'iprint': 101, 'maxls': 20}

    # read in the perturbations
    perturbation_file = open("circle.txt",'r')
    lines = perturbation_file.readlines()
    perturbation_file.close()
    perturbations = np.zeros((30,7))
    count = 0
    for line in lines:
        perturbations[count,:] = np.fromstring(line, sep = " ")
        count += 1

    for i in range(start_index,end_index):

        if i == 0:
            xi = x0
        else:
            xi = x0 + 0.1*perturbations[i-1,:dim]

        # instantiate the Funcgradmon object
        filename = problem + "_run_starting_at_x" + str(i) + solver + ".npz"
        if problem == "snm2":
            fg = Funcgradmon(snm_2d_objective, snm_2d_objective_der, filename, verbose=1)
        elif problem == "snm5":
            fg = Funcgradmon(snm_5d_objective, snm_5d_objective_der, filename, verbose=1)
        elif problem == "pnm4":
            fg = Funcgradmon(pnm_4d_objective, pnm_4d_objective_der, filename, verbose=1)
        elif problem == "pnm5":
            fg = Funcgradmon(pnm_5d_objective, pnm_5d_objective_der, filename, verbose=1)

        # remove old output and temp dat files
        #os.system("rm out* temp*")

        # call LBFGS
        # bounds: these are arbitrary, but seem more than reasonable:
        lb = -9.9999*np.ones(dim)
        ub = 9.9999*np.ones(dim)
        bounds = np.vstack((lb,ub))
        bounds = bounds.T

        if solver == "lbfgs":
            #opt = nlopt.opt(nlopt.LD_LBFGS, dim)
            #opt.set_min_objective(fg)
            #opt.set_lower_bounds(lb)
            #opt.set_upper_bounds(ub)
            ##opt.set_ftol_rel(1e-16)
            ##opt.set_ftol_abs(1e-16)
            ##opt.set_xtol_rel(1e-16)
            ##opt.set_xtol_abs(1e-16)
            #opt.set_maxeval(100)
            #opt.set_vector_storage(dim)
            #res = opt.optimize(xi)

            res = minimize(fg, xi, method='L-BFGS-B', jac = True, bounds = bounds, tol = tolerance, options=options)
        elif solver == "scipy_neldermead":
            res = minimize(fg, xi, method='Nelder-Mead', jac = True, bounds=bounds, tol=tolerance)
        elif solver == "neldermead":
            xstr = ["%.17f" % elem for elem in xi]
            absxstr = ["%.17f" % elem for elem in np.abs(xi)]
            if problem == "snm2":
                callstr = "".join(["./script_dfo_snm.sh ", xstr[0], " ", xstr[1]])
                output_file = "out_snm"
            elif problem == "pnm4":
                callstr = "".join(["./script_dfo_pnm.sh ", xstr[0], " ", xstr[1], " ", xstr[2], " ", xstr[3]])
                output_file = "out_pnm"
            os.system(callstr)

        # write the output file
        if solver == "lbfgs":
            filename = problem + "_run_starting_at_x" + str(i) + "lbfgs.npz"
            fg.savez(filename, paramstr="...")
        elif solver == "scipy_neldermead":
            filename = problem + "_run_starting_at_x" + str(i) + "scipy_neldermead.npz"
            fg.savez(filename, paramstr="...")
        elif solver == "neldermead":
            for j in range(len(xstr)):
                output_file = output_file + "_" + absxstr[j]
            output_file = output_file + ".txt"

            # open the output file, and split the data
            file = open(output_file)
            line = file.read().replace(",", " ")
            file.close()
            line = line.split()

            num_entries = len(line)

            # parse line to get arrays
            if problem == "snm2":
                num_evals = np.int(num_entries/6)
                f = np.zeros(num_evals)
                x = np.zeros((num_evals, 2))
                g = np.zeros((num_evals,2))
                for ctr in range(num_evals):
                    f[ctr] = np.float(line[6*ctr+1])
                    x[ctr,0] = np.float(line[6*ctr+2])
                    x[ctr,1] = np.float(line[6*ctr+3])
                    g[ctr,0] = np.float(line[6*ctr+4])
                    g[ctr,1] = np.float(line[6*ctr+5])
            elif problem == "pnm4":
                num_evals = np.int(num_entries/10)
                f = np.zeros(num_evals)
                x = np.zeros((num_evals, 4))
                g = np.zeros((num_evals,4))
                for ctr in range(num_evals):
                    f[ctr] = np.float(line[10*ctr+1])
                    x[ctr,0] = np.float(line[10*ctr+2])
                    x[ctr,1] = np.float(line[10*ctr+3])
                    x[ctr,2] = np.float(line[10*ctr+4])
                    x[ctr,3] = np.float(line[10*ctr+5])
                    g[ctr, 0] = np.float(line[10 * ctr + 6])
                    g[ctr, 1] = np.float(line[10 * ctr + 7])
                    g[ctr, 2] = np.float(line[10 * ctr + 8])
                    g[ctr, 3] = np.float(line[10 * ctr + 9])
            filename = problem + "_run_starting_at_x" + str(i) + "neldermead.npz"
            np.savez(filename, f=f, x=x, g=g)


if __name__ == '__main__':
    main()
