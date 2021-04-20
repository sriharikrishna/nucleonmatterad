import os
import os.path
import numpy as np
from scipy.optimize import minimize
import time

class Timer:
    def __init__(self):
        self.start = time.time()

    def __call__(self):
        return time.time() - self.start

# this class wraps functions and gradient functions to write all evaluations to a file
class Funcmon(object):

    def __init__( self, func, filename, verbose=1 ):
        self.func = func
        self.verbose = verbose
        self.filename = filename
        self.x, self.f, self.time = [], [], []  # growing lists
        self.t = 0
        self.time_passed = Timer()

    def __call__( self, x): #( self, x, grad):
        """ f, g = func(x), gradfunc(x); save them; return f, g """
        x = np.asarray_chkfinite( x )  # always
        f = self.func(x)
        secs = self.time_passed()
        self.x.append( np.copy(x) )
        self.f.append( _copy( f ))
        self.time.append(_copy(secs))
        if self.verbose:
            print("%3d:" % self.t)
            fmt = "%-12g" if np.isscalar(f)  else "%s\t"
            print(fmt % f)
            print("x: %s" % x)  # with user's np.set_printoptions
                # better df dx dg
        # callback: plot
        self.savez(self.filename)
        self.t += 1
        #grad[:] = g
        return f

    def restart( self, n ):
        """ x0 = fg.restart( n )  returns x[n] to minimize( fg, x0 )
        """
        x0 = self.x[n]  # minimize from here
        del self.x[:n]
        del self.f[:n]
        self.t = n
        if self.verbose:
            print("Funcmon: restart from x[%d] %s" % (n, x0))
        return x0

    def savez( self, npzfile, **kw ):
        """ np.savez( npzfile, x= f= g= ) """
        x, f, time = map( np.array, [self.x, self.f, self.time] )
        if self.verbose:
            asum = "f: %s \nx: %s" % (
                _asum(f), _asum(x) )
            print("Funcmon: saving to %s: \n%s \n" % (npzfile, asum))
        np.savez( npzfile, x=x, f=f, time=time, **kw )

    def load( self, npzfile ):
        load = np.load( npzfile )
        x, f = load["x"], load["f"]
        if self.verbose:
            asum = "f: %s \nx: %s " % (
                _asum(f), _asum(x) )
            print("Funcmon: load %s: \n%s \n" % (npzfile, asum))
        self.x = list( x )
        self.f = list( f )
        self.loaddict = load
        return self.restart( len(x) - 1 )

class Gradmon(object):

    def __init__( self, gradfunc, filename, verbose=1 ):
        self.gradfunc = gradfunc
        self.verbose = verbose
        self.filename = filename
        self.x,  self.g, self.time = [], [], []  # growing lists
        self.t = 0
        self.time_passed = Timer()

    def __call__( self, x): #( self, x, grad):
        """ f, g = func(x), gradfunc(x); save them; return f, g """
        x = np.asarray_chkfinite( x )  # always
        g = self.gradfunc(x)
        g = np.asarray_chkfinite( g )
        secs = self.time_passed()
        self.time.append(_copy(secs))
        self.x.append( np.copy(x) )
        self.g.append( np.copy(g) )
        if self.verbose:
            print("%3d:" % self.t)
            print("x: %s" % x)  # with user's np.set_printoptions
            print("\tgrad: %s" % g)
                # better df dx dg
        # callback: plot
        self.savez(self.filename)
        self.t += 1
        #grad[:] = g
        return g

    def restart( self, n ):
        """ x0 = fg.restart( n )  returns x[n] to minimize( fg, x0 )
        """
        x0 = self.x[n]  # minimize from here
        del self.x[:n]
        del self.g[:n]
        self.t = n
        if self.verbose:
            print("Funcgradmon: restart from x[%d] %s" % (n, x0))
        return x0

    def savez( self, npzfile, **kw ):
        """ np.savez( npzfile, x= f= g= ) """
        x, g, time = map( np.array, [self.x, self.g, self.time] )
        if self.verbose:
            asum = "f: %s \nx: %s \ng: %s" % (
                _asum(f), _asum(x), _asum(g) )
            print("Funcgradmon: saving to %s: \n%s \n" % (npzfile, asum))
        np.savez( npzfile, x=x, g=g, time=time, **kw )

    def load( self, npzfile ):
        load = np.load( npzfile )
        x, f, g = load["x"], load["f"], load["g"]
        if self.verbose:
            asum = "x: %s \ng: %s" % ( _asum(x), _asum(g) )
            print("Funcgradmon: load %s: \n%s \n" % (npzfile, asum))
        self.x = list( x )
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
    callstr = "".join(["./script_nm_snm_f_only.sh ",xstr[0]," ",xstr[1]])

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
    elif problem == "pnm4" and solver == "lbfgs":
        os.system("make -f MakefileTapf clean; make -f Makefile clean; make -f Makefile CASE=pnm prep ; make -f MakefileTapf ALL=1 NUCMAT=1 CASE=pnm")
        print("ok")
    elif problem == "pnm4" and solver == "scipy_neldermead"):
        os.system("mkdir -p pnm/obj/; make clean; make -f MakefileTapf clean; make prep CASE=pnm; make CASE=pnm CUSTOM_INPUTS=1 NUCMAT=1")
        #print("ok")
    elif problem == "pnm4" and solver == "neldermead":
        os.system( "make -f MakefileTapf clean; make -f Makefile clean; make -f Makefile CASE=pnm prep ; make -f MakefileTapf ALL=1 CASE=pnm CUSTOM_INPUTS=1")
        #print("ok")
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

        # instantiate the Func/Gradmon object
        f_filename = problem + "f_timing_run_starting_at_x" + str(i) + solver + ".npz"
        g_filename = problem + "g_timing_run_starting_at_x" + str(i) + solver + ".npz"
        if problem == "snm2":
            f = Funcmon(snm_2d_objective,f_filename,verbose=1)
            g = Gradmon(snm_2d_objective_der,g_filename,verbose=1)
        elif problem == "pnm4":
            f = Funcmon(pnm_4d_objective, f_filename, verbose=1)
            g = Gradmon(pnm_4d_objective_der,g_filename, verbose=1)

        # remove old output and temp dat files
        #os.system("rm out* temp*")

        # call LBFGS
        # bounds: these are arbitrary, but seem more than reasonable:
        lb = -9.9999*np.ones(dim)
        ub = 9.9999*np.ones(dim)
        bounds = np.vstack((lb,ub))
        bounds = bounds.T

        if solver == "lbfgs":
            res = minimize(f, xi, method='L-BFGS-B', jac = g, bounds = bounds, tol = tolerance, options=options)
        elif solver == "scipy_neldermead":
            res = minimize(f, xi, method='Nelder-Mead', bounds=bounds, tol=tolerance)
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

if __name__ == '__main__':
    main()
