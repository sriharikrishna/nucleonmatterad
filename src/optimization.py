import os
import os.path
import numpy as np
from scipy.optimize import minimize

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

    def __init__( self, func, gradfunc, verbose=1 ):
        self.func = func
        self.gradfunc = gradfunc
        self.verbose = verbose
        self.x, self.f, self.g = [], [], []  # growing lists
        self.t = 0

    def __call__( self, x ):
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
        self.t += 1
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

    return g

def snm_7d_objective(x):
    # convert x into strings
    xstr = ["%.17f" % elem for elem in x]
    absxstr = ["%.17f" % elem for elem in np.abs(x)]

    # write call string
    callstr = "".join(["./script_nm_snm_fullx.sh ",xstr[0]," ",xstr[1]," ",xstr[2]," ",xstr[3]," ",xstr[4]," ",xstr[5]," ",xstr[6]])

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

def snm_7d_objective_der(x):
    xstr = ["%.17f" % elem for elem in x]
    absxstr = ["%.17f" % elem for elem in np.abs(x)]

    # write call string
    callstr = "".join(
        ["./script_nm_snm_fullx.sh ", xstr[0], " ", xstr[1], " ", xstr[2], " ", xstr[3], " ", xstr[4], " ", xstr[5],
         " ", xstr[6]])

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
    g = np.zeros(7)
    g[0] = np.float(line[8])
    g[1] = np.float(line[9])
    g[2] = np.float(line[10])
    g[3] = np.float(line[11])
    g[4] = np.float(line[12])
    g[5] = np.float(line[13])
    g[6] = np.float(line[14])

    # SAFETY:
    if np.isnan(f) or np.isinf(f):
        g = np.zeros(4)

    remove_string = "rm "
    remove_string = remove_string.join([output_file])
    os.system(remove_string)

    return g

def pnm_4d_objective(x):
    # convert x into strings
    xstr = ["%.17f" % elem for elem in x]
    absxstr = ["%.17f" % elem for elem in np.abs(x)]

    # write call string
    callstr = "".join(["./script_nm_pnm.sh ",xstr[0]," ",xstr[1]," ",xstr[2]," ",xstr[3]])

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

    # we read in f and g from the .xt:
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
        g = np.zeros(4)

    remove_string = "rm " + output_file
    os.system(remove_string)

    return g

def pnm_7d_objective(x):
    # convert x into strings
    xstr = ["%.17f" % elem for elem in x]
    absxstr = ["%.17f" % elem for elem in np.abs(x)]

    # write call string
    callstr = "".join(["./script_nm_pnm_fullx.sh ",xstr[0]," ",xstr[1]," ",xstr[2]," ",xstr[3]," ",xstr[4]," ",xstr[5]," ",xstr[6]])

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

def pnm_7d_objective_der(x):
    xstr = ["%.17f" % elem for elem in x]
    absxstr = ["%.17f" % elem for elem in np.abs(x)]

    # write call string
    callstr = "".join(
        ["./script_nm_pnm_fullx.sh ", xstr[0], " ", xstr[1], " ", xstr[2], " ", xstr[3], " ", xstr[4], " ", xstr[5],
         " ", xstr[6]])

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
    g = np.zeros(7)
    g[0] = np.float(line[8])
    g[1] = np.float(line[9])
    g[2] = np.float(line[10])
    g[3] = np.float(line[11])
    g[4] = np.float(line[12])
    g[5] = np.float(line[13])
    g[6] = np.float(line[14])

    # SAFETY:
    if np.isnan(f) or np.isinf(f):
        g = np.zeros(7)

    remove_string = "rm "
    remove_string = remove_string.join([output_file])
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
    if problem == "snm2":
        os.system("make -f MakefileTapf clean; make -f Makefile clean; make -f Makefile CASE=snm prep ; make -f MakefileTapf ALL=1 NUCMAT=1 CASE=snm")
    elif problem == "snm7":
        os.system("make -f MakefileTapf clean; make -f Makefile clean; make -f Makefile CASE=snm prep ; make -f MakefileTapf ALL=1 FULLX=1 NUCMAT=1 CASE=snm")
    elif problem == "pnm4":
        os.system("make -f MakefileTapf clean; make -f Makefile clean; make -f Makefile CASE=pnm prep ; make -f MakefileTapf ALL=1 NUCMAT=1 CASE=pnm")
        #print("ok")
    elif problem == "pnm7":
        os.system("make -f MakefileTapf clean; make -f Makefile clean; make -f Makefile CASE=pnm prep ; make -f MakefileTapf ALL=1 FULLX=1 NUCMAT=1 CASE=pnm")
    else:
        raise Exception('invalid problem name {} supplied'.format(problem))

    # initial point
    if problem == "snm2":
        x0 = np.array([3.997, 0.796])
        dim = 2
    if problem == "snm7":
        x0 = np.array([3.997, 0.796, 0.796, 0.796, 1.000, 1.000, 0.000])
        dim = 7
    if problem == "pnm4":
        x0 = np.array([1.750, 0.531, 1.564, 3.747])
        dim = 4
    if problem == "pnm7":
        x0 = np.array([1.750, 0.531, 0.531, 0.531, 1.564, 3.747, 1.564])
        dim = 7
    # tolerance (definition changes with solver)
    tolerance = 1e-8

    # some options to play with
    options = {'disp': True, 'maxcor': 10, 'ftol': 1e-16, 'gtol': 1e-8, 'eps': 1e-08, 'maxfun': 1000,
               'maxiter': 1000, 'iprint': 101, 'maxls': 20}

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
            xi = x0 + perturbations[i-1,:dim]

        # instantiate the Funcgradmon object
        if problem == "snm2":
            fg = Funcgradmon(snm_2d_objective, snm_2d_objective_der, verbose=1)
        elif problem == "snm7":
            fg = Funcgradmon(snm_7d_objective, snm_7d_objective_der, verbose=1)
        elif problem == "pnm4":
            fg = Funcgradmon(pnm_4d_objective, pnm_4d_objective_der, verbose=1)
        elif problem == "pnm7":
            fg = Funcgradmon(pnm_7d_objective, pnm_7d_objective_der, verbose=1)

        # remove old output and temp dat files
        #os.system("rm out* temp*")

        # call LBFGS
        # bounds: these are arbitrary, but seem more than reasonable:
        lb = -9.9999*np.ones(dim)
        ub = 9.9999*np.ones(dim)
        bounds = np.vstack((lb,ub))
        bounds = bounds.T

        if solver == "lbfgs":
            res = minimize(fg, xi, method='L-BFGS-B', jac = True, bounds = bounds, tol = tolerance, options=options)
        elif solver == "neldermead":
            res = minimize(fg, xi, method='Nelder-Mead', jac=True, bounds=bounds, tol=tolerance, options=options)

        # write the output file
        if solver == "lbfgs":
            filename = problem + "_run_starting_at_x" + str(i) + ".npz"
        elif solver == "neldermead":
            filename = problem + "_run_starting_at_x" + str(i) + "neldermead.npz"
        fg.savez(filename, paramstr="...")

if __name__ == '__main__':
    main()