import os
import os.path
import numpy as np
from scipy.optimize import minimize
import nlopt
import glob

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

    def __init__( self, func, filename, computing_grads=True, verbose=1 ):
        self.func = func
        self.verbose = verbose
        self.filename = filename
        self.x, self.f = [], []
        self.computing_grads = computing_grads
        if computing_grads:
            self.g = []
        self.t = 0

    def __call__( self, x):
        """ f, g = func(x), gradfunc(x); save them; return f, g """
        x = np.asarray_chkfinite( x )  # always
        if self.computing_grads:
            f,g = self.func(x)
            g = np.asarray_chkfinite( g )
            self.g.append( np.copy(g) )
        else:
            f = self.func(x)
        self.x.append( np.copy(x) )
        self.f.append( _copy( f ))
        if self.verbose:
            print("%3d:" % self.t)
            fmt = "%-12g" if np.isscalar(f)  else "%s\t"
            print(fmt % f)
            print("x: %s" % x)  # with user's np.set_printoptions
            if self.computing_grads:
                print("\tgrad: %s" % g)
        self.savez(self.filename)
        self.t += 1
        if self.computing_grads:
            return f, g
        else:
            return f

    def restart( self, n ):
        """ x0 = fg.restart( n )  returns x[n] to minimize( fg, x0 )
        """
        x0 = self.x[n]  # minimize from here
        del self.x[:n]
        del self.f[:n]
        if self.computing_grads:
            del self.g[:n]
        self.t = n
        if self.verbose:
            print("Funcgradmon: restart from x[%d] %s" % (n, x0))
        return x0

    def savez( self, npzfile, **kw ):
        """ np.savez( npzfile, x= f= g= ) """
        if self.computing_grads:
            x, f, g = map( np.array, [self.x, self.f, self.g] )
            np.savez( npzfile, x=x, f=f, g=g, **kw )
        else:
            x, f = map(np.array, [self.x, self.f] )
            np.savez( npzfile, x=x, f=f, **kw )

    def load( self, npzfile ):
        load = np.load( npzfile )
        if self.computing_grads:
            x, f, g = load["x"], load["f"], load["g"]
            self.g = list( g )
        else:
            x, f = load["x"], load["f"]
        self.x = list( x )
        self.f = list( f )
        self.loaddict = load
        return self.restart( len(x) - 1 )

def _copy( x ):
    return x if x is None  or np.isscalar(x) \
        else np.copy( x )

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

    # SAFETY
    f = np.float(line[0])
    if np.isinf(f) or np.isnan(f):
        g = np.zeros(2)

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

def pnm_4d_objective(x,rho,lc,ls,lt):
    # convert x into strings
    xstr = ["%.4f" % elem for elem in x]
    absxstr = ["%.4f" % elem for elem in np.abs(x)]
    rhostr = "%.4f" % rho

    # write call string
    callstr = "".join(["./script_nm_pnm_f_only_4.sh ",xstr[0]," ",xstr[1]," ",xstr[2]," ",xstr[3], " ", rhostr, " ", str(lc), " ", str(ls), " ", str(lt)])

    # call the call string
    os.system(callstr)

    # the output was written to out_snm_abs(x(i))....txt
    output_file = "out_pnm"
    for j in range(len(xstr)):
        #output_file = output_file.join(["_", absxstr[j]])
        output_file = output_file + "_" + absxstr[j]
    output_file = output_file + "_" + rhostr + "_" + str(lc) + "_" + str(ls) + "_" + str(lt) + ".txt"

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
    print(remove_string)
    os.system(remove_string)

    output_file = "out_tap_all_nucmat_pnm"
    for j in range(len(xstr)):
        output_file = output_file + "_" + absxstr[j]
    output_file = output_file + "_" + rhostr + "_" + str(lc) + "_" + str(ls) + "_" + str(lt)
    remove_string = "rm " + output_file
    print(remove_string)
    os.system(remove_string)

    output_file = "temp"
    for j in range(len(xstr)):
        output_file = output_file + "_" + absxstr[j]
    output_file = output_file + "_" + rhostr + "_" + str(lc) + "_" + str(ls) + "_" + str(lt) + ".dat"
    remove_string = "rm " + output_file
    print(remove_string)
    os.system(remove_string)

    return f

def pnm_4d_objective_der(x,rho,lc,ls,lt):
    xstr = [str(elem) for elem in x]
    rhostr = "%.3f" % rho

    # write call string
    callstr = "".join(["./script_nm_pnm_4.sh ", xstr[0], " ", xstr[1], " ", xstr[2], " ", xstr[3], " ", rhostr, " ", str(lc), " ", str(ls), " ", str(lt)])

    os.system(callstr)

    # we read in f and g from the .txt:
    searchstr = "out_pnm*" + rhostr + "_" + str(lc) + "_" + str(ls) + "_" + str(lt) + ".txt"
    output_file = glob.glob(searchstr)
    pnm_file = output_file[0]
    file = open(pnm_file)
    line = file.read().replace(","," ")
    file.close()
    line = line.split()

    # write g
    f = np.float(line[0])
    g = np.zeros(4)
    g[0] = np.float(line[6])
    g[1] = np.float(line[7])
    g[2] = np.float(line[8])
    g[3] = np.float(line[9])

    # SAFETY:
    if np.isnan(f):
        f = 1e3
        g = np.zeros(4)
    
    # CLEANUP:
    remove_string = "rm " + pnm_file
    os.system(remove_string)

    searchstr = "out_tap*" + rhostr + "_" + str(lc) + "_" + str(ls) + "_" + str(lt)
    output_file = glob.glob(searchstr)
    remove_string = "rm " + output_file[0]
    os.system(remove_string)

    searchstr = "temp*" + rhostr + "_" + str(lc) + "_" + str(ls) + "_" + str(lt) + ".dat"
    output_file = glob.glob(searchstr)
    remove_string = "rm " + output_file[0]
    os.system(remove_string)

    return f,g

def main():
    # args:
    # arg1: dor
    # arg2: alpha
    # arg3: betas
    # arg4: betat
    # arg5: problem (e.g., "pnm4")
    # arg6: solver (e.g., "lbfgs")
    # arg7: rho
    # arg8: lc
    # arg9: ls
    # arg10: lt

    import sys
    args = sys.argv
    dor = args[1]
    alpha = args[2]
    betas = args[3]
    betat = args[4]
    problem = args[5]
    solver = args[6]
    rho = float(args[7])
    lc = int(args[8])
    ls = int(args[9])
    lt = int(args[10])

    # initial point
    # WARNING, CURRENTLY HARDCODED TO READ FIRST FOUR ARGUMENTS AS INITIAL POINT FOR PNM, MUST FIX LATER
    if problem == "snm2":
        x0 = np.array([3.997, 0.796])
        dim = 2
    if problem == "pnm4":
        x0 = np.array([float(dor), float(alpha), float(betas), float(betat)])
        dim = 4 
    # tolerance (definition changes with solver)
    tolerance = 1e-8

    # some options to play with
    options = {'disp': True, 'maxcor': dim, 'ftol': 1e-16, 'gtol': 1e-8, 'eps': 1e-08, 'maxfun': 100,
               'maxiter': 100, 'iprint': 101, 'maxls': 20}

    # read in the perturbations
    #perturbation_file = open("circle.txt",'r')
    #lines = perturbation_file.readlines()
    #perturbation_file.close()
    #perturbations = np.zeros((30,7))
    #count = 0
    #for line in lines:
    #    perturbations[count,:] = np.fromstring(line, sep = " ")
    #    count += 1

    start_index = 0
    end_index = 1
    for i in range(start_index,end_index):

        if i == 0:
            xi = x0
        else:
            xi = x0 + 0.1*perturbations[i-1,:dim]

        # instantiate the Funcgradmon object
        filename = problem + "_run_starting_at_x" + str(i) + solver + "_rho=" + str(rho) + "_lc_" + str(lc) + "_ls_" + str(ls) + "_lt_" + str(lt) + ".npz"

        if problem == "snm2":
            raise Exception('snm2 is currently not implemented, pending Krishna changes to allow for rho, lc, ls, lt inputs')
            #fg = Funcgradmon(snm_2d_objective, snm_2d_objective_der, filename, verbose=1)
        elif problem == "pnm4":
            if solver == "lbfgs":
                objective = lambda x: pnm_4d_objective_der(x,rho,lc,ls,lt)
                fg = Funcgradmon(objective, filename, computing_grads = True,verbose=1)
            else:
                objective = lambda x: pnm_4d_objective(x,rho,lc,ls,lt)
                fg = Funcgradmon(objective, filename, computing_grads = False,verbose=1)

        # call LBFGS
        # bounds: these are arbitrary, but seem more than reasonable:
        lb = -9.999*np.ones(dim)
        ub = 9.999*np.ones(dim)
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
            res = minimize(fg, xi, method='Nelder-Mead', tol=tolerance)
        elif solver == "fd_lbfgs":
            res = minimize(fg, xi, method='L-BFGS-B', tol=tolerance)
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
            filename = problem + "_rho=" + str(rho) + "_lc=" + str(lc) + "_ls=" + str(ls) + "_lt=" + str(lt) + "_run_starting_at_x" + str(i) + "lbfgs.npz"
            fg.savez(filename, paramstr="...")
        elif solver == "scipy_neldermead":
            filename = problem + "_run_starting_at_x" + str(i) + "scipy_neldermead.npz"
            fg.savez(filename, paramstr="...")
        elif solver == "fd_lbfgs":
            filename = problem + "_run_starting_at_x" + str(i) + "fd_lbfgs.npz"
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
