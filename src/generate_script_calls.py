import os
import sys
import argparse
from multiprocessing import Process

def scriptgen(ndim, case, dor, ast, als, atn,
              bst, btn, bls, rho, solver="lbfgs"):
    cmdlist = []
    count = 0
    for ls in range (10,33,2):
        for ll in range (10,ls+1,2):
            for lc in range (10,ll+1,2):
                count+=1
                if ndim == 2:
                    cmd = "python3 optimization.py {} {} {} {} {} {} {} {} > opt_out_{}_{}_{}_{}_{}_{}".format(dor, ast, case, solver, rho, lc, ls, ll, dor, ast, rho, lc, ls, ll)
                if ndim == 4:
                    cmd = "python3 optimization.py {} {} {} {} {} {} {} {} {} {} > opt_out_{}_{}_{}_{}_{}_{}_{}_{}".format(dor, ast, bst, btn, case, solver, rho, lc, ls, ll, dor, ast, bst, btn, rho, lc, ls, ll)
                if ndim == 7:
                    cmd = "python3 optimization.py {} {} {} {} {} {} {} {} {} {} {} {} {} > opt_out_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}".format(dor, ast, atn, als, bst, btn, bls, case, solver, rho, lc, ls, ll, dor, ast, atn, als, bst, btn, bls, rho, lc, ls, ll)
                cmdlist.append(cmd)
    print(count)
    return cmdlist

def call_script(cmd):
    print(cmd)
    #os.system(cmd)

def runbatch(batch=0, cmdlist=[]):
    N = min(os.cpu_count(), 48)
    print("batch = ", batch, "startingat = ",(batch-1)*N)
    for count in range((batch-1)*N, batch*N):
      cmd = cmdlist[count-1]
      p = Process(target=call_script, args=(cmd,))
      p.start()
      os.system("taskset -p -c %d %d" % ((count % os.cpu_count()), p.pid))

if __name__ == "__main__":
    print("os.cpu_count()", os.cpu_count(), (364/os.cpu_count()))
    # parse arguments
    parser = argparse.ArgumentParser(
        description="Generate scripts of SNM/PNM optimization calls.",
        formatter_class=argparse.MetavarTypeHelpFormatter,
    )
    parser.add_argument("--batch", type=int, required=True, choices=range(1,int(364/os.cpu_count()+1)), help='''There are 364 combinations of lc,ls,lt.
                        They can run batches of number of CPUs available in your system. 
                        This argument requires selects 364/os.cpu_count() batch you want to run''')
    parser.add_argument("--ndim", required=True, type=int, choices=[2,4,7])
    parser.add_argument("--case", type=str, default="snm2", choices=["snm2", "pnm4"])
    parser.add_argument("--dor", type=float, help = "default: %(default)s", default=5.528)
    parser.add_argument("--ast", type=float, help = "default: %(default)s", default=1.557)
    parser.add_argument("--atn", type=float, help = "default: %(default)s", default=0.569)
    parser.add_argument("--als", type=float, help = "default: %(default)s", default=3.365)
    parser.add_argument("--bst", type=float, help = "default: %(default)s", default=1.527)
    parser.add_argument("--btn", type=float, help = "default: %(default)s", default=1.276)
    parser.add_argument("--bls", type=float, help = "default: %(default)s", default=0.116)
    parser.add_argument("--rho", type=float, help = "default: %(default)s", default=0.7)
    args = parser.parse_args(sys.argv[1:])
    cmdlist = scriptgen(args.ndim, args.case, args.dor, args.ast, args.als, args.atn,
              args.bst, args.btn, args.bls, args.rho)
    runbatch(args.batch, cmdlist)
