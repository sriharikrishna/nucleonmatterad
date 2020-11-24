import numpy as np
import os
import sys

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
    output_file_2 = "out_tap_all_nucmat_pnm"
    temp_file = "temp"

    for j in range(len(xstr)):
        #output_file = output_file.join(["_", absxstr[j]])
        output_file = output_file + "_" + absxstr[j]
        output_file_2 = output_file_2 + "_" + absxstr[j]
        temp_file = temp_file + "_" + absxstr[j]
    output_file = output_file + ".txt" #output_file.join([".txt"])
    temp_file = temp_file + ".dat"

    # we read in f and g from the .txt:

    file = open(output_file)
    line = file.read().replace(","," ")
    file.close()
    line = line.split()

    # write f
    f = np.float(line[0])
    g = np.zeros(4)
    g[0] = np.float(line[5])
    g[1] = np.float(line[6])
    g[2] = np.float(line[7])
    g[3] = np.float(line[8])

    # SAFETY
    if np.isnan(f) or np.isinf(f):
        f = 1e9
        g = 1e9*np.ones(4)

    remove_string = "rm " + output_file
    os.system(remove_string)
    remove_string = "rm " + output_file_2
    os.system(remove_string)
    remove_string = "rm " + temp_file
    os.system(remove_string)

    return f, g

def main():
    args = sys.argv
    i = int(args[1])

    center = np.array([1.7637,0.5027,1.6206,4.1998])
    width = 0.05
    num_entries = 15
    first_coord_array = np.linspace(center[1]-width,center[1]+width,num_entries)
    second_coord_array = np.linspace(center[2]-width,center[2]+width,num_entries)

    f_values = np.zeros(num_entries)
    g_values = np.zeros((num_entries,4))

    # fix the first and last coordinates
    for j in range(num_entries):
        f,g = pnm_4d_objective(np.array([center[0],first_coord_array[i],second_coord_array[j],center[3]]))
        f_values[j] = f
        g_values[j,:] = g.transpose()

    filename = "contours_row_" + str(i) + ".npz"
    np.savez(filename, f_values=f_values,  g_values=g_values)

if __name__ == '__main__':
    main()