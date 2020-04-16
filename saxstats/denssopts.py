#!/usr/bin/env python

from ._version import __version__
import os, argparse
import imp
try:
    imp.find_module('matplotlib')
    matplotlib_found = True
    import matplotlib.pyplot as plt
    from  matplotlib.colors import colorConverter as cc
    import matplotlib.gridspec as gridspec
except ImportError:
    matplotlib_found = False

def param_file_handler(line):
    '''
    Small function that helps with handling of the parameter file input JAS
    '''
    line = line.strip().strip(";")
    file, file_args = line.split(':')
    file_args = file_args.strip().split(';')
    file_arg_dict = {val.split("==")[0]:val.split("==")[1] for val in file_args}
    for ele in file_arg_dict.keys():
        try:
            if ele == 'minimum_density' or ele == 'maximum_density' or ele == 'buffer_scattering_length_densities': #should be float
                file_arg_dict[ele] = float(file_arg_dict[ele])
            elif file_arg_dict[ele][0] == "[": #should be a list
                temp = file_arg_dict[ele].strip("[").strip("]").split(",")
                file_arg_dict[ele] = [int(num.strip()) for num in temp]
            elif file_arg_dict[ele][0] == 't' or file_arg_dict[ele][0] == 'T': #value is true
                file_arg_dict[ele] = True
            elif file_arg_dict[ele][0] == 'f' or file_arg_dict[ele][0] == 'F': #value is true
                file_arg_dict[ele] = False
            else: #should be integer
                file_arg_dict[ele] = int(file_arg_dict[ele])
        except ValueError: #if the value isn't an integer/float/bool/list of integers
            continue
    return file, file_arg_dict

def parse_arguments(parser,gnomdmax=None):

    parser.add_argument("--version", action="version",version="%(prog)s v{version}".format(version=__version__))
    parser.add_argument("-f", "--file", type=str, help="SAXS data file for input (either .dat or .out)")
    parser.add_argument("-fm", "--filemultiple", default = None, nargs = '+', type=str, help="Input multiple SAXS/SANS scattering profiles for complex analysis (.dat or .out). Include spaces between each file.")
    parser.add_argument("-paramf", "--parameter_file", default = None, type=argparse.FileType('r'), help="Specify input parameters for each file per line and keep the same ordering for the files as they are entered on the command line.")
    parser.add_argument("-avg_steps", "--avg_steps", default = 1, type=int, help="Number of steps before all Neutron Contrast data is averaged together. Only applies to -fm or --filemultiple")
    parser.add_argument("-bsld", "--buffer_scattering_length_densities", default = None, nargs="+", type=float, help="Input buffer scattering length densities that correspond to files entered with -fm. Add spaces between each number." )
    parser.add_argument("-avg_w", "--average_weights", default = None, nargs="+", type=float, help="Argument where you can input the weights used for density averaging parallel to the --filemultiple input.")
    parser.add_argument("-cboost_start", "--contrast_boosting_start", default = None, type=int, help="Add a step number when the contrast boosting begins for your DENSS Multiple reconstruction.")
    parser.add_argument("-pct_noise", "--reference_percent_noise", default = 0.1, type=float, help="Percent of values in the reference .mrc file to be altered randomly")
    parser.add_argument("-noise_on", "--noise_on", dest="noisy", action="store_true", help="Turns on the feature to add noise to a given reference .mrc file")
    parser.add_argument("-u", "--units", default="a", type=str, help="Angular units (\"a\" [1/angstrom] or \"nm\" [1/nanometer]; default=\"a\")")
    parser.add_argument("-d", "--dmax", default=None, type=float, help="Estimated maximum dimension")
    parser.add_argument("-v", "--voxel", default=None, type=float, help="Set desired voxel size, setting resolution of map")
    parser.add_argument("-os","--oversampling", default=3., type=float, help="Sampling ratio")
    parser.add_argument("-n", "--nsamples", default=None, type=int, help="Number of samples, i.e. grid points, along a single dimension. (Sets voxel size, overridden by --voxel. Best optimization with n=power of 2. Default=64)")
    parser.add_argument("--ne", default=10000, type=float, help="Number of electrons in object")
    parser.add_argument("-s", "--steps", type=int, default=None, help="Maximum number of steps (iterations)")
    parser.add_argument("-ncs", "--ncs", default=0, type=int, help="Rotational symmetry")
    parser.add_argument("-ncs_steps","--ncs_steps", default=[3000,5000,7000,9000], type=int, nargs='+', help="List of steps for applying NCS averaging (default=3000,5000,7000,9000)")
    parser.add_argument("-ncs_axis", "--ncs_axis", default=1, type=int, help="Rotational symmetry axis (options: 1, 2, or 3 corresponding to xyz principal axes)")
    parser.add_argument("-o", "--output", default=None, help="Output map filename")
    parser.add_argument("-m", "--mode", default="SLOW", type=str, help="Mode. F(AST) sets default options to run quickly for simple particle shapes. S(LOW) useful for more complex molecules. M(EMBRANE) mode allows for negative contrast. (default SLOW)")
    parser.add_argument("--seed", default=None, help="Random seed to initialize the map")
    parser.add_argument("-ld_on","--limit_dmax_on", dest="limit_dmax", action="store_true", help="Limit electron density to sphere of radius 0.6*Dmax from center of object.")
    parser.add_argument("-ld_off","--limit_dmax_off", dest="limit_dmax", action="store_false", help="Do not limit electron density to sphere of radius 0.6*Dmax from center of object. (default)")
    parser.add_argument("-ld_steps","--limit_dmax_steps", default=None, type=int, nargs='+', help="List of steps for limiting density to sphere of Dmax (default=500)")
    parser.add_argument("-rc_on", "--recenter_on", dest="recenter", action="store_true", help="Recenter electron density when updating support. (default)")
    parser.add_argument("-rc_off", "--recenter_off", dest="recenter", action="store_false", help="Do not recenter electron density when updating support.")
    parser.add_argument("-rc_steps", "--recenter_steps", default=None, type=int, nargs='+', help="List of steps to recenter electron density.")
    parser.add_argument("-rc_mode", "--recenter_mode", default="com", type=str, help="Recenter based on either center of mass (com, default) or maximum density value (max)")
    parser.add_argument("-p_on","--positivity_on", dest="positivity", action="store_true", help="Enforce positivity restraint inside support. (default)")
    parser.add_argument("-p_off","--positivity_off", dest="positivity", action="store_false", help="Do not enforce positivity restraint inside support.")
    parser.add_argument("-neg_on","--negativity_on", dest="negativity", action="store_true", help="Enforce negativity restraint inside support. (default)")
    parser.add_argument("-neg_off","--negativity_off", dest="negativity", action="store_false", help="Do not enforce negativity restraint inside support.")
    parser.add_argument("-fld_on","--flatten_low_density_on", dest="flatten_low_density", action="store_true", help="Density values near zero (.01 e-/A3) will be set to zero. (default)")
    parser.add_argument("-fld_off","--flatten_low_density_off", dest="flatten_low_density", action="store_false", help="Density values near zero (.01 e-/A3) will be not set to zero.")
    parser.add_argument("-min","--minimum_density", default=None, type=float, help="Minimum density value in e-/angstrom^3 (must also set --ne to be meaningful)")
    parser.add_argument("-max","--maximum_density", default=None, type=float, help="Maximum density value in e-/angstrom^3 (must also set --ne to be meaningful)")
    parser.add_argument("-rho", "--rho_start", default=None, type=str, help="Starting electron density map filename (for use in denss.refine.py only)")
    parser.add_argument("-e_on","--extrapolate_on", dest="extrapolate", action="store_true", help="Extrapolate data by Porod law to high resolution limit of voxels. (default)")
    parser.add_argument("-e_off","--extrapolate_off", dest="extrapolate", action="store_false", help="Do not extrapolate data by Porod law to high resolution limit of voxels.")
    parser.add_argument("-sw_on","--shrinkwrap_on", dest="shrinkwrap", action="store_true", help="Turn shrinkwrap on (default)")
    parser.add_argument("-sw_off","--shrinkwrap_off", dest="shrinkwrap", action="store_false", help="Turn shrinkwrap off")
    parser.add_argument("-sw_start","--shrinkwrap_sigma_start", default=3, type=float, help="Starting sigma for Gaussian blurring, in voxels")
    parser.add_argument("-sw_end","--shrinkwrap_sigma_end", default=1.5, type=float, help="Ending sigma for Gaussian blurring, in voxels")
    parser.add_argument("-sw_decay","--shrinkwrap_sigma_decay", default=0.99, type=float, help="Rate of decay of sigma, fraction")
    parser.add_argument("-sw_threshold","--shrinkwrap_threshold_fraction", default=0.20, type=float, help="Minimum threshold defining support, in fraction of maximum density")
    parser.add_argument("-sw_iter","--shrinkwrap_iter", default=20, type=int, help="Number of iterations between updating support with shrinkwrap")
    parser.add_argument("-sw_minstep","--shrinkwrap_minstep", default=None, type=int, help="First step to begin shrinkwrap")
    parser.add_argument("-ec_on","--enforce_connectivity_on", dest="enforce_connectivity", action="store_true", help="Enforce connectivity of support, i.e. remove extra blobs (default)")
    parser.add_argument("-ec_off","--enforce_connectivity_off", dest="enforce_connectivity", action="store_false", help="Do not enforce connectivity of support")
    parser.add_argument("-ec_steps","--enforce_connectivity_steps", default=None, type=int, nargs='+', help="List of steps to enforce connectivity")
    parser.add_argument("-cef", "--chi_end_fraction", default=0.001, type=float, help="Convergence criterion. Minimum threshold of chi2 std dev, as a fraction of the median chi2 of last 100 steps.")
    parser.add_argument("--write_xplor_format", default=False, action="store_true", help="Write out XPLOR map format (default only write MRC format).")
    parser.add_argument("--write_freq", default=100, type=int, help="How often to write out current density map (in steps, default 100).")
    parser.add_argument("--cutout_on", dest="cutout", action="store_true", help="When writing final map, cut out the particle to make smaller files.")
    parser.add_argument("--cutout_off", dest="cutout", action="store_false", help="When writing final map, do not cut out the particle to make smaller files (default).")
    parser.add_argument("--plot_on", dest="plot", action="store_true", help="Create simple plots of results (requires Matplotlib, default if module exists).")
    parser.add_argument("--plot_off", dest="plot", action="store_false", help="Do not create simple plots of results. (Default if Matplotlib does not exist)")
    parser.add_argument("-q", "--quiet", action="store_true", help="Do not display running statistics. (default False)")
    parser.add_argument("--force_run", action="store_true", help="Force denss to run even if q=0 does not exist. (default False)")
    parser.set_defaults(limit_dmax=False)
    parser.set_defaults(shrinkwrap=True)
    parser.set_defaults(recenter=True)
    parser.set_defaults(positivity=True)
    parser.set_defaults(negativity=False)
    parser.set_defaults(flatten_low_density=False)
    parser.set_defaults(extrapolate=True)
    parser.set_defaults(enforce_connectivity=True)
    parser.set_defaults(cutout=False)
    parser.set_defaults(quiet = False)
    parser.set_defaults(force_run = False)
    if matplotlib_found:
        parser.set_defaults(plot=True)
    else:
        parser.set_defaults(plot=False)
    args = parser.parse_args()

    if args.output is None:
        basename, ext = os.path.splitext(args.file)
        args.output = basename
    else:
        args.output = args.output

    #A bug appears to be present when disabling the porod extrapolation.
    #for now that option will be disabled until I come up with a fix
    if args.extrapolate is False:
        print ("There is currently a bug when disabling the Porod "
               "extrapolation (the -e_off option). \n For now, extrapolation "
               "has been re-enabled until a bug fix is released. ")
        args.extrapolate = True

    #for FAST or SLOW modes, set some default values for a few options
    if args.mode[0].upper() == "F":
        args.mode = "FAST"
        nsamples = 32
        shrinkwrap_minstep = 1000
        enforce_connectivity_steps = [2000]
        recenter_steps = list(range(501,2502,500))
    elif args.mode[0].upper() == "S":
        args.mode = "SLOW"
        nsamples = 64
        shrinkwrap_minstep = 5000
        enforce_connectivity_steps = [6000]
        recenter_steps = list(range(501,8002,500))
    elif args.mode[0].upper() == "M":
        args.mode = "MEMBRANE"
        nsamples = 64
        args.positivity = False
        shrinkwrap_minstep = 0
        shrinkwrap_threshold_fraction = 0.1
        enforce_connectivity_steps = [300]
        recenter_steps = list(range(501,8002,500))
    else:
        args.mode = "None"

    if args.filemultiple != None:

        #all the arguments that can be unique for each file
        include_args = [
        'positivity', 'negativity', 'flatten_low_density', 'minimum_density', 'maximum_density', 'ncs',
        'ncs_axis', 'ncs_steps', 'recenter', 'recenter_steps', 'recenter_mode', 'shrinkwrap',
        'shrinkwrap_minstep', 'shrinkwrap_iter', 'enforce_connectivity', 'enforce_connectivity_steps',
        'limit_dmax', 'limit_dmax_steps', 'mode', 'average_weights', 'buffer_scattering_length_densities']
        filemultiple = args.filemultiple
        #initialize all modifiable parameters as lists of length 'k' (number of contrast files)
        for temp_param in include_args:
            #couldn't generalize around the dot operator so I caved and listed them
            if temp_param == "positivity":
                args.positivity = [args.positivity]*len(filemultiple)
            if temp_param == "negativity":
                args.negativity = [args.negativity]*len(filemultiple)
            if temp_param == "flatten_low_density":
                args.flatten_low_density = [args.flatten_low_density]*len(filemultiple)
            if temp_param == "minimum_density":
                args.minimum_density = [args.minimum_density]*len(filemultiple)
            if temp_param == "maximum_density":
                args.maximum_density = [args.maximum_density]*len(filemultiple)
            if temp_param == "ncs":
                args.ncs = [args.ncs]*len(filemultiple)
            if temp_param == "ncs_axis":
                args.ncs_axis = [args.ncs_axis]*len(filemultiple)
            if temp_param == "ncs_steps":
                args.ncs_steps = [args.ncs_steps]*len(filemultiple)
            if temp_param == "recenter":
                args.recenter = [args.recenter]*len(filemultiple)
            if temp_param == "recenter_steps":
                args.recenter_steps = [args.recenter_steps]*len(filemultiple)
            if temp_param == "recenter_mode":
                args.recenter_mode = [args.recenter_mode]*len(filemultiple)
            if temp_param == "shrinkwrap":
                args.shrinkwrap = [args.shrinkwrap]*len(filemultiple)
            if temp_param == "shrinkwrap_minstep":
                args.shrinkwrap_minstep = [args.shrinkwrap_minstep]*len(filemultiple)
            if temp_param == "shrinkwrap_iter":
                args.shrinkwrap_iter = [args.shrinkwrap_iter]*len(filemultiple)
            if temp_param == "enforce_connectivity":
                args.enforce_connectivity = [args.enforce_connectivity]*len(filemultiple)
            if temp_param == "enforce_connectivity_steps":
                args.enforce_connectivity_steps = [args.enforce_connectivity_steps]*len(filemultiple)
            if temp_param == "limit_dmax":
                args.limit_dmax = [args.limit_dmax]*len(filemultiple)
            if temp_param == "limit_dmax_steps":
                args.limit_dmax_steps = [args.limit_dmax_steps]*len(filemultiple)
            if temp_param == "mode":
                args.mode = [args.mode]*len(filemultiple)
            if temp_param == "average_weights" and args.average_weights == None: #hasn't been initialized yet
                args.average_weights = [1.0]*len(filemultiple) #gives equal weights per file by default
            if temp_param == "buffer_scattering_length_densities" and args.buffer_scattering_length_densities == None:
                args.buffer_scattering_length_densities = [0.0]*len(filemultiple)
        if args.parameter_file != None:
            ## There is a parameter file to edit individual parameters
            args.parameter_file = [foo.strip() for foo in args.parameter_file]
            temp_cfile_param_list = []
            for k in range(len(args.parameter_file)):
                cfname, cf_args = param_file_handler(args.parameter_file[k])
                temp_cfile_param_list.append(cf_args)
            for k in range(len(temp_cfile_param_list)): #list of dictionaries with indices corresponding to contrast files
                for para, argval in temp_cfile_param_list[k].items():
                    if para == "positivity":
                        args.positivity[k] = argval
                    if para == "negativity":
                        args.negativity[k] = argval
                    if para == "flatten_low_density":
                        args.flatten_low_density[k] = argval
                    if para == "minimum_density":
                        args.minimum_density[k] = argval
                    if para == "maximum_density":
                        args.maximum_density[k] = argval
                    if para == "ncs":
                        args.ncs[k] = argval
                    if para == "ncs_axis":
                        args.ncs_axis[k] = argval
                    if para == "ncs_steps":
                        args.ncs_steps[k] = argval
                    if para == "recenter":
                        args.recenter[k] = argval
                    if para == "recenter_steps":
                        args.recenter_steps[k] = argval
                    if para == "recenter_mode":
                        args.recenter_mode[k] = argval
                    if para == "shrinkwrap":
                        args.shrinkwrap[k] = argval
                    if para == "shrinkwrap_minstep":
                        args.shrinkwrap_minstep[k] = argval
                    if para == "shrinkwrap_iter":
                        args.shrinkwrap_iter[k] = argval
                    if para == "enforce_connectivity":
                        args.enforce_connectivity[k] = argval
                    if para == "enforce_connectivity_steps":
                        args.enforce_connectivity_steps[k] = argval
                    if para == "limit_dmax":
                        args.limit_dmax[k] = argval
                    if para == "limit_dmax_steps":
                        args.limit_dmax_steps[k] = argval
                    if para == "average_weights":
                        args.average_weights[k] = argval
                    if para == "buffer_scattering_length_densities":
                        args.buffer_scattering_length_densities[k] = argval
                    if para == "mode":
                        if argval[0].upper() == "F":
                            args.mode[k] = "FAST"
                            nsamples = 32
                            args.shrinkwrap_minstep[k] = 1000
                            args.enforce_connectivity_steps[k] = [2000]
                            args.recenter_steps[k] = range(501,2502,500)
                        elif argval[0].upper() == "S":
                            args.mode[k] = "SLOW"
                            nsamples = 64
                            args.shrinkwrap_minstep[k] = 5000
                            args.enforce_connectivity_steps[k] = [6000]
                            args.recenter_steps[k] = range(501,8002,500)
                        elif argval[0].upper() == "M":
                            args.mode[k] = "MEMBRANE"
                            nsamples = 64
                            args.positivity[k] = False
                            args.shrinkwrap_minstep[k] = 0
                            shrinkwrap_threshold_fraction = 0.1
                            args.enforce_connectivity_steps[k] = [300]
                            args.recenter_steps[k] = range(501,8002,500)
                        else:
                            args.mode[k] = "None"

        else: #no parameter file specified - default input
            for k in range(len(filemultiple)):
                if args.mode[k][0].upper() == "F":
                    args.mode[k] = "FAST"
                    nsamples = 32
                    args.shrinkwrap_minstep[k] = 1000
                    args.enforce_connectivity_steps[k] = [2000]
                    args.recenter_steps[k] = range(501,2502,500)
                elif args.mode[k][0].upper() == "S":
                    args.mode[k] = "SLOW"
                    nsamples = 64
                    args.shrinkwrap_minstep[k] = 5000
                    args.enforce_connectivity_steps[k] = [6000]
                    args.recenter_steps[k] = range(501,8002,500)
                elif args.mode[k][0].upper() == "M":
                    args.mode[k] = "MEMBRANE"
                    nsamples = 64
                    args.positivity[k] = False
                    args.shrinkwrap_minstep[k] = 0
                    shrinkwrap_threshold_fraction = 0.1
                    args.enforce_connectivity_steps[k] = [300]
                    args.recenter_steps[k] = range(501,8002,500)
                else:
                    args.mode[k] = "None"

    #allow user to explicitly modify those values by resetting them here to the user defined values
    if args.nsamples is not None:
        nsamples = args.nsamples

    if args.shrinkwrap_minstep is not None:
        shrinkwrap_minstep = args.shrinkwrap_minstep

    if args.shrinkwrap_threshold_fraction is not None:
        shrinkwrap_threshold_fraction = args.shrinkwrap_threshold_fraction

    if args.enforce_connectivity_steps is not None:
        enforce_connectivity_steps = args.enforce_connectivity_steps
    if not isinstance(enforce_connectivity_steps, list):
        enforce_connectivity_steps = [ enforce_connectivity_steps ]

    if args.recenter_steps is not None:
        recenter_steps = args.recenter_steps
    if not isinstance(recenter_steps, list):
        recenter_steps = [ recenter_steps ]

    if args.limit_dmax_steps is not None:
        limit_dmax_steps = args.limit_dmax_steps
    else:
        limit_dmax_steps = [502]
    if not isinstance(limit_dmax_steps, list):
        limit_dmax_steps = [ limit_dmax_steps ]

    if args.steps is not None:
        steps = args.steps

    if args.dmax is not None and args.dmax >= 0:
        dmax = args.dmax
    elif gnomdmax is not None:
        dmax = gnomdmax
    else:
        dmax = 100.0

    if args.voxel is None and nsamples is None:
        voxel = 5.
    elif args.voxel is None and nsamples is not None:
        voxel = dmax * args.oversampling / nsamples
    else:
        voxel = args.voxel

    #now recollect all the edited options back into args
    args.nsamples = nsamples
    args.shrinkwrap_minstep = shrinkwrap_minstep
    args.shrinkwrap_threshold_fraction = shrinkwrap_threshold_fraction
    args.enforce_connectivity_steps = enforce_connectivity_steps
    args.recenter_steps = recenter_steps
    args.limit_dmax_steps = limit_dmax_steps
    args.dmax = dmax
    args.voxel = voxel

    return args
