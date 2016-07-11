import ast
import eccpy.settings as eccpysettings
import eccpy.tools as tools
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import sys

def compare_rawdata(settings_excel_file, sample_names, **kwargs):
    """ Compare raw dose-response curves between selected samples, for selected experiments.

    Processes the datafiles in settings excel file marked as "TRUE" for "run gatherer."
    Collects output from the "run_curvefit" program, but only for the selected samples.
    Recreates the fitted curves from the four-parameter Hill equation, with the previously calculated hill_constants.

    The output of the compare_rawdata is saved in the same folder as the "run_gatherer":
    ORIGINAL_SUBFOLDER_WITH_SETTINGS_EXCEL_FILE/analysed/todays_date/

    Running this script will overwrite any previous files with the same name (i.e., analysed on the same day)

    Parameters
    ----------
    settings_excel_file : String
        Path to the settings file containing the list of datafiles for analysis, and also chosen parameters
    sample_names : list of strings, or list of tuples
        For the output of one figure comparing the raw data of samples:
            sample_names will be a list of sample names
            e.g. sample_names=["control_sample","sample1", "sample2]"
        For the output of multiple figures comparing the raw data of samples:
            sample_names will be a tuple of lists of of sample names
            e.g. sample_names=(["control_sample","sample1"], ["control_sample","sample2]")
            Note: for many output figures, creating a list_output_fig_names is recommended.
        The strings (e.g. "sample1") must resemble the sample names as written in the original data files,
        rather than shortened names (e.g. "s1") according to the list in the settings file.

    Keyword Arguments
    -----------------
    list_output_fig_names : list of strings
        List of output figure names that is used for a multiple analyses, where the sample_names is a list of tuples.
        e.g.
        sample_names=(["control_sample","sample1"], ["control_sample","sample2]", ["sample1","sample2]")
        list_output_fig_names=["control_vs_sample1", "control_vs_sample2", "sample1_vs_sample2"]

    Saved Files and Figures
    -------
    Dose_Response_Curve_Comparison : Scattergram (orig datapoints), line chart (fitted curves)
        Dose-response curves for the selected samples (from sample_names), extracted from the selected experiments
        as labelled "True" in the settings file.
        Each unique sample is given a different colour.
        x-axis : dose units
        y-axis : response units
        scattergram : original datapoints from experiment
        line chart : fitted curve used to calculate EC50, recreated from saved hill_constants

    Note
    -------
    The compare_rawdata function is best used for the comparison of 2-3 samples, with 5-10 experiments.
    It currently can accept 8 different samples.
    Increasing the number of samples or experiments is likely to result in a very cluttered graph.
    """
    print("\nStarting compare_rawdata program")

    # if there is a list of output figure names in the keyword arguments, check that the length matches the sample_names
    if "list_output_fig_names" in kwargs.keys():
        list_output_fig_names = kwargs["list_output_fig_names"]
        # check that the list of sample tuples and list of names has the same length
        if len(list_output_fig_names) != len(sample_names):
            raise IndexError("The length of sample_names does not match the "
                             "list_output_fig_names. Please check your list of samples.")

    # create output folder and output file basename
    outpath, basename = eccpysettings.setup_output_folder(settings_excel_file, "compare_raw")
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    # setup tableau20 colour list
    t20 = tools.setup_t20_colour_list()
    # add black (k) to the front of the list
    t20.insert(0,"0.5")
    # add the relevant paths to the data files to the dataframe for files (dff)
    settings, dff, df_samplenames = eccpysettings.read_settings_file(settings_excel_file)
    # create a list of unique markers for the scattergram
    markerlist = [".",",","o","v","^","<",">","1","2","3","4","8","s","p","*","h","H","+","x","D","d","|","_"]
    # extend the list, in the unlikely case that someone has many replicates
    markers = markerlist + markerlist + markerlist
    # set transparency of datapoints
    alpha = 0.5
    # set default fontsize
    plt.rcParams['font.size'] = 6
    # define xycoordinates for later annotations
    xyc = "axes fraction"
    # extract list of adjusted datasets for analysis
    datasets = ast.literal_eval(settings["datasets"])

    # create boolean
    at_least_one_sample_found_in_selected_datafiles = False

    if not isinstance(sample_names[0], tuple):
        # the list of sample names contains only strings, and therefore only a single raw analysis is performed
        # convert to a list containing only one tuple, with the sample names for comparison
        sample_names = [tuple(sample_names)]

    for s_num, sample_tuple in enumerate(sample_names):
        if True in list(dff.loc[:, "run gatherer"]):
            n_files_to_analyse = dff.loc[dff["run gatherer"] == True].shape[0]
            for d in datasets:
                # change the dataset name (e.g. "_orig" to "") to an empty string if there is only one dataset for analysis
                d_name = "" if len(datasets) == 1 else d
                # close any open plots
                plt.close("all")
                # create new canvas (figure) containing a single plot (ax)
                fig, ax = plt.subplots()
                # create a counter for the number of files
                file_counter = 0
                # iterate through all of the data files labeled for analysis
                for fn in dff.loc[dff["run gatherer"] == True].index:
                    file_counter += 1
                    # print a dot for each file analysed, for each sample name
                    sys.stdout.write(".")
                    sys.stdout.flush()
                    # open output summary file with LD50 values as pandas dataframe
                    ofd_EC50_eval_csv = dff.loc[fn,"ofd_EC50_eval_csv"]
                    if os.path.isfile(ofd_EC50_eval_csv):
                        filename = os.path.split(ofd_EC50_eval_csv)[1]
                        df = pd.read_csv(ofd_EC50_eval_csv)
                        # set the index as the sample_name (long name)
                        df.set_index("sample_name", inplace=True)
                        # redefine to only include data that is labelled as "data_seems_okay"
                        df = df.loc[df["data_seems_okay{}".format(d)] == True]
                        sample_counter = 0
                        for sample_name in sample_tuple:
                            # counter = 0
                            if sample_name in df.index:
                                at_least_one_sample_found_in_selected_datafiles = True
                                # obtain the bool, or series of bools that say if the data is okay
                                data_seems_okay_X = df.loc[sample_name,"data_seems_okay{}".format(d)]
                                # if it's not a series, the sample name was only found once in that experiment
                                if not isinstance(data_seems_okay_X, pd.Series):
                                    # counter += 1
                                    # convert the x_orig data from a stringlist to a numpy array
                                    x = np.array(ast.literal_eval(df.loc[sample_name,"x{}".format(d)]))
                                    # convert the y_orig data from sample_name stringlist to a numpy array
                                    y = np.array(ast.literal_eval(df.loc[sample_name,"y{}".format(d)]))
                                    # plot the datapoints for that set of data
                                    if sample_counter == 0:
                                        # if it's the first datapoint from that file, set a label for the legend
                                        ax.scatter(x, y, color = t20[sample_counter], s=15, alpha=alpha,
                                                   marker=markers[file_counter], label=filename[:16])
                                    else:
                                        # otherwise, do not write another legend label
                                        ax.scatter(x, y, color = t20[sample_counter], s=15, alpha=alpha,
                                                   marker=markers[file_counter], label="_nolabel_")
                                    # retrieve the hill constants for the curve
                                    hill_constants = ast.literal_eval(df.loc[sample_name,"hill_constants{}".format(d)])
                                    # create 500 datapoints on the x-axis to plot the curve
                                    x_fitted_norm = np.linspace(0, 1, 500)
                                    # create the y datapoints using the sigmoid equation
                                    y_fitted_norm = tools.hill_eq(hill_constants, x_fitted_norm)
                                    # denormalise the x datapoints to the original concentrations
                                    x_fitted = tools.denormalise_0_1(x_fitted_norm, x.min(), x.max())
                                    # denormalise the y datapoints to the original concentrations
                                    y_fitted = tools.denormalise_0_1(y_fitted_norm, y.min(), y.max())
                                    # plot the curve of the fitted data, using the same colours as the datapoints
                                    ax.plot(x_fitted, y_fitted, color = t20[sample_counter], alpha=alpha)
                                    # sample_counter += 1

                                # if it is a series, the sample name was found more than once in that experiment
                                elif isinstance(data_seems_okay_X, pd.Series):
                                    # retrieve the list of x datapoints, y datapoints, and hill constants from curve
                                    x_list_replicates = list(df.loc[sample_name,"x{}".format(d)])
                                    y_list_replicates = list(df.loc[sample_name,"y{}".format(d)])
                                    hill_constants_reps = list(df.loc[sample_name,"hill_constants{}".format(d)])
                                    for i in range(len(x_list_replicates)):
                                        # counter += 1
                                        # convert the x, y and hill constants from a stringlists to numpy arrays
                                        x = np.array(ast.literal_eval(x_list_replicates[i]))
                                        y = np.array(ast.literal_eval(y_list_replicates[i]))
                                        hill_constants = np.array(ast.literal_eval(hill_constants_reps[i]))
                                        # plot the datapoints for that set of data
                                        if sample_counter == 0:
                                            # if it's the first datapoint from that file, set a label for the legend
                                            ax.scatter(x, y, color = t20[sample_counter], s=15, alpha=alpha,
                                                       marker=markers[file_counter], label=filename[:8])
                                        else:
                                            # otherwise, do not write another legend label
                                            ax.scatter(x, y, color = t20[sample_counter], s=15, alpha=alpha,
                                                       marker=markers[file_counter], label="_nolabel_")
                                        # create 500 datapoints on the x-axis to plot the curve
                                        x_fitted_norm = np.linspace(0, 1, 500)
                                        # create the y datapoints using the sigmoid equation
                                        y_fitted_norm = tools.hill_eq(hill_constants, x_fitted_norm)
                                        # denormalise the x datapoints to the original concentrations
                                        x_fitted = tools.denormalise_0_1(x_fitted_norm, x.min(), x.max())
                                        # denormalise the y datapoints to the original concentrations
                                        y_fitted = tools.denormalise_0_1(y_fitted_norm, y.min(), y.max())
                                        # plot the curve of the fitted data, using the same colours as the datapoints
                                        ax.plot(x_fitted, y_fitted, color = t20[sample_counter], alpha=alpha)
                                else:
                                    raise TypeError("data_seems_okay_X is neither bool nor series")
                            sample_counter += 1


                if not at_least_one_sample_found_in_selected_datafiles:
                    raise ValueError("No samples found in the selected datasets!\nSamples: {}".format(sample_names))
                xaxis_pos = 0.02
                yaxis_pos = np.linspace(0.95,0.7,8)
                for n, sample_name in enumerate(sample_tuple):
                    ax.annotate(s=sample_name,  xy=(xaxis_pos,yaxis_pos[n]),
                                xycoords=xyc,
                                color = t20[n])
                ymin, ymax = ax.get_ylim()
                ax.set_ylim(ymin,ymax*1.3)
                xmin, xmax = ax.get_xlim()
                # ax.set_xlim(-10,xmax*1.1)
                ax.set_xlim(xmin - xmax * 0.1, xmax * 1.1)
                # ax.set_xlim(-10, 200)
                ax.legend(ncol=2, scatterpoints=1)
                if "list_output_fig_names" in kwargs.keys():
                    # set figure name "fig_name" for saving.
                    fig_name = list_output_fig_names[s_num]
                else:
                    # If a list of tuple names is not given, use the sample_tuple number, "n"
                    fig_name = s_num
                ax.set_title("comparison of raw data for selected samples ({e} experiments),  "
                             "{b} {c}".format(b=d_name,c=os.path.split(settings_excel_file)[1],e=n_files_to_analyse))
                # set xlabel, ylabel
                ax.set_xlabel("{a} ({b})".format(a=settings["x-axis (dose) label"],b=settings["x-axis (dose) units"]))
                ax.set_ylabel(settings["y-axis (response) label"],rotation='vertical')
                # save the figure in png format
                figpath = os.path.join(outpath, "{b}_{n}{d}.png".format(b=basename,n=fig_name,d=d_name))
                fig.savefig(figpath, format = "png", dpi = 150)
                plt.close("all")
    print('\nComparison of raw data is finished.\nOutput files are saved in the following directory:\n{}'.format(outpath))