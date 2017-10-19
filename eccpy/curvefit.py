#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ECCpy is a user-friendly program for the high-throughput calculation of EC50 values.

Copyright (C) 2016  Mark George Teese

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq, leastsq
import os
import csv
import sys
import ast
import eccpy.settings as eccpysettings
import eccpy.judgefit as judgefit
import eccpy.tools as tools

def run_curvefit(settings_excel_file):
    """ Prepares dose-response curves and calculates EC50 values for samples in a single experiment.

    Processes the datafiles in settings excel file marked as "TRUE" for "run curvefit."
    Fits sigmoidal curves to data using the Hill equation.
    Calculates the EC50 value.

    Running this script will overwrite any previous output files with the same names (i.e., analysed from data
    in the same folder)

    Parameters
    ----------
    settings_excel_file : str
        Path of settings file containing the list of datafiles for analysis, and also chosen user parameters.

    Saved Files and Figures
    -----------------------
    fig0_single_sample_png : png (also as pdf)
        Fitted sigmoidal curve for an individual sample. Contains details regarding the
        automatic judement as to whether the "data_needs_checking" or whether the data seems okay. Contains curves
        for original data, and also adjusted datasets (e.g. "fixed upper limit", or in future "outliers removed")
    EC50_analysis_fig_01 : png (also as pdf)
        Scattergram with all data and sigmoidal curves for that experiment.
        Bar-chart with all calculated EC50 values for that experiment, sample letter only on x-axis.
    EC50_analysis_fig_02 : png (also as pdf)
        Bar-chart with all calculated EC50 values for that experiment,
        full name on x-axis
    ofd_EC50_eval_csv : csv
        Output summary csv file with EC50 values.
        Comma separation, English numbering system, with values in quotation marks to increase compatibility.
    ofd_EC50_eval_excel : excel file
        Output summary excel file with EC50 values

    Usage
    -----
    import eccpy
    settings = r"C:\path\to\your\settings\file.xlsx"
    eccpy.run_curvefit(settings)

    Note
    -------
    Automatically runs the "judge_fit" program, which tries to automatically determine if the data is good,
    identify problems, and decide whether the EC50 value is accurate.
        - EC50 values judged to be okay will NOT be labeled as "data_needs_checking", and will be collected in the
          analysis scripts, which compare results between different days or experiments
        - EC50 values labeled as "data_needs_checking" will NOT be collected in the analysis scripts.

    EC50 values are calculated after fitting the data to a four-parameter sigmoid equation (Hill equation).
    For details regarding the EC50 calculation and equation, see Wikipedia (https://en.wikipedia.org/wiki/EC50).
    The EC50 is calculated by finding the dose (x value) at the position halfway between the min and max response.
    This is calculated by root-finding using the brent equation, which is more reliable than the EC50 parameter
    in the Hill equation.
    """

    # add the relevant paths to the data files to the dataframe for files (dff)
    settings, dff, df_samplenames = eccpysettings.read_settings_file(settings_excel_file)
    # create t20 colour list
    t20 = tools.setup_t20_colour_list()

    """
    ANALYSE THE RAW DATA
    """
    # if any of the files are labelled True for "run curvefit"
    if True in list(dff.loc[:, "run curvefit"]):
        print("Starting run_curvefit program for selected samples.\n")
        for fn in dff.loc[dff["run curvefit"] == True].index:
            calc_EC50(fn, dff, settings, t20)
    else:
        print("None of the datafiles are marked TRUE for 'run curvefit'. Suggest checking the excel settings file.")
    print("\nrun_curvefit program is finished.")


def calc_EC50(fn, dff, settings, t20):
    """ Calculates multiple EC50 values from a single input data file.

    Parameters
    ----------
    fn : int
        File number in list of files to analyse.
    dff : pandas DataFrame
        Dataframe for Files. Contains all the paths for input and output files.
        Created from the "files" tab of the settings excel file.
        Contains the "True" / "False" list of input files to analyse in that run.
    settings : pandas Series
        Parameter settings for fitting, and data analysis
    t20 : list
        Tableau20 list of colours

    Saved Files
    -------
    fig0_single_sample_png : png (also as pdf)
        Fitted sigmoidal curve for an individual sample. Contains details regarding the
        automatic judement as to whether the "data_needs_checking" or whether the data seems okay. Contains curves
        for original data, and also adjusted datasets (e.g. "fixed upper limit", or in future "outliers removed")
    EC50_analysis_fig_01 : png (also as pdf)
        Scattergram with all data and sigmoidal curves for that experiment.
        Bar-chart with all calculated EC50 values for that experiment, sample letter only on x-axis.
    EC50_analysis_fig_02 : png (also as pdf)
        Bar-chart with all calculated EC50 values for that experiment,
        full name on x-axis
    ofd_EC50_eval_csv : csv
        Output summary csv file with EC50 values.
        Comma separation, English numbering system, with values in quotation marks to increase compatibility.
    ofd_EC50_eval_excel : excel file
        Output summary excel file with EC50 values
    """
    # define datasets that have been adjusted before attempting fitting ("_orig" is default, "_ful", "fixed upper limit"
    #  is for specific LD50 analyses
    datasets = ast.literal_eval(settings["datasets"])

    # extract the data file path, method (e.g. EC50), expected curveshape (S or Z), etc. from settings file
    data_file = dff.loc[fn, "response data file"]
    print(data_file)
    method = "{ct}{pr}".format(ct=settings["calculation_type"],
                               pr=str(settings["percentage_response"]))
    dose_response_curveshape = settings["dose_response_curveshape"]
    doselabel = settings["x-axis (dose) label"]
    doseunits = settings["x-axis (dose) units"]

    # create new output file directories, if they don't exist already
    dir_columns = ["output_folder", "ofd_csv", "ofd_curves"]
    if settings["save_as_pdf"] in (True, "TRUE"):
        dir_columns.append("ofd_pdfs")
    for column in dir_columns:
        if not os.path.exists(dff.loc[fn, column]):
            os.makedirs(dff.loc[fn, column])

    # examine the input file to confirm integrity, correct datatype, etc
    dff = examine_input_datafile(fn, dff)
    if dff.loc[fn, "resp_datafile_ok"] != True:
        raise tools.DatafileError("The response datafile is not readable, or in the incorrect format."
                            "\nFile affected : \n{}".format(data_file))

    # obtain the details regarding the data format
    eccpy, resp_datafileformat, resp_machinetype, resp_assaytype = dff.loc[fn, "response dataformat"].split("|")

    # for the 96-well samples, obtain path for file with the dose concentrations and sample names
    dose_conc_excelfile = dff.loc[fn, "dose conc file"]
    # replace np.nan with an empty string, if the sample is from the 12-well platereader and no excel file is given
    if isinstance(dose_conc_excelfile, float):
        if np.isnan(dose_conc_excelfile):
            dose_conc_excelfile = ""
    # define path to file with dose concentrations
    dose_conc_excel_path = os.path.join(dff.loc[fn, "input file directory"], dose_conc_excelfile)

    if resp_machinetype == "versamax" and dff.loc[fn, "resp_datafile_ok"] == True:
        # read versamax text file, convert to dataframe
        df_resp_orig, df_resp_all, df_dose_orig = read_versamax_txt_datafile(fn, dff, resp_assaytype)
    else:
        df_dose_orig = "not yet created"
        df_resp_all = "not yet created"
        # do nothing. the response data file is not from VersaMax, is probably from the 12-well platereader
        pass
    data_file_path = dff.loc[fn,"data_file_path"]
    # standardise the dose concentration data so that samples are listed in index(vertical), and doses horizontal
    df_dose_all, df_resp_all = standardise_doseconc_data(fn, dff, df_dose_orig, df_resp_all, data_file_path)

    if resp_machinetype == "versamax":

        # create empty dataframe to hold the distributed datapoints
        df_resp_all = pd.DataFrame()

        if resp_assaytype in ["8dose12sample", "12dose8sample", "24dose4sample"]:
            # double-check that the dose file is really an excel file
            if dose_conc_excel_path[-4:] not in [".xls", "xlsx"]:
                raise ValueError("File with dose concentrations does not end in .xls or .xlsx. ({}) Please check settings file.".format(dose_conc_excel_path))
            # define the relevant column name with the response values
            col_resp = "MeanAbsorb"
        elif "ampconc" in resp_assaytype:
            # DEPRECATED. Use instead ["8dose12sample", "12dose8sample", "24dose4sample"]
            # define the relevant column name with the response values
            col_resp = "MeanOD600"

        # iterate through the rows, dispersing the datapoints throughout the new dataframe df_resp_all
        for row in df_resp_orig.index:
            # extract the sample name (e.g. A01, or AA01)
            Sample_Name = df_resp_orig.loc[row, "Sample"]
            # split the Sample_Name into components, sample letter and dose number (e.g. AA, 01)
            sLet_versamax = Sample_Name[:-2]
            dosenum = Sample_Name[-2:]
            # now distribute each response datapoint in a new dataframe,
            # where the sLet is the index and the dosenumber the columns
            # index_in_df_resp_all = df_resp_orig.loc[row, "index"]
            # column_in_df_resp_all = df_resp_orig.loc[row, "column"]
            response_datapoint = df_resp_orig.loc[row, col_resp]
            # df_resp_all.loc[index_in_df_resp_all,column_in_df_resp_all] = y_value_response
            df_resp_all.loc[sLet_versamax, dosenum] = response_datapoint

        # create new DataFrame for Sample names (dfS)
        dfS = pd.read_excel(dose_conc_excel_path, sheetname="samples", index_col=0)
        # replace original index with the resp_assaytype
        dfS["orig_index"] = dfS.index
        assert resp_assaytype in dfS.columns
        dfS.set_index(resp_assaytype, drop=False, inplace=True)

        # Match "Contains_Data" between the samples tab (dfS), and the XxdoseYsample tab (df_dose_all)
        # First, convert all "true-like" to python bool. (and anything else to False)
        dfS["Contains_Data"] = dfS["Contains_Data"].apply(tools.convert_truelike_to_bool)
        df_dose_all["Contains_Data"] = df_dose_all["Contains_Data"].apply(tools.convert_truelike_to_bool)
        # The lengths won't match, but this doesn't mean much. Find the shortest list.
        len_df_dose_all, len_dfS = df_dose_all.shape[0], dfS.shape[0]
        min_num_datapoints = np.min([len_df_dose_all, len_dfS])
        # extract list of bool from both lists,truncated to have the same length
        df_dose_all_Contains_Data = list(df_dose_all["Contains_Data"])[:min_num_datapoints]
        dfS_Contains_Data = list(dfS["Contains_Data"])[:min_num_datapoints]
        # Third, check if these two lists match
        matching = df_dose_all_Contains_Data == dfS_Contains_Data
        if matching:
            # add the "Contains_Data" col from dose excel file to dataframe containing response data (from txt datafile)
            df_resp_all["Contains_Data"] = df_dose_all_Contains_Data
        else:
            raise tools.DataMismatchError("\n\nThe 'Contains_Data' columns/rows are not matching in {a} & samples tabs."
                                    "\n\nDouble-check excel file.\n\nAffected file :\n{p}\n\n"
                                    "{b} tab:\n{c}\n\nsamples tab:\n{d}\n\n".format(a=resp_assaytype,
                                                                                p=dose_conc_excel_path,
                                                                                b=resp_assaytype,
                                                                                c=df_dose_all_Contains_Data,
                                                                                d=dfS_Contains_Data))
        # create a list of the samples that is the same number of rows as df_dose_all (original dfS index is ignored)
        n_rows_df_dose_all = df_dose_all.shape[0]
        # transfer sample names to dataframes with dose concentrations & response values
        series_of_sample_names_with_standard_AA_AB_index = dfS.loc[dfS.Contains_Data]["samples"]
        df_dose_all['samples'] = series_of_sample_names_with_standard_AA_AB_index
        df_resp_all['samples'] = series_of_sample_names_with_standard_AA_AB_index

    # create a view on the dataframes, so that it only shows the microplate data manually marked as "Contains_Data" = True
    dfdose = df_dose_all[df_dose_all.Contains_Data == True].copy()
    dfresp = df_resp_all[df_resp_all.Contains_Data == True].copy()
    dict_dfe = {}

    # determine the longest sample name
    sample_name_len_max = df_dose_all.samples.str.len().max()

    #set the fontsize for the figure
    fig_fontsize = 6
    #set the default font for the figures
    plt.rc('font', family='sans-serif')
    plt.rc('font', serif='Helvetica Neue')
    plt.rc('text', usetex='false')
    plt.rcParams.update({'font.size': fig_fontsize})

    # iterate through all of the samples marked for analysis within the excel file with dose concentrations
    for sNum, sLet in enumerate(dfdose.index):
        # create a new dataframe for the evaluation of the EC50 calculations
        # Create a DataFrame for EC50 data (dfe)
        # There will be a new dfe for each sample number. Each dfe will be added to a dictionary, which is then
        # converted to df_eval
        dfe = pd.DataFrame()
        dfe.loc["sLet", sLet] = sLet
        dfe.loc["sNum", sLet] = sNum
        # obtain the name for that sample
        sample_name = str(dfdose.loc[sLet, "samples"])
        dfe.loc["sample_name", sLet] = sample_name

        # set up the path for the image files to be saved in
        fig0_single_sample_png = os.path.join(dff.loc[fn,"ofd_curves"], "%s " % sLet + sample_name) + ".png"
        fig0_single_sample_pdf  = os.path.join(dff.loc[fn,"ofd_pdfs"], "%s " % sLet + sample_name) + ".pdf"

        #reindex so that only rows that contain data in both dataframes (i.e. dose and response data) are kept for analysis
        #take the index of both dataframes after NaN is removed.
        index_dfs = dfresp.loc[sLet,:].dropna().index
        #Find common elements using the set intersection function.
        cols_with_data_in_both_x_and_y = index_dfs.intersection(dfdose.loc[sLet,:].dropna().index)
        # reindex to drop the columns with text or boolean values
        x_orig = dfdose.loc[sLet,:].reindex(index = cols_with_data_in_both_x_and_y)
        # drop the two text columns, Contains_Data and samples (sample names)
        x_orig.drop(["Contains_Data", "samples"], inplace=True)
        # convert dtype to float
        x_orig = x_orig.astype(float)
        # add the original x values to the output dataframe
        dfe.loc["x",sLet] = list(x_orig)
        # select the original y (response) values
        y_orig = dfresp.loc[sLet,:].reindex(index = cols_with_data_in_both_x_and_y)
        # drop the two text columns, Contains_Data and samples (sample names)
        y_orig.drop(["Contains_Data", "samples"], inplace=True)
        # convert dtype to float
        y_orig = y_orig.astype(float)

        # add to output dataframe
        dfe.loc["y",sLet] = list(y_orig)

        for d in datasets:
            # currently add the ynorm_orig and xnorm_orig to all datasets. Other datasets might require changes later.
            # normalise the x datapoints between 0 and 1 to improve curve fitting
            dfe.loc["xnorm{}".format(d),sLet], dfe.loc["xmin{}".format(d),sLet], dfe.loc["xmax{}".format(d),sLet] = tools.normalise_0_1(x_orig)
            # normalise the y datapoints between 0 and 1 to improve curve fitting
            dfe.loc["ynorm{}".format(d),sLet], dfe.loc["ymin{}".format(d),sLet], dfe.loc["ymax{}".format(d),sLet] = tools.normalise_0_1(y_orig)

        #make an array of >250 datapoints representing the x-axis of the curve
        min = 0
        max = settings["fitted_curve_xaxis_max"]
        n_datapoints = settings["fitted_curve_n_datapoints"]
        dfe.loc["x_fitted_norm", sLet] = np.linspace(min, max, n_datapoints)
        dfe.loc["n_doseconc_tested", sLet] = len(x_orig)

        #######################################################################################################
        #                                                                                                     #
        #                      Create fixed upper limit (ful) dataset for LD50 data                           #
        #          (lowdose horizontal datapoints assumed to have high variation, normalised between          #
        #               two smaller yvalues so data can be fitted to a Z-shaped sigmoidal curve)              #
        #                                                                                                     #
        #######################################################################################################

        if "_ful" in datasets:
            # copy x values to the ful dataset, which utilise the same x-datapoints
            dfe.loc["x_ful", sLet] = dfe.loc["x", sLet]
            dfe.loc["x_fitted_norm_ful", sLet] = dfe.loc["x_fitted_norm", sLet]

            ful_max = settings["ful.yaxis_fixed_upper_limit_max"]
            ful_min = settings["ful.yaxis_fixed_upper_limit_min"]

            # find where the first datapoint drops below the settings["ful.yaxis_fixed_upper_limit_min"]
            if y_orig.min() < ful_min:
                index_y_ful = np.min(np.where(y_orig < ful_min))
            else:
                # there are insufficient lowresponse datapoints(insuff_lowresp_dp). E.g. cells are overgrown.
                # All datapoints are above ful. np.where returns a tuple. Replace index with len(array).
                index_y_ful = len(y_orig)
            # select values above index_y_ful, which are above the fixed upper limit
            y_orig_data_above_ful = y_orig[:index_y_ful]
            y_orig_data_below_ful = y_orig[index_y_ful:]

            if len(y_orig_data_above_ful) > 1:
                # normalise these datapoints between 0 and 1
                y_orig_data_above_ful_norm_0_1 = tools.normalise_0_1(y_orig_data_above_ful)[0]
                # calculate the width for the normalisation of the adjusted datapoints
                normalisation_width = ful_max - ful_min
                # convert the data which is normalised from 0-1, to normalised between the fixed_upper_limit_min and max
                y_orig_data_above_ful_norm = y_orig_data_above_ful_norm_0_1 * normalisation_width + ful_min
                y_ful = np.append(y_orig_data_above_ful_norm, y_orig_data_below_ful)
                # add to output dataframe
                dfe.loc["y_ful",sLet] = list(y_ful)
            else:
                # if there is only one datapoint above the limit, ful adjustment is useless, y_ful = y_orig
                dfe.loc["y_ful",sLet] = dfe.loc["y",sLet]

        # normalise the y (response) data between 0 and 1 for the fitting algorithm
        for d in datasets:
            ynorm, ymin, ymax = tools.normalise_0_1(np.array(dfe.loc["y{}".format(d), sLet]))
            dfe.loc["ynorm{}".format(d), sLet] = ynorm
            dfe.loc["ymin{}".format(d), sLet] = ymin
            dfe.loc["ymax{}".format(d), sLet] = ymax

        #######################################################################################################
        #                                                                                                     #
        #              Setup figure showing the dose-response curves for a single sample                      #
        #                                                                                                     #
        #######################################################################################################

        # close any open figures
        plt.close("all")
        # Create a new figure (i.e. a new 2x2 canvas). If EC50_calculable is false, some errors will be printed.
        fig, axarr = plt.subplots(nrows=2, ncols=2, dpi=300)
        # set an annotation fontsize
        af = 10
        # create a colour list for the various datasets, with selected colours from the tableau20
        co = [t20[0], t20[2], t20[3], t20[4], t20[5], t20[6]]
        # create dictionary to hold the size of the plotted datapoints for each dataset (gradually smaller)
        sd = {}
        # create dictionary to hold the linestyles for the EC50 locators
        ls = {}
        linestylelist = ["dashed", "dotted", "-.", "dashed", "dotted", "-.","dashed", "dotted", "-.",]
        for n,d in enumerate(datasets):
            size_start = 15
            size_increment = -7
            sd[d] = size_start + size_increment * n
            ls[d] = linestylelist[n]
        # set alpha (transparency) level for plotted lines and data in general
        al = 0.8
        # define xycoordinates for later annotations
        xyc = "axes fraction"

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        #                                                                                                     #
        #                        Dose-Respose Curve Fig01: Raw data                                           #
        #                    (Datapoints only - fitting is plotted later)                                     #
        #                                                                                                     #
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        #            _________
        #           |XXXX|    |      Subplot 0, axarr[0, 0]
        #           |XXXX|____|
        #           |    |    |
        #           |____|____|

        # Determine which datasets are going to be plotted
        datasets = ast.literal_eval(settings["datasets"])

        # set the subplot number on the canvas
        Plot_Nr = 1
        for n, d in enumerate(datasets):
            # change the dataset name (e.g. "_orig" to "") to an empty string if there is only one dataset for analysis
            d_name = "" if len(datasets) == 1 else d
            # add original datapoints as a scattergram
            axarr[0,0].scatter(x_orig, dfe.loc["y{}".format(d), sLet], color=co[n], s=sd[d], label=d_name[1:], alpha=al)
        # set xlabel, ylabel, title, etc
        axarr[0,0].set_xlabel("{a} ({b})".format(a=doselabel, b=doseunits), fontsize = fig_fontsize)
        axarr[0,0].set_ylabel(settings["y-axis (response) label"],rotation='vertical', fontsize = fig_fontsize)
        axarr[0,0].set_title("%s   %s" %(sLet, sample_name), fontsize = fig_fontsize)
        axarr[0,0].grid(True, color = '0.75')
        # set the limit of the y-axis to 1 + 0.1, if all datapoints are very low
        ymin, ymax = axarr[0,0].get_ylim()
        # xmin, xmax = axarr[0,0].get_xlim()
        #set y-axis intercept
        # ylim_min_raw = -0.05
        # define the y-limit for the non-normalised (raw) data as the minimum value minus a percentage of max value
        ylim_min_raw = ymin - ymax * 0.05
        # set the x-axis limits
        # xlim_min_raw = -10
        # xlim_min_raw = xmin# - xmax * 0.01
        # define the x-limit for the non-normalised (raw) data as the minimum value minus a percentage of max value
        xlim_min_raw = x_orig.min() - x_orig.max() * 0.1
        ylim_max_raw = y_orig.max() + 0.1 if y_orig.max() > 1.0 else 1.0
        axarr[0,0].set_ylim(ylim_min_raw, ylim_max_raw)
        # set the x-axis limit so that the legend does not hide too many data points
        # find the maximum dose conc. in the whole experiment for that day
        # maxAC = x_orig.max()
        # obtain the variable altering the extension of the x-axis
        # x_axis_extension_after_ACmax_in_plot1 = dff.loc[fn, "x-axis extension in summary fig_0"]
        # define 110% of the limit of the x-axis as the maximum dose conc.
        # xlim_max_plot1 = maxAC + x_axis_extension_after_ACmax_in_plot1
        # define 110% of the limit of the x-axis as the maximum dose conc.
        xlim_max_plot1 = x_orig.max() * 1.1
        axarr[0,0].set_xlim(xlim_min_raw,xlim_max_plot1)
        # axarr[0, 0].set_xlim(0, xlim_max_plot1)
        if "_ful" in datasets:
            axarr[0,0].annotate(s="original data", xy=(0.71,0.9), fontsize=af, xycoords=xyc, color = t20[0])
            axarr[0,0].annotate(s="fixed upper limit (ful) data", xy=(0.43,0.8), fontsize=af, xycoords=xyc, color = t20[2])

        #######################################################################################################
        #                                                                                                     #
        #                     Preliminary judgements whether EC50 is calculable                               #
        #                                                                                                     #
        #######################################################################################################

        # count the number of orig response datapoints above and below the "yaxis upper-lower cutoff" (yulc) value
        dfe.loc["n_resp_dp_below_yulc", sLet] = len(np.where(y_orig < settings["yaxis upper-lower cutoff"])[0])
        dfe.loc["n_resp_dp_above_yulc", sLet] = len(np.where(y_orig > settings["yaxis upper-lower cutoff"])[0])

        if dfe.loc["n_resp_dp_below_yulc", sLet] < settings["min_num_dp_above&below_yulc"]:
            """The cells are "insuff_lowresp_dp" if any of the following are true
                - there are less than two datapoints above the live-dead cutoff value
            """
            for d in datasets:
                dfe.loc["EC50{}".format(d), sLet], dfe.loc["rsquared{}".format(d), sLet] = "insuff_lowresp_dp", 0
                dfe.loc["EC50_hill_eq{}".format(d), sLet], dfe.loc["n_highdose_datapoints{}".format(d), sLet] = 0, 0
                dfe.loc["n_lowdose_datapoints{}".format(d), sLet], dfe.loc["EC50_calculable{}".format(d), sLet] = 0, False
                dfe.loc["residuals_mean{}".format(d),sLet], dfe.loc["y_fitted{}".format(d),sLet] = np.nan, np.nan
                dfe.loc["y_fitted_norm{}".format(d),sLet], dfe.loc["indices_lowdose_datapoints{}".format(d),sLet] = np.nan, np.nan
                dfe.loc["indices_lowdose_datapoints_excl_nearest_EC50{}".format(d),sLet] = np.nan
                dfe.loc["x_fitted{}".format(d), sLet], dfe.loc["response_lowdose_datapoints{}".format(d), sLet] = np.nan, np.nan
                # label the sample as not okay
                dfe.loc["ymin{}".format(d), "%s_okay" % sLet] = False
                dfe.loc["data_seems_okay{}".format(d),sLet] = False

        elif any([y_orig[1] < settings["min_resp_at_2nd_doseconc"],
                 dfe.loc["n_resp_dp_above_yulc", sLet] < settings["min_num_dp_above&below_yulc"],
                 np.array(y_orig)[-2] < settings["min_resp_at_2ndlast_doseconc"]]):
            """For high-throughput LD50 calculations, the cells have "insuff_highresp_dp" if any of the following are true:
                - the y-value of the second datapoint (second dose) is smaller than a fixed minimum value (min_resp_at_2nd_doseconc)
                - there are less than two datapoints above a fixed value (yaxis upper-lower cutoff)
            """
            for d in datasets:
                dfe.loc["EC50{}".format(d), sLet], dfe.loc["rsquared{}".format(d), sLet] = "insuff_highresp_dp", 0
                dfe.loc["EC50_hill_eq{}".format(d), sLet], dfe.loc["n_highdose_datapoints{}".format(d), sLet] = 0, 0
                dfe.loc["n_lowdose_datapoints{}".format(d), sLet], dfe.loc["EC50_calculable{}".format(d), sLet] = 0, False
                dfe.loc["residuals_mean{}".format(d),sLet], dfe.loc["y_fitted{}".format(d),sLet] = np.nan, np.nan
                dfe.loc["y_fitted_norm{}".format(d),sLet], dfe.loc["indices_lowdose_datapoints{}".format(d),sLet] = np.nan, np.nan
                dfe.loc["indices_lowdose_datapoints_excl_nearest_EC50{}".format(d),sLet] = np.nan
                dfe.loc["x_fitted{}".format(d), sLet], dfe.loc["response_lowdose_datapoints{}".format(d), sLet] = np.nan, np.nan
                # label the sample as not okay
                dfe.loc["ymax{}".format(d), "%s_okay" % sLet] = False
                dfe.loc["data_seems_okay{}".format(d),sLet] = False
        else:
            #######################################################################################################
            #                                                                                                     #
            #                               Fit sigmoidal curve to the data                                       #
            #                                                                                                     #
            #######################################################################################################

            #as a starting point, guess the sigmoidal constants
            if dose_response_curveshape == "S":
                hill_constants_guess = (0.0,1.0,0.5,10.0)
            elif dose_response_curveshape == "Z":
                hill_constants_guess = (1.0,0.0,0.5,10.0)

            for d in datasets:
                # change the dataset name (e.g. "_orig" to "") to an empty string if there is only one dataset for analysis
                d_name = "" if len(datasets) == 1 else d
                ynorm = np.array(dfe.loc["ynorm{}".format(d),sLet])
                xnorm = np.array(dfe.loc["xnorm{}".format(d),sLet])
                #use the scipy optimise function to fit a curve to the data points
                hill_constants, cov, infodict, mesg, ier = leastsq(tools.residuals, hill_constants_guess,
                                                                   args=(tools.hill_eq, xnorm, ynorm),
                                                                   full_output=1)

                # save the hill constants for later use
                dfe.loc["hill_constants{}".format(d),sLet] = list(hill_constants)

                # obtain the rsquared value for the fit of the curve to the data
                # code is from http://stackoverflow.com/questions/7588371/scipy-leastsq-goodness-of-fit-estimator
                ss_err = np.sum(np.array(infodict['fvec'])**2)
                ss_tot = np.sum((ynorm-ynorm.mean())**2)
                rsquared = 1 - (ss_err/ss_tot)
                # add the rsquared value to the output dataframe
                dfe.loc["rsquared{}".format(d),sLet] = rsquared

                # also calculate the average residual as a rough "goodness of fit"
                # first apply the optimised function to the original datapoints
                y_fitted_xnorm = tools.hill_eq(hill_constants, xnorm)
                # calculate the residuals, as the (observed y-values)-(fitted y-values). Convert to positive floats.
                residuals_norm = abs(ynorm - y_fitted_xnorm)
                # calculate mean residual
                residuals_norm_mean = residuals_norm.mean()
                # denormalise to the original y-value scale
                residuals_mean = tools.denormalise_0_1(residuals_norm_mean, dfe.loc["ymin{}".format(d), sLet],
                                                 dfe.loc["ymax{}".format(d), sLet])
                # add to output dataframe
                dfe.loc["residuals_mean{}".format(d),sLet] = residuals_mean
                # use mean residual as rough yerr when plotting the EC50 in a barchart. Multiply by 100 to make visible.
                dfe.loc["yerr{}".format(d),sLet] = residuals_mean * 100

                # obtain the constants of the optimised sigmoidal Hill function
                upper, lower, EC50_hill_eq_norm, hillslope = hill_constants

                # denormalise
                xmin = dfe.loc["xmin{}".format(d),sLet]
                xmax = dfe.loc["xmax{}".format(d),sLet]
                dfe.loc["EC50_hill_eq{}".format(d),sLet] = tools.denormalise_0_1(hill_constants[2], xmin, xmax)

                # add to output dataframe
                dfe.loc["EC50_hill_eq_norm{}".format(d), sLet] = EC50_hill_eq_norm
                dfe.loc["hillslope{}".format(d), sLet] = hillslope

                #calculate the value for y for the 1500 points
                x_fitted_norm = np.array(dfe.loc["x_fitted_norm{}".format(d), sLet])
                dfe.loc["y_fitted_norm{}".format(d), sLet] = tools.hill_eq(hill_constants, x_fitted_norm)

                #######################################################################################################
                #                                                                                                     #
                #                          Calculate the EC50 using the fitted curve                                  #
                #                                                                                                     #
                #######################################################################################################

                # obtain the calculation method from the settings file
                method_calc_y50 = settings["method_calc_y50"]
                # obtain the percentage response used for calculation (e.g. 50 ,for EC50)
                percentage_response = settings["percentage_response"]
                # calculate fraction response (i.e. convert EC50 to 0.5)
                fract_response = percentage_response / 100

                if method_calc_y50 == "y50 = (curve_max - curve_min)*0.5 + curve_min":
                    if dose_response_curveshape == "S":
                        # define the x-value for curvemin as the first x-value
                        curvemin = dfe.loc["xnorm{}".format(d),sLet][0]
                        # define the x-value for curvemax as the last x-value
                        curvemax = dfe.loc["xnorm{}".format(d),sLet][-1]
                    if dose_response_curveshape == "Z":
                        # define the x-value for curvemin as the last x-value
                        curvemin = dfe.loc["xnorm{}".format(d),sLet][-1]
                        # define the x-value for curvemax as the first x-value
                        curvemax = dfe.loc["xnorm{}".format(d),sLet][0]
                    # use the hill equation to find the y-value of the curve at these positions
                    y50_curvemin = tools.hill_eq(dfe.loc["hill_constants{}".format(d),sLet], curvemin)
                    y50_curvemax = tools.hill_eq(dfe.loc["hill_constants{}".format(d),sLet], curvemax)
                    # define y50 (yvalue at top of curve - yvalue at bottom of the curve) * 0.5 [or 0.9 for EC90, etc.]
                    y50_norm = (y50_curvemax - y50_curvemin) * fract_response + y50_curvemin

                # detect extended curve formula (e.g. "y50 = (extendedcurve|0.2|_max - extendedcurve|0.2|_min)*0.5")
                elif "extendedcurve" in method_calc_y50 and "|" in method_calc_y50:
                    if dose_response_curveshape == "S":
                        # define the x-value for curvemin as the first x-value
                        curvemin = dfe.loc["xnorm{}".format(d),sLet][0]
                        # define the x-value for curvemax as the last x-value
                        curvemax = dfe.loc["xnorm{}".format(d),sLet][-1]
                    if dose_response_curveshape == "Z":
                        # define the x-value for curvemin as the last x-value
                        curvemin = dfe.loc["xnorm{}".format(d),sLet][-1]
                        # define the x-value for curvemax as the first x-value
                        curvemax = dfe.loc["xnorm{}".format(d),sLet][0]
                    x_range = curvemax - curvemin
                    # extract the extension (e.g. 0.2, 20% from the text string in the settings file)
                    extension_curvemax = float(method_calc_y50.split("|")[1])
                    extension_curvemin = float(method_calc_y50.split("|")[3])
                    if dose_response_curveshape == "S":
                        # define the x-value for curvemin as the first x-value, minus the xrange * extension
                        curvemin = dfe.loc["xnorm{}".format(d),sLet][0] - x_range * extension_curvemin
                        # define the x-value for curvemax as the last x-value
                        curvemax = dfe.loc["xnorm{}".format(d),sLet][-1] + x_range * extension_curvemax
                    if dose_response_curveshape == "Z":
                        # define the x-value for curvemin as the last x-value
                        curvemin = dfe.loc["xnorm{}".format(d),sLet][-1] + x_range * extension_curvemax
                        # define the x-value for curvemax as the first x-value
                        curvemax = dfe.loc["xnorm{}".format(d),sLet][0] - x_range * extension_curvemin
                    # use the hill equation to find the y-value of the curve at these positions
                    y50_curvemin = tools.hill_eq(dfe.loc["hill_constants{}".format(d),sLet], curvemin)
                    y50_curvemax = tools.hill_eq(dfe.loc["hill_constants{}".format(d),sLet], curvemax)
                    # define y50 (yvalue at top of curve - yvalue at bottom of the curve) * 0.5 [or 0.9 for EC90, etc.]
                    y50_norm = (y50_curvemax - y50_curvemin) * fract_response

                elif method_calc_y50 == "y50 = (resp_max - resp_min)*0.5":
                    # define y50 as the centre between min and max datapoints
                    # currently the data is normalised between zero and one, so actually y50 = fract_response
                    y50_norm = (ynorm.max() - ynorm.min()) *  fract_response

                elif method_calc_y50 == "y50 = (resp_end - resp_start)*0.5":
                    y50_norm = (ynorm[-1] - ynorm[0]) * fract_response

                elif method_calc_y50 == "y50 = (resp_start - resp_end)*0.5":
                    y50_norm = (ynorm[0] - ynorm[-1]) * fract_response

                else:
                    raise ValueError("method_calc_y50 ({}) is not recognised. "
                                     "Please check the excel settings file.".format(method_calc_y50))

                # add value to output series
                dfe.loc["y50_norm{}".format(d),sLet] = y50_norm

                #the y-value of 50% cell density is calculated as the middle position in the curve
                #if the curve is perfectly symmetrical, the EC50 should equal the constant 'k' from the hill_constants
                dfe.loc["curve_max_norm{}".format(d),sLet] = dfe.loc["y_fitted_norm{}".format(d),sLet].max()
                dfe.loc["curve_min_norm{}".format(d),sLet] = dfe.loc["y_fitted_norm{}".format(d),sLet].min()

                #dfe.loc["EC50_norm_bq{}".format(d),"%s_okay" % sLet]
                brentq_out_tuple = calc_EC50_brent_eq(sLet, sample_name, dfe.loc["hill_constants{}".format(d), sLet],
                                                      dfe.loc["y50_norm{}".format(d), sLet])

                dfe.loc["EC50_norm_bq{}".format(d), sLet], dfe.loc["EC50_calculable{}".format(d), sLet] = brentq_out_tuple

                # add if the EC50 was calculable to the summary dataframe "okay" column
                if dfe.loc["EC50_calculable{}".format(d), sLet] == True:
                    dfe.loc["EC50_norm_bq{}".format(d),"%s_okay" % sLet] = True
                else:
                    dfe.loc["EC50_norm_bq{}".format(d),"%s_okay" % sLet] = False

                # denormalise the EC50
                dfe.loc["EC50{}".format(d),sLet] = float(tools.denormalise_0_1(brentq_out_tuple[0], xmin, xmax))

                # denormalise the fitted y-values for the curve
                dfe.loc["y_fitted{}".format(d),sLet] = tools.denormalise_0_1(dfe.loc["y_fitted_norm{}".format(d),sLet],
                                                                       dfe.loc["ymin{}".format(d),sLet],
                                                                       dfe.loc["ymax{}".format(d),sLet])

                # denormalise the y50, the y-value used to calculated the EC50
                dfe.loc["y50{}".format(d),sLet] = tools.denormalise_0_1(dfe.loc["y50_norm{}".format(d),sLet],
                                                                  dfe.loc["ymin{}".format(d),sLet],
                                                                  dfe.loc["ymax{}".format(d),sLet])

                # denormalise the fitted data, back to original dose and response concentrations
                dfe.loc["x_fitted{}".format(d),sLet] = tools.denormalise_0_1(x_fitted_norm, xmin, xmax)

            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            #                                                                                                     #
            #                         Dose-Respose Curve Fig01: Fitted curve                                      #
            #                             (raw datapoints plotted earlier)                                        #
            #                                                                                                     #
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            #            _________
            #           |XXXX|    |      Subplot 0, axarr[0, 0]
            #           |XXXX|____|
            #           |    |    |
            #           |____|____|

            # Determine which datasets are going to be plotted
            datasets = ast.literal_eval(settings["datasets"])

            for n, d in enumerate(datasets):
                # change the dataset name (e.g. "_orig" to "") to an empty string if there is only one dataset for analysis
                d_name = "" if len(datasets) == 1 else d
                # add original datapoints as a scattergram
                axarr[0,0].scatter(dfe.loc["x{}".format(d), sLet], dfe.loc["y{}".format(d), sLet], color=co[n], s=sd[d],
                                   label=d_name[1:], alpha=al)
                # plot fitted curve on the subplot with the original y values
                axarr[0, 0].plot(dfe.loc["x_fitted{}".format(d), sLet], dfe.loc["y_fitted{}".format(d), sLet],
                                 '-', color = co[n], alpha=al)
                # extract y50 and EC50 from dataframe
                y50 = dfe.loc["y50{}".format(d), sLet]
                EC50 = dfe.loc["EC50{}".format(d),sLet]
                # draw horizontal line from y50 to EC50
                axarr[0, 0].hlines(y=y50, xmin=xlim_min_raw, xmax=EC50, colors = co[n], linestyles=ls[d], label='', alpha=al)
                # draw vertical line at EC50 from y50
                axarr[0, 0].vlines(x=EC50, ymin=ylim_min_raw, ymax=y50, colors=co[n], linestyles=ls[d])
                """Plot vertical and horizontal lines showing location of EC50 according to the hill equation.
                Not recommended. Hill Eq constant is less reliable than the brent root-finding method.
                EC50_hill = dfe.loc["EC50_hill_eq{}".format(d),sLet]
                axarr[0, 0].hlines(y=y50, xmin=xlim_min_raw, xmax=EC50_hill, colors = co[n], linestyles=ls[d])
                axarr[0, 0].vlines(x=EC50_hill, ymin=ylim_min_raw, ymax=y50, colors=co[n], linestyles=ls[d])"""

            if len(datasets) > 1:
                lg = axarr[0,0].legend(loc=7, ncol=1, scatterpoints=1)
                lg.draw_frame(False)

            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            #                                                                                                     #
            #             Dose-Respose Curve Fig02: Normalised data, with fitted curve.                           #
            #                                                                                                     #
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            #            _________
            #           |    |    |          Subplot 1, axarr[1, 0]
            #           |____|____|
            #           |XXXX|    |
            #           |XXXX|____|

            #set y-axis intercept
            ymin_norm = -0.2
            xmin_norm = -0.05
            for n, d in enumerate(datasets):
                # add normalised datapoints as a scattergram
                axarr[1,0].scatter(dfe.loc["xnorm{}".format(d), sLet], dfe.loc["ynorm{}".format(d), sLet],
                                   color=co[n], s=sd[d], label=d_name[1:])
                # add the fitted curve as a line plot
                axarr[1,0].plot(dfe.loc["x_fitted_norm{}".format(d),sLet], dfe.loc["y_fitted_norm{}".format(d), sLet],
                                '-', color=co[n], alpha=0.8)
                # add horizontal line at y50
                axarr[1,0].hlines(y=dfe.loc["y50_norm{}".format(d), sLet], xmin=xmin_norm, colors = co[n],
                                             xmax=dfe.loc["EC50_norm_bq{}".format(d), sLet], linestyles=ls[d])


                # add vertical line at EC50
                axarr[1,0].vlines(x=dfe.loc["EC50_norm_bq{}".format(d), sLet], ymin=ymin_norm, colors = co[n],
                                             ymax=dfe.loc["y50_norm{}".format(d), sLet], linestyles=ls[d])

            # set xlabel, ylabel, title, grid, etc
            axarr[1,0].set_xlabel("dose concentration (normalised)", fontsize = fig_fontsize)
            axarr[1,0].set_ylabel("response concentration (normalised)",rotation='vertical', fontsize = fig_fontsize)
            axarr[1, 0].text(0.6, 1.1, "normalised data", horizontalalignment='center', fontsize=fig_fontsize)
            axarr[1,0].grid(True, color = '0.75')
            axarr[1,0].set_ylim(ymin_norm, 1.2)
            axarr[1,0].set_xlim(xmin_norm, 1.2)

            if len(datasets) > 1:
                lg = axarr[1,0].legend(loc=7, scatterpoints=1)
                lg.draw_frame(False)
            if "_ful" in datasets:
                # set annotation in top right
                axarr[1,0].annotate(s="normalised (ful) data", xy=(0.53,0.9), fontsize=af, xycoords=xyc, color=t20[2])
                axarr[1,0].annotate(s="normalised data", xy=(0.63,0.8), fontsize=af, xycoords=xyc, color=t20[0])

            # analyse the curve fit and data to judge whether the EC50 value is accurate
            dfe = judgefit.judge_fit(dfe, sLet, settings)
            # dfe_index = pd.Series(dfe.index)
            # dfe_index_ful = dfe_index[dfe_index.apply(lambda x : "_ful" in x)]
            # dfe_index_orig = dfe_index[dfe_index.apply(lambda x : x[-5:] == "_orig")]
            # dfe_ful = dfe.loc[dfe_index_ful]
            # dfe_orig = dfe.loc[dfe_index_orig]

            dict_data_okay = {}
            for d in datasets:
                # create empty list to hold all row indices in dfe that relate to that dataset (based on suffix)
                list_rows_d = []
                for row in dfe.index:
                    # if the suffix in the column matches the dataset (e.g. x_orig[-5:] == "_orig"]
                    if row[-5:] == d:
                        # add the row index label to the list
                        list_rows_d.append(row)
                # reindex the dataframe to contain only the judgement data relevant to that dataset
                dfe_d = dfe.reindex(index=list_rows_d)
                if False in list(dfe_d.loc[:,"{}_okay".format(sLet)]):
                    data_seems_okay = False
                else:
                    # none of the judge_fit judgements suggest the data is not good, label as okay
                    data_seems_okay = True
                # add to a dictionary for quick access
                dict_data_okay[d] = data_seems_okay
                # add the final judgement to the list, for use in the dataframe later
                dfe.loc["data_seems_okay{}".format(d),sLet] = data_seems_okay

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        #                                                                                                     #
        #          Dose-Respose Curve Fig03: Notes on automatic judgement of data quality                     #
        #                                                                                                     #
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        #            _________
        #           |    |XXXX|          Subplot 3, axarr[0, 1]
        #           |____|XXXX|
        #           |    |    |
        #           |____|____|

        Plot_Nr = 3
        # write text on the empty figure
        yaxis_pos = np.linspace(0.9,0.1,10)
        xaxis_left = 0.05
        # set the xaxis position of the annotation for the ful samples
        xful = 0.6
        # set the xaxis position of the annotation for the orig samples
        xori = 0.8
        # # data_evaluation = "data seems good" if data_seems_okay_ful else "ful data needs checking"
        # title_colour = "k" if data_seems_okay_ful else "r"
        title_summ = "Sample %s" % (sLet)
        axarr[0,1].set_title(title_summ, fontsize = af, color = "k", alpha=0.75)

        if resp_machinetype == "versamax" and resp_assaytype in ["8_ampconc", "12_ampconc", "24_ampconc"]:
            if "Uniprot#1" in dfS.columns:
                # combine columns to create the N-terminal vector name
                Nvect = "".join(dfS.loc[sLet,"Uniprot#1":"Mutant#1"])
                # combine columns to create the C-terminal vector name
                Cvect = "".join(dfS.loc[sLet,"Uniprot#2":"Mutant#2"])
                # extract assay version (1.1, 1.2, etc)
                Version = dfS.loc[sLet,"notes"].strip('_')
                axarr[0,1].annotate(s="Nvect: %s" % Nvect, xy=(xaxis_left,yaxis_pos[0]), fontsize=af, xycoords=xyc, alpha=0.75)
                axarr[0,1].annotate(s="Cvect: %s" % Cvect, xy=(xaxis_left,yaxis_pos[1]), fontsize=af, xycoords=xyc, alpha=0.75)
                # axarr[0,1].annotate(s=data_evaluation, xy=(0.6,yaxis_pos[0]), fontsize=af, xycoords=xyc, alpha=0.75)
                axarr[0,1].annotate(s=Version, ha='right', xy=(0.98,yaxis_pos[0]), fontsize=af, xycoords=xyc, alpha=0.75)

        # create dictionary to hold x-axis locations of annotations, for each dataset tested
        xd = {}
        for n,d in enumerate(datasets):
            # start x position
            start = 0.85
            # distance between x positions (will need to make smaller if datasets tested > 2)
            dist = 0.2
            # calculate position for that annotation
            xpos = start - dist * n
            # add to dictionary
            xd[d] = xpos

        # for d in datasets:
        #     #add headers to table showing the rsquared and other aspects of the fit and dataset
        #     axarr[0,1].annotate(s=d_name[1:], xy=(xd[d],yaxis_pos[2]), fontsize=af, xycoords=xyc, alpha=0.75)

        axarr[0,1].annotate(s="{a} ({b})".format(a=method, b=doseunits), xy=(xaxis_left,yaxis_pos[3]), fontsize=af, xycoords=xyc, alpha=0.75)
        axarr[0,1].annotate(s="rsquared", xy=(xaxis_left,yaxis_pos[4]), fontsize=af, xycoords=xyc, alpha=0.75)
        axarr[0,1].annotate(s="hillslope", xy=(xaxis_left,yaxis_pos[5]), fontsize=af, xycoords=xyc, alpha=0.75)
        axarr[0,1].annotate(s="n_lowdose_datapoints", xy=(xaxis_left,yaxis_pos[6]), fontsize=af, xycoords=xyc, alpha=0.75)
        axarr[0,1].annotate(s="std_lowdose_datapoints", xy=(xaxis_left,yaxis_pos[7]), fontsize=af, xycoords=xyc, alpha=0.75)
        axarr[0,1].annotate(s="n_highdose_datapoints", xy=(xaxis_left,yaxis_pos[8]), fontsize=af, xycoords=xyc, alpha=0.75)
        axarr[0,1].annotate(s="std_highdose_datapoints", xy=(xaxis_left,yaxis_pos[9]), fontsize=af, xycoords=xyc, alpha=0.75)

        for d in datasets:
            # change the dataset name (e.g. "_orig" to "") to an empty string if there is only one dataset for analysis
            d_name = "" if len(datasets) == 1 else d
            #add headers to table showing the rsquared and other aspects of the fit and dataset
            axarr[0,1].annotate(s=d_name[1:], xy=(xd[d],yaxis_pos[2]), fontsize=af, xycoords=xyc, alpha=0.75)
            if dfe.loc["EC50_calculable{}".format(d),sLet]:
                EC50colour = "k" if dfe.loc["data_seems_okay{}".format(d),sLet] == True else "r"
                axarr[0,1].annotate(s="%0.2f" % dfe.loc["EC50{}".format(d),sLet], xy=(xd[d],yaxis_pos[3]), fontsize=af, xycoords=xyc, alpha=0.75, color=EC50colour)
                # rsquared of the fit to the data
                axarr[0,1].annotate(s="%0.2f"% dfe.loc["rsquared{}".format(d),sLet], xy=(xd[d],yaxis_pos[4]), fontsize=af, xycoords=xyc, alpha=0.75,color=dfe.loc["rsquared{}".format(d),"%s_colour" % sLet])
                # hillslope of the fit to the data
                axarr[0,1].annotate(s="%0.1f"% dfe.loc["hillslope{}".format(d),sLet], xy=(xd[d],yaxis_pos[5]), fontsize=af, xycoords=xyc, alpha=0.75, color = dfe.loc["hillslope{}".format(d),"%s_colour" % sLet])
                # number of lowdose datapoints
                axarr[0,1].annotate(s="%i"% dfe.loc["n_lowdose_datapoints{}".format(d),sLet], xy=(xd[d],yaxis_pos[6]), fontsize=af, xycoords=xyc, alpha=0.75, color = dfe.loc["n_lowdose_datapoints{}".format(d),"%s_colour" % sLet])
                # std of lowdose datapoints
                axarr[0,1].annotate(s="%0.2f"% dfe.loc["std_resp_lowdose_datapoints{}".format(d),sLet], xy=(xd[d],yaxis_pos[7]), fontsize=af, xycoords=xyc, alpha=0.75, color = dfe.loc["std_resp_lowdose_datapoints{}".format(d),"%s_colour" % sLet])
                 # number of highdose datapoints
                axarr[0,1].annotate(s="%i"% dfe.loc["n_highdose_datapoints{}".format(d),sLet], xy=(xd[d],yaxis_pos[8]), fontsize=af, xycoords=xyc, alpha=0.75, color = dfe.loc["n_highdose_datapoints{}".format(d),"%s_colour" % sLet])
                # std of highdose datapoints
                axarr[0,1].annotate(s="%0.2f"% dfe.loc["std_resp_highdose_datapoints{}".format(d),sLet], xy=(xd[d],yaxis_pos[9]), fontsize=af, xycoords=xyc, alpha=0.75, color = dfe.loc["std_resp_highdose_datapoints{}".format(d),"%s_colour" % sLet])

            else:
                if isinstance(dfe.loc["EC50{}".format(d),sLet], str):
                    dfe.loc["EC50_to_insert{}".format(d),sLet] = "N/A"
                elif isinstance(dfe.loc["EC50{}".format(d),sLet], float):
                    dfe.loc["EC50_to_insert{}".format(d),sLet] = "0.0f".format(dfe.loc["EC50{}".format(d),sLet])
                # insert error string or EC50, coloured red to indicate likely poor data
                axarr[0,1].annotate(s=dfe.loc["EC50_to_insert{}".format(d),sLet], xy=(xd[d],yaxis_pos[3]), fontsize=af, xycoords=xyc, alpha=0.75, color = "r")

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        #                                                                                                     #
        #         Dose-Respose Curve Fig04: More notes on automatic judgement of data quality                 #
        #                                                                                                     #
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        #            _________
        #           |    |    |          Subplot 4, axarr[1, 1]
        #           |____|____|
        #           |    |XXXX|
        #           |____|XXXX|

        Plot_Nr = 4

        # #add a table showing the rsquared and other aspects of the fit and dataset
        # axarr[1,1].annotate(s="ful", xy=(xful,yaxis_pos[2]), fontsize=af, xycoords=xyc, alpha=0.75)
        # axarr[1,1].annotate(s="orig", xy=(0.8,yaxis_pos[2]), fontsize=af, xycoords=xyc, alpha=0.75)

        axarr[1,1].annotate(s="dose conc. stepsize", xy=(xaxis_left,yaxis_pos[4]), fontsize=af, xycoords=xyc, alpha=0.75)
        axarr[1,1].annotate(s="slope at lowdose", xy=(xaxis_left,yaxis_pos[5]), fontsize=af, xycoords=xyc, alpha=0.75)
        axarr[1,1].annotate(s="slope at highdose", xy=(xaxis_left,yaxis_pos[6]), fontsize=af, xycoords=xyc, alpha=0.75)

        for n,d in enumerate(datasets):
            # change the dataset name (e.g. "_orig" to "") to an empty string if there is only one dataset for analysis
            d_name = "" if len(datasets) == 1 else d
            #add a table showing the rsquared and other aspects of the fit for adjusted datasets
            axarr[1,1].annotate(s=d_name[1:], xy=(xd[d],yaxis_pos[2]), fontsize=af, xycoords=xyc, alpha=0.75)
            EC50_calculable = dfe.loc["EC50_calculable{}".format(d),sLet]
            if EC50_calculable:
                # add the stepsize near the EC50, which determines whether more dose concentrations are necessary
                stepsize = "{:0.2f}".format(dfe.loc["doseconc_stepsize_at_EC50{}".format(d), sLet])
                stepcolour = dfe.loc["doseconc_stepsize_at_EC50{}".format(d), "%s_colour" % sLet]

                # if dfe.loc["doseconc_stepsize_at_EC50{}".format(d),"%s_okay" % sLet] == True:
                axarr[1,1].annotate(s=stepsize,xy=(xd[d],yaxis_pos[4]), fontsize=af, xycoords=xyc, alpha=0.75, color=stepcolour)
                saxe_lowdose = "{:0.2f}".format(dfe.loc["saxe_lowdose{}".format(d), sLet])
                saxe_highdose = "{:0.2f}".format(dfe.loc["saxe_highdose{}".format(d), sLet])
                axarr[1,1].annotate(s=saxe_lowdose,xy=(xd[d],yaxis_pos[5]), fontsize=af, xycoords=xyc, alpha=0.75,
                                    color=dfe.loc["saxe_lowdose{}".format(d),"%s_colour" % sLet])
                axarr[1,1].annotate(s=saxe_highdose,xy=(xd[d],yaxis_pos[6]), fontsize=af, xycoords=xyc, alpha=0.75,
                                    color=dfe.loc["saxe_highdose{}".format(d),"%s_colour" % sLet])
                if stepcolour != "k":
                    doseconc_steps_at_EC50 = dfe.loc["doseconc_steps_at_EC50{}".format(d),sLet]
                    axarr[1,0].plot(doseconc_steps_at_EC50,(0,0), color=stepcolour, linestyle="-", lw=2)

                # if the slope at the lowdose is above the chosen cutoff, draw a line on the normalised plot, axarr[1,0]
                if dfe.loc["saxe_lowdose{}".format(d), sLet] > settings["max_lowdose_slope"]:
                    saxe_lowdose_values = dfe.loc["saxe_lowdose_values{}".format(d), sLet]
                    # draw red vertical line showing the slope at the lowdose datapoint
                    axarr[1,0].plot(saxe_lowdose_values[0], saxe_lowdose_values[1], 'r-', lw=2)
                # if the slope at the highdose is higher than the chosen cutoff, draw a line on tho normalised plot
                if dfe.loc["saxe_highdose{}".format(d), sLet] > settings["max_highdose_slope"]:
                    saxe_highdose_values = dfe.loc["saxe_highdose_values{}".format(d), sLet]
                    # draw red vertical line showing the slope at the lowdose datapoint
                    axarr[1,0].plot(saxe_highdose_values[0], saxe_highdose_values[1], 'r-', lw=2)

            else:
                # optional: print the dataframe showing which parameters are not acceptable
                EC50 = dfe.loc["EC50{}".format(d),sLet]
                if isinstance(EC50, str):
                    # define x_position of annotation. Place the orig at around 0.6, & 2nd dataset(_ful) at around 0.3
                    xd_wide = xd[d]-0.3-0.3*n
                    axarr[1,1].annotate(s=d_name[1:],xy=(xd_wide,yaxis_pos[5]), fontsize=af, xycoords=xyc, color="r")
                    axarr[1,1].annotate(s=EC50,xy=(xd_wide,yaxis_pos[6]), fontsize=af, xycoords=xyc, color="r")

        # create a dictionary with the formatted EC50 values for printing on the figure
        EC50_str_dict = {}
        for d in datasets:
            EC50 = dfe.loc["EC50{}".format(d),sLet]
            EC50_str_dict[d] =  "%.01f" % EC50 if isinstance(EC50, float) else EC50
        # use % formatting to add the padding for the sample name
        padding = sample_name_len_max if sample_name_len_max <= 44 else 44
        samplestring = "\n{sLet} {s:>%s} {m} : " % int(padding)
        sys.stdout.write(samplestring.format(sLet=sLet, s=sample_name[:44], m=method))
        for n, d in enumerate(datasets):
            data_seems_okay = dfe.loc["data_seems_okay{}".format(d),sLet]
            if data_seems_okay:
                data_evaluation = "(data seems good)"
            else:
                if EC50_calculable:
                    data_evaluation = "(needs checking)"
                else:
                    if EC50_str_dict[d] == "insuff_lowresp_dp":
                        data_evaluation = "(EC50 not calculable, insufficient low response datapoints)"
                    elif EC50_str_dict[d] == "insuff_highresp_dp":
                        data_evaluation = "(EC50 not calculable, insufficient high response datapoints)"
            sys.stdout.write("{d:>4} = {EC:>5} {eval:>17}".format(d=d_name[1:], EC=EC50_str_dict[d], eval=data_evaluation))
            if n < len(datasets):
                sys.stdout.write(", ")
            sys.stdout.flush()

        #save figure with the fitted curve and calculated EC50 value
        fig.tight_layout()
        fig.savefig(fig0_single_sample_png, format='png', dpi=140)
        if settings["save_as_pdf"] in (True,"TRUE"):
            fig.savefig(fig0_single_sample_pdf, format='pdf')
        plt.close('all')
        dict_dfe[sLet] = dfe

    # create a DataFrame for EVALuation (df_eval), which combines all the dfe dataframes from the previous loop
    # note that df_eval and dfe are transposed, which is somewhat confusing.
    df_eval = pd.DataFrame()
    for sLet in dict_dfe:
        # extract each individual dfe separately, and add to df_eval
        indiv_dfe = dict_dfe[sLet]
        df_eval = pd.concat([df_eval,indiv_dfe], axis=1)

    # transpose, (exchange column and rows)
    df_eval = df_eval.T
    # create a column containing the sample letter and sample name combined
    df_eval['sLet_plus_sample_name'] = df_eval['sLet'] + " " + df_eval['sample_name']
    # create columns to be visible first in the output csv and excel files
    sel_col = ['sLet_plus_sample_name']
    for d in datasets:
        sel_col_d = ['EC50{}'.format(d), 'data_seems_okay{}'.format(d), 'rsquared{}'.format(d), 'doseconc_stepsize_at_EC50{}'.format(d), 'EC50_hill_eq{}'.format(d), 'curve_min_norm{}'.format(d),
                     'hillslope{}'.format(d), 'std_resp_highdose_datapoints{}'.format(d), 'n_highdose_datapoints{}'.format(d)]
        sel_col += sel_col_d
    # sel_col += ['sNum', 'n_highdose_datapoints_orig','n_lowdose_datapoints_orig']
    # reindex so columns are first
    df_eval = tools.reindex_df_so_selected_cols_are_first(df_eval, sel_col, raise_error=False)

    # divide the dataframe into two new dataframes, one for the values (_val) and another for the booleans
    # related to whether the data is okay or not. simply copy orig dataframe and drop unwanted columns.
    df_eval_values = df_eval.copy()
    for row in df_eval_values.index:
        if "_okay" in row:
            df_eval_values.drop(row, inplace=True)
        if "_colour" in row:
            df_eval_values.drop(row, inplace=True)
    # copy orig dataframe and drop unwanted columns.
    df_eval_bool = df_eval.copy()
    for row in df_eval_bool.index:
        if "_okay" not in row:
            df_eval_bool.drop(row, inplace=True)
    # drop empty columns
    df_eval_bool.dropna(how="all", axis=1, inplace=True)

    # sort by index (sLet)
    df_eval_values.sort_index(inplace=True)
    # # reindex so that selected columns are visible first
    # df_eval_values_selected_cols = ['sample_name', 'EC50_ful', 'data_seems_okay_ful', 'data_seems_okay_orig', 'rsquared_orig', 'EC50_hill_eq_ful','n_highdose_datapoints_orig','n_lowdose_datapoints_orig','sNum']
    # df_eval_values = reindex_df_so_selected_cols_are_first(df_eval_values, df_eval_values_selected_cols)
    # sort by index (sLet)
    df_eval_bool.sort_index(inplace=True)
    # # reindex so that selected columns are visible first
    # df_eval_bool = reindex_df_so_selected_cols_are_first(df_eval_bool, df_eval_values_selected_cols)

    #set up the summary figure to contain 2 subplots
    n_plots_per_fig = 2
    nrows = 2
    dict_organising_subplots = tools.create_dict_organising_subplots(n_plots_per_fig,n_rows=nrows)
    #set the fontsize for the figure
    fig_fontsize = 6

    # for d in datasets:
    #     # go through all of the data, and set the label colour to red if the data_seems_okay is false
    #     df_eval_values["xlabel_colour{}".format(d)] = df_eval_values["data_seems_okay{}".format(d)].apply(lambda c: "k" if c == True else "r")

    for d in datasets:
        # change the dataset name (e.g. "_orig" to "") to an empty string if there is only one dataset for analysis
        d_name = "" if len(datasets) == 1 else d
        # go through all of the data, and set the label colour to red if the data_seems_okay is false
        df_eval_values["xlabel_colour{}".format(d)] = df_eval_values["data_seems_okay{}".format(d)].apply(lambda c: "k" if c == True else "r")

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        #                                                                                                     #
        #             Summ_Plot01, summ_Fig00: scattergram, original data with fitted curve                   #
        #                                                                                                     #
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        #            _________
        #           |XXXXXXXXX|          Subplot 0, axarr[0]
        #           |XXXXXXXXX|
        #           |         |
        #           |_________|

        # set the axes array position
        r = 0
        # set the figure number (canvas number, for summmary figures)
        fig_nr = 1
        fig, axarr = plt.subplots(nrows=2, ncols=1, dpi=300)
        EC50_calculable_list = list(df_eval_values.loc[:,"EC50_calculable{}".format(d)])

        if True in EC50_calculable_list:
            # plot the curves first, as this is used in the legend
            for sLet in df_eval_values.index:
                sNum = df_eval_values.loc[sLet,'sNum']
                axarr[r].plot(df_eval_values.loc[sLet,"x_fitted{}".format(d)],
                                           df_eval_values.loc[sLet,"y_fitted{}".format(d)],
                                           '-',
                                           color = t20[sNum],
                                           alpha = 0.8,
                                           label = sLet)
            # set the legend. Note that this is based on the last samples plotted
            lg = axarr[r].legend(df_eval_values['sLet_plus_sample_name'], loc='upper right', scatterpoints=1)
            lg.draw_frame(False)
        # set xlabel, ylabel, title, grid, etc
        axarr[r].set_xlabel("{a} ({b})".format(a=doselabel, b=doseunits), fontsize = fig_fontsize)
        axarr[r].set_ylabel(settings["y-axis (response) label"],rotation='vertical', fontsize = fig_fontsize)
        axarr[r].set_title('{a} data                              {b}'.format(a=d_name[1:],b=data_file),
                                fontsize = fig_fontsize, x = 0.22)
        axarr[r].grid(True, color = '0.75')
        # plot the raw datapoints last, as these interfere with the legend
        ylim_min_list = []
        ylim_max_list = []
        for sLet in df_eval_values.index:
            sNum = df_eval_values.loc[sLet,'sNum']
            axarr[r].scatter(df_eval_values.loc[sLet,"x{}".format(d)],
                                df_eval_values.loc[sLet,"y{}".format(d)],
                                color = t20[sNum],
                                alpha = 0.8,
                                s = 15,
                                label = sLet)
            # grab the min and max for the y-values. This is used for ax.set_ylim.
            ylim_min_list.append(np.array(df_eval_values.loc[sLet,"y{}".format(d)]).min())
            ylim_max_list.append(np.array(df_eval_values.loc[sLet, "y{}".format(d)]).max())

        lowest_y_datapoint = np.array(ylim_min_list).min()
        highest_y_datapoint = np.array(ylim_max_list).max()
        axarr[r].set_ylim(lowest_y_datapoint, highest_y_datapoint)

        if not True in list(df_eval_values.loc[:,"EC50_calculable{}".format(d)]):
            # if none of the curves could be fitted, base the legend on the datapoints rather than the curves
            lg = axarr[r].legend(df_eval_values['sLet_plus_sample_name'], loc='upper right', scatterpoints=1)
            lg.draw_frame(False)
        # set the x-axis limit so that the legend does not hide too many data points
        # find the maximum dose concentration in the whole experiment for that day
        maxAC = x_orig.max().max()
        # # obtain the variable altering the extension of the x-axis
        # x_axis_extension_after_dosemax_in_summ_plot = dff.loc[fn, "x-axis extension in summary fig_0"]
        # #define the limit of the x-axis as the maximum dose conc.
        # xlim_max = maxAC + x_axis_extension_after_dosemax_in_summ_plot
        # # set the x-axis limits
        # axarr[r].set_xlim(-10,xlim_max)

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        #                                                                                                     #
        #                      Summ_Plot02, summ_Fig01: barchart EC50 (sLet on x-axis)                        #
        #                                                                                                     #
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        #            _________
        #           |         |          Subplot 1, axarr[1]
        #           |_________|
        #           |XXXXXXXXX|
        #           |XXXXXXXXX|

        # Plot_Nr = Plot_Nr + 1
        # newfig, savefig, fig_nr, plot_nr_in_fig, r, c = dict_organising_subplots[Plot_Nr]
        # set the axes array position
        r = 1
        # use the sample letter alone as the name on the x-axis
        x_names = df_eval_values.sLet
        # the number of boxes in the bar-chart is the length of the initial dataset
        x_n_boxes = x_names.shape[0]
        # the position of the boxes in the bar-chart is the range(n_boxes)
        box_indices = range(x_n_boxes)
        if True in list(df_eval_values.loc[:,"EC50_calculable{}".format(d)]):
            # define the y-axis data
            y_EC50 = pd.to_numeric(df_eval_values["EC50{}".format(d)], errors="coerce")
            # obtain yerr [[REMOVED: NOT SUITABLE FOR DATASETS WITH SMALL DOSE VALUES]]
            #yerr = df_eval_values["yerr{}".format(d)].fillna(0)
            # create a new object (usually not used) that contains a bar-chart on figure (ax) #1
            bar_container = axarr[r].bar(box_indices, y_EC50, color = t20, align="center",
                                         error_kw=dict(ecolor='k', lw=1, capsize=2, capthick=1))#yerr=yerr,
        # set the xticks
        axarr[r].set_xticks(box_indices)
        # set the labels of the x-axis
        axarr[r].set_xticklabels(x_names)
        for xtick, colour in zip(axarr[r].get_xticklabels(), df_eval_values["xlabel_colour{}".format(d)]):
            xtick.set_color(colour)
        # set the limits of the x-axis
        axarr[r].set_xlim([-1, x_n_boxes])
        # set the x-axis title
        axarr[r].set_xlabel("sample letter")
        # set the y-axis title
        axarr[r].set_ylabel("{a} ({b})".format(a=method, b=doseunits))
        #save figure
        fig.tight_layout()
        fig.savefig(dff.loc[fn,"EC50_analysis_fig_basename"] + d_name + '.png', format='png', dpi=150)
        if settings["save_as_pdf"] in (True,"TRUE"):
            fig.savefig(dff.loc[fn,"EC50_analysis_fig_basename_pdf"] + d_name + '.pdf', format='pdf')

        if True in list(df_eval_values.loc[:,"EC50_calculable{}".format(d)]):
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            #                                                                                                     #
            #                             Barchart01: barchart EC50 (full name on x-axis)                         #
            #                                                                                                     #
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            #            _________
            #           |XXXXXXXXX|          Full Canvas
            #           |XXXXXXXXX|
            #           |XXXXXXXXX|
            #           |XXXXXXXXX|
            '''
            Barchart01: barchart original data, EC50_orig (full name on x-axis)
            '''
            bar_fig_nr = 1
            # create new figure (i.e., canvas)
            fig, ax = plt.subplots()
            # use the sample letter plus the sample name as the name on the x-axis
            x_names = df_eval_values.sLet_plus_sample_name
            # yerr_ful = df_eval_values.residuals_mean_ful.fillna(0)*100
            axarr[r].set_title('{a} data                              {b}'.format(a=d_name[1:],b=data_file),
                                fontsize = fig_fontsize, x = 0.22)
            # create a new object (usually not used) that contains a bar-chart on figure (ax) #1
            bar_container = ax.bar(box_indices, y_EC50, color=t20, align="center",
                                   error_kw=dict(ecolor='k', lw=1, capsize=2, capthick=1))#yerr=yerr,
            # set the xticks. apply the appropriate colour, as decided by the "judge" script
            ax.set_xticks(box_indices)
            # df_eval_values["xlabel_colour_orig"] = df_eval_values.data_seems_okay_orig.apply(lambda c: "k" if c == True else "r")
            for xtick, colour in zip(ax.get_xticklabels(), df_eval_values["xlabel_colour{}".format(d)]):
                xtick.set_color(colour)
            # set the labels of the x-axis
            ax.set_xticklabels(x_names, rotation = 90)
            # set the limits of the x-axis
            ax.set_xlim([-1, x_n_boxes])
            # set the y-axis title
            ax.set_ylabel("{a} ({b})".format(a=method, b=doseunits))

            #save figure
            try:
                fig.tight_layout()
            except ValueError:
                sys.stdout.write("Sample names may need to be truncated for barchart. Current length = {}".format(x_names.str.len().max()))
                sys.stdout.flush()
                x_names = x_names.str[0:70]
                ax.set_xticklabels(x_names, rotation=90)
                fig.tight_layout()
            fig.savefig(dff.loc[fn,"EC50_analysis_fig_basename"] + d_name + "_bar" + '.png', format='png', dpi=150)
            if settings["save_as_pdf"] in (True, "TRUE"):
                fig.savefig(dff.loc[fn,"EC50_analysis_fig_basename_pdf"] + d_name + "_bar" + '.pdf', format='pdf')
            plt.close('all')

    # drop the columns with a large number of datapoints, to reduce size of the output files
    #list_fitted_cols = "y_fitted_norm_ful", "x_fitted_norm","y_fitted_norm_orig","x_fitted_orig","y_fitted_orig", "y_fitted_ful",

    list_arraylike_cols = []
    list_fitted_cols = []
    for d in datasets:
        # define list of large columns with fitted data, for that dataset
        fitted_d = ["x_fitted_norm{}".format(d), "y_fitted_norm{}".format(d), "x_fitted{}".format(d), "y_fitted{}".format(d)]
        # add to full list for all datasets
        list_fitted_cols = list_fitted_cols + fitted_d
        # define list of arraylike columns, for that dataset
        arraylike_d = ["x{}".format(d), "y{}".format(d), "indices_lowdose_datapoints{}".format(d),
                       "indices_lowdose_datapoints_excl_nearest_EC50{}".format(d),"response_lowdose_datapoints{}".format(d)]
        # add to full list for all datasets
        list_arraylike_cols = list_arraylike_cols + arraylike_d

    # drop columns with fitted data, to reduce filesize
    df_eval_values.drop(list_fitted_cols, axis=1, inplace=True)

    # convert listlike to stringlists
    df_eval_values = tools.convert_listlike_cols_to_str(df_eval_values, list_arraylike_cols)

    # df_eval_bool = convert_listlike_cols_to_str(df_eval_bool, list_arraylike_cols)
    # save evaluation dataframe to csv
    df_eval_values.to_csv(dff.loc[fn,"ofd_EC50_eval_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC)
    df_eval_values.to_csv(dff.loc[fn,"ofd_EC50_eval_tabsep_csv"], sep="\t", quoting=csv.QUOTE_NONNUMERIC)
    # save evaluation dataframe to excel
    writer = pd.ExcelWriter(dff.loc[fn,"ofd_EC50_eval_excel"])#engine='xlsxwriter'
    df_eval_values.to_excel(writer, sheet_name="v_" + data_file[:20])
    df_eval_bool.to_excel(writer, sheet_name="b_" + data_file[:20])
    settings.to_frame().to_excel(writer, sheet_name="settings")
    writer.save()
    writer.close()
    # save the settings in the csv folder, so that a permanent record of the settings is kept
    settings_csv_path = os.path.join(dff.loc[fn, "ofd_csv"], "settings.csv")
    settings.to_csv(settings_csv_path)
    print('\n-------------------------------------\n')

def calc_EC50_brent_eq(sLet, sample_name, hill_constants, y50_norm):
    """ Calculates the EC50 from the fitted curve to normalised data.

    Parameters
    ----------
        sLet : string
            Sample letter
        sample_name : string
            Full sample name
        hill_constants: tuple
            Hill constants calculated from the fitting algorithm
        y50_norm : float
            The y-value at 50% response, used for the EC50 calculation. Normalised between 0 and 1, as with the rest
            of the datapoints.

    Returns
    -------
        EC50_norm_bq : float
            EC50 value derived from normalised data. Will always be a float between 0.0 and 1.0.
            Note that EC50_norm_bq will be later "denormalised" to EC50{}".format(dataset), e.g. EC50_orig.
        EC50_calculable : bool
            Whether EC50 was calculable using the brentq function.

    Notes
    -------
    The brent equation(in python the brentq function) is used here to find the position on the x-axis at a
    particular y-value. This tends to be more robust than the EC50 parameter derived from the Hill equation.

    In some cases, the brentq function will result in ValueError: f(a) and f(b) must have different signs.
    This is due to a value outside of the initial range (between 0 and 1), and usually indicates a very poor fit
    http://stackoverflow.com/questions/22277982/how-to-find-50-point-after-curve-fitting-using-numpy/22288826#22288826
    """

    # try to determine the EC50 between 0.0 and 1.0, which works for most data
    try:
        EC50_norm_bq = brentq(tools.hill_eq_brentq, 0.0, 1.0,args = (hill_constants, y50_norm))
        EC50_calculable = True
    except ValueError:
        # widen the scope of the EC50 to outside of the actual range of datapoints.
        # This generally indicates poor data quality, however there may be exceptions. The need for a wider range
        # is therefore printed in the console, however this wider scope is not a feature in the judge_fit algorithm,
        # and therefore is not a factor used to determine EC50 data quality.
        try:
            print("ValueError encountered in brentq. Attempting wider scope for EC50. "
                  "This is unusual, and may indicate poor data quality.")
            EC50_norm_bq = brentq(tools.hill_eq_brentq, -1.0, 2.0,args = (hill_constants, y50_norm))
            EC50_calculable = True
        except ValueError:
            EC50_calculable = False
            print("ValueError encountered in brentq, even with wider scope! EC50 is not calculable. "
                  "Sample : %s_%s" % (sLet, sample_name))
    return EC50_norm_bq, EC50_calculable

def standardise_doseconc_data(fn, dff, df_dose_orig, df_resp_all, data_file_path):
    """ Standardises the dose and response data to a uniform pandas dataframes.

    For 96-well formats, the df_dose_all is created, and df_resp_all is simply returned.
    For 12-well formats where dose and response data both come from a single excel file, both df_dose_all and
    df_resp_all are created.

    Parameters
    ----------
    fn : int
        File number in list of files to analyse.
    assay_type : string
        Assay type, for example indicating a 96-well assay, with 12 dose concentrations and 8 samples
    df_dose_orig : pandas DataFrame
        The original dataframe with the dose data for standardisation. It will, for example, show a horizontal
        increase in dose for 96-well assay, with 12 dose concentrations and 8 samples, but a vertical increase
        in dose for 8 dose concentrations and 12 samples.
    df_resp_all : pandas DataFrame
        Dataframe with the response data, which is usually more uniform than the dose data and requires less
        standardisation. For 96-well formats, this is already created and will simply be returned. For 12-well
        formats containing dose and response data in excel format, this object will be created.
    data_file_path : filepath
        Path to the file containing the response data. Note that some input files will contain both dose and
        response data (e.g. Excel inputs). Others only contain response data, with doses provided in a separate file.

    Returns
    -------
    df_dose_all : pandas DataFrame
        Dose concentrations in a standardised format, matching the response dataframe.
    df_resp_all : pandas DataFrame
        Response concentrations in a standardised format, matching the dose dataframe.

    Notes
    -----
    Currently compatible with the 96-well data from the versamax microplate reader (8, 12 and 24 doseconc
    variations) and also the data from the Fluostar, which an read 12-well plates directly.  For the 12-well
    format, the excel file contains both dose and response data. For 96-well formats, a separate excel file contains
    the dose data.
    """
    eccpy, resp_datafileformat, resp_machinetype, resp_assaytype = dff.loc[fn, "response dataformat"].split("|")

    # create a list of double-letters (AA, AB, AC, AD etc) to act as the index for the samples (sample letters, sLet)
    let_list = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
    sLet_list = []
    for let1 in let_list[0:8]:
        for let2 in let_list:
            joined_letters = "{a}{b}".format(a=let1, b=let2)
            sLet_list.append(joined_letters)

    #######################################################################################################
    #                                                                                                     #
    #                          Formatting Dose Concentrations for xlsx|generic                            #
    #                                                                                                     #
    #######################################################################################################
    if resp_machinetype == "generic":
            #load the data into two dataframes, 1) dose concentrations, 2) response data
            df_dose_data = pd.read_excel(data_file_path, sheetname='dose')
            df_resp_data = pd.read_excel(data_file_path, sheetname='response')

            # drop any columns that are completely empty
            df_dose_data = df_dose_data.dropna(how="all", axis=1)
            df_resp_data = df_resp_data.dropna(how="all", axis = 1)

            # drop any rows that are completely empty
            df_dose_data = df_dose_data.dropna(how="all")
            df_resp_data = df_resp_data.dropna(how="all")

            # currently read_excel doesn't allow duplicate columns, adding .1, .2, .3 at the end of the column names
            # see https://github.com/pydata/pandas/issues/10523
            # in this case, if the secondlast character of the column name (sample name) is ".", remove the .1, .2 etc
            # this allows the same sample to be analysed multiple times in the same experiment
            # but will give errors if someone actually has a period (".") in the sample name
            remove_dup_numbers = lambda x : ".".join(x.split(".")[:-1]) if "." in x[-3:] else x
            df_resp_data.columns = pd.Series(df_resp_data.columns).apply(remove_dup_numbers)
            df_dose_data.columns = pd.Series(df_dose_data.columns).apply(remove_dup_numbers)

            # check if the column names are exactly the same
            # People also often have names that are inconsistent between the dose and the response tabs.
            if list(df_resp_data.columns) != list(df_dose_data.columns):
                dose_cols = list(df_dose_data.columns)
                response_cols = list(df_resp_data.columns)
                diff = list(set(dose_cols) - set(response_cols))
                print("The sample names are different between the dose and response tabs. "
                      "Strongly recommend double-checking the names, and deleting any remaining "
                      "junk in empty columns. The response names will be taken as TRUE, and dose names ignored. "
                      "    dose samples: {d}\nresponse samples: {r} \ndifferences: {dr}".format(d=dose_cols,
                                                                                                r=response_cols,
                                                                                                dr=diff))
                # only a warning is shown. Uncomment the next line if you want an error, which stops the program.
                #raise DataMismatchError()
            # temp_index = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ123456789")
            # The excel files often contain junk, with undeleted data from previous experiments.
            # People also often have names that are inconsistent between the dose and the response tabs.
            # Assume the df_resp_data has the "correct" sample names
            # count the number of dose concentrations (according to the number of rows in the response tab)
            n_doseconc = df_resp_data.shape[0]
            # count the number of samples (according to the A600 tab)
            n_samples = df_resp_data.shape[1]
            # slice the sample data
            df_resp_all = df_resp_data.iloc[:n_doseconc + 1,:n_samples + 1]
            df_dose_all = df_dose_data.iloc[:n_doseconc + 1,:n_samples + 1]
            # transpose (exchange columns and rows) to standardise format with 96-well data
            df_resp_all = df_resp_all.T
            df_dose_all = df_dose_all.T
            # assume response data has the correct names, and label the columns for dose and resp dataframes likewise
            df_resp_all['samples'] = df_resp_data.columns
            # label all columns (after the dropna step) as containing data
            df_resp_all['Contains_Data'] = True
            # label the index with letters (sLet)
            index_letters = sLet_list[0:n_samples]
            df_resp_all.index = index_letters
            df_dose_all.index = index_letters
            # add the samples columnn to the dose conc. dataframe, based on the "df_samples_data" dataframe.
            # Ignore any columns in the samples dataframe longer than the doseconc dataframe.
            df_dose_all['samples'] = df_resp_data.columns
            # label all columns as containing data (in the generic template, only useful data should be pasted)
            df_dose_all['Contains_Data'] = True

    #######################################################################################################
    #                                                                                                     #
    #          Formatting Dose Concentrations for txt|versamax  (96-well plate reader)                    #
    #                                                                                                     #
    #######################################################################################################
    elif resp_machinetype == "versamax":
        if resp_assaytype == "8dose12sample":
            # transpose (flip index and column), and create a copy of dataframe
            df_dose_orig = df_dose_orig.T
            # number of plates in excel sheet
            num_plates = 4
            # number of columns per plate (from sLet to Contains_Data)
            cols_per_plate = 10
            # create empty dataframe to hold info from all plates
            df_dose_all = pd.DataFrame()
            for i in range(num_plates):
                # isolate the dose_conc data for a single plate. Create copy.
                dfDC_plate = df_dose_orig.iloc[:,cols_per_plate * i:cols_per_plate * (i + 1)].copy()
                # add data to existing dataframe
                df_dose_all = pd.concat([df_dose_all,dfDC_plate])
            # move sLet to the index
            df_dose_all.set_index("sLet", inplace=True)
            # reformat the columns from ints (1,2) to formatted 2-digit strings (01,02)
            df_dose_all.columns = tools.format_cols_2digit(df_dose_all)
            # change the "Contains_Data" column to bool, rather than object
            df_dose_all['Contains_Data'] = df_dose_all.Contains_Data.apply(lambda xx: True if xx == "TRUE" else False)

        if resp_assaytype == "8_ampconc":
            # DEPRECATED! USE "8dose12sample" INSTEAD FOR FUTURE DATA
            # transpose (flip index and column), and create a copy of dataframe
            df_dose_orig = df_dose_orig.T
            # isolate data from plate 1
            dfDC_plate1 = df_dose_orig.iloc[:,:9].copy()
            # isolate data from plate 2
            dfDC_plate2 = df_dose_orig.iloc[:,9:].copy()
            # join dataframes containing AmpConc for plates 1 and 2
            df_dose_all = pd.concat([dfDC_plate1,dfDC_plate2])
            # convert the index to the sample letters, sLet
            # NOTE in this deprecated code limited to only 2 plates, only single letters are used.
            df_dose_all.index = let_list[0: df_dose_all.shape[0]]
            # now convert to simple columns
            # create numbers for new index (simple A01 style)
            DC_numbers_8DC = ["%02d" % n for n in range(1,9)]
            DC_numbers_8DC.append('Contains_Data')
            #DC_numbers_8DC.append('orig_index')
            df_dose_all.columns = DC_numbers_8DC
            # because the columns in the orig excel file contain a mix of int and bool, the dtype is "object"
            # change the "Contains_Data" column to bool, rather than object
            df_dose_all['Contains_Data'] = df_dose_all.Contains_Data.apply(lambda xx: True if xx == "TRUE" else False)

        elif resp_assaytype == "12dose8sample" or resp_assaytype == "12_ampconc":
            # 12_ampconc is DEPRECATED. USE "12dose8sample" INSTEAD FOR FUTURE DATA
            # reformat the columns from ints (1,2) to formatted 2-digit strings (01,02)
            df_dose_orig.columns = tools.format_cols_2digit(df_dose_orig)
            # create copy of orig dataframe. The df_dose_all should now be formatted to match the df_resp_all.
            df_dose_all = df_dose_orig.copy()

        elif resp_assaytype == "24dose4sample":
            # reformat the columns from ints (1,2) to formatted 2-digit strings (01,02)
            df_dose_orig.columns = tools.format_cols_2digit(df_dose_orig)
            # select indices for even and odd rows (the 24 dose conc cover 2 rows, starting at A01, and ending at B12)
            even_rows = df_dose_orig.index[0::2]
            odd_rows = df_dose_orig.index[1::2]
            # select data of even and odd rows
            dfDC_even_rows = df_dose_orig.loc[even_rows].copy()
            dfDC_odd_rows = df_dose_orig.loc[odd_rows].copy()
            # create a copy of the "Contains_Data" column from the even rows, and then drop from dataframe
            series_Contains_Data = dfDC_even_rows.Contains_Data.copy()
            dfDC_even_rows.drop("Contains_Data", axis=1, inplace=True)
            # copy the index from the dataframe with even rows to the dataframe with odd rows
            dfDC_odd_rows.index = dfDC_even_rows.index
            # reformat the columns from ints (1,2) to formatted 2-digit strings (01,02)
            dfDC_odd_rows.columns = tools.format_cols_2digit(dfDC_odd_rows)
            # create df_dose_all by concatenating dataframes along axis 1 (joining columns)
            df_dose_all = pd.concat([dfDC_even_rows, dfDC_odd_rows], axis=1)
            df_dose_all["Contains_Data"] = series_Contains_Data

        elif resp_assaytype == "24_ampconc":
            # DEPRECATED! USE "24dose4sample" INSTEAD FOR FUTURE DATA
            # reformat the columns from ints (1,2) to formatted 2-digit strings (01,02)
            df_dose_orig.columns = tools.format_cols_2digit(df_dose_orig)
            temp_index = list("ABCDEFGHIJKLMNOP")
            df_dose_orig.index = temp_index
            rows_containing_01to12_AC = df_dose_orig.index[0::2]
            rows_containing_13to24_AC = df_dose_orig.index[1::2]
            final_index = list("ABCDEFGH")
            dfDC_01to12 = df_dose_orig.loc[rows_containing_01to12_AC]
            dfDC_01to12.index = final_index
            dfDC_01to12 = dfDC_01to12.drop("Contains_Data", axis=1)
            dfDC_13to24 = df_dose_orig.loc[rows_containing_13to24_AC]
            dfDC_13to24.index = final_index
            dfDC_13to24_cols = ["%02d" % r for r in range(13,25)]
            dfDC_13to24_cols.append(df_dose_orig.columns[-1])
            dfDC_13to24.columns = dfDC_13to24_cols
            df_dose_all = pd.concat([dfDC_01to12, dfDC_13to24], axis=1)
        else:
            tools.DatafileError('response assaytype not identified')

    #######################################################################################################
    #                                                                                                     #
    #         Formatting Dose Concentrations for xlsx|fluostar  (12-well & 96-well plate reader)          #
    #                                                                                                     #
    #######################################################################################################
    elif resp_machinetype == "fluostar":
        if resp_assaytype == "12wellplate":
            #load the data into two dataframes, 1) dose concentrations, 2) response data
            df_resp = pd.read_excel(data_file_path, sheetname='A600', skiprows = 2)
            df_amp_conc = pd.read_excel(data_file_path, sheetname='amp_conc', skiprows = 2)
            number_of_initial_A600_columns = df_resp.shape[1]
            #select only the useful data, assuming that the first 4 columns are ignored
            df_resp_data = df_resp.iloc[:,3:]
            # select useful data. ignore any columns longer than the A600 data
            df_amp_conc_data = df_amp_conc.iloc[:,3:number_of_initial_A600_columns]
            # drop any columns that are completely empty
            df_resp_data = df_resp_data.dropna(how="all", axis = 1)
            df_amp_conc_data = df_amp_conc_data.dropna(how="all", axis=1)
            # drop any rows that are completely empty
            df_resp_data = df_resp_data.dropna(how="all")
            df_amp_conc_data = df_amp_conc_data.dropna(how="all")
            # currently read_excel doesn't allow duplicate columns, adding .1, .2, .3 at the end of the column names
            # see https://github.com/pydata/pandas/issues/10523
            # in this case, if the secondlast character of the column name (sample name) is ".", remove the .1, .2 etc
            # this allows the same sample to be analysed multiple times in the same experiment
            # but will give errors if someone actually has a period (".") in the sample name
            remove_dup_numbers = lambda x : ".".join(x.split(".")[:-1]) if "." in x[-3:] else x
            df_resp_data.columns = pd.Series(df_resp_data.columns).apply(remove_dup_numbers)
            df_amp_conc_data.columns = pd.Series(df_amp_conc_data.columns).apply(remove_dup_numbers)
            # check if the column names are exactly the same
            if list(df_resp_data.columns) != list(df_amp_conc_data.columns):
                print("The sample names are different between the A600 and amp_conc tabs. "
                      "The 'samples' tab in the excel file is no longer used in the script, but is still "
                      "useful for generating accurate and consistent sample names. The names should be the same for both the "
                      "A600 and amp_conc tabs. Strongly recommend double-checking the names, and deleting any remaining "
                      "junk in empty columns. The OD600 names will be taken as TRUE, and amp_conc names ignored. "
                      "A600\n%sAmpconc\n%s " % (list(df_resp_data.columns), list(df_amp_conc_data.columns)))
            # df_samples_data = df_samples.iloc[:,0:].dropna(how="all", axis=1)
            temp_index = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ123456789")
            # The excel files often contain junk, with undeleted data from previous experiments
            # Assume the df_resp_data has the "correct" sample names
            # count the number of dose concentrations (according to the number of rows in the A600 tab)
            n_doseconc = df_resp_data.shape[0]
            # count the number of samples (according to the A600 tab)
            n_samples = df_resp_data.shape[1]
            # slice the sample data, transpose (exchange columns and rows) to standardise format weth 96-well data
            df_resp_all = df_resp_data.iloc[:n_doseconc + 1,:n_samples + 1].T
            # assume the A600 data has the correct names, and label the columns likewise
            df_resp_all['samples'] = df_resp_data.columns
            # label all columns (after the dropna step) as containing data
            df_resp_all['Contains_Data'] = True
            # label the index with letters (sLet)
            index_letters = temp_index[0:n_samples]
            df_resp_all.index = index_letters
            # slice the sample data, transpose (exchange columns and rows) to standardise format weth 96-well data
            df_dose_all = df_amp_conc_data.iloc[:n_doseconc + 1,:n_samples + 1].T
            # add the samples columnn to the dose conc. dataframe, based on the "df_samples_data" dataframe.
            # Ignore any columns in the samples dataframe longer than the doseconc dataframe.
            df_dose_all['samples'] = df_resp_data.columns
            # label all columns (after the dropna step) as containing data
            df_dose_all['Contains_Data'] = True
            df_dose_all.index = index_letters
        else:
            tools.DatafileError('response assaytype not identified')
    else:
        raise tools.DatafileError('response machinetype not identified')

    return df_dose_all, df_resp_all

def examine_input_datafile(fn, dff):
    """ Checks the file integrity of the input datafiles.

    If the datafile is .xlsx, it assumes that the input is from the fluostar 12-well platereader.
    For 96-well samples, where the input is a text file, it first looks for the keyword
    identifying as a VERSAmax file from the langosch lab (ident_txtstr), and then looks for the keyword
    indicating the 96-well format(either 8,12,24 dose concentrations, from list_96_well_doseconc_types).

    Parameters
    ----------
    fn : int
        File number in list of files to analyse.
    dff : pandas DataFrame
        Dataframe for Files. Contains all the paths for input and output files.
        Created from the "files" tab of the settings excel file.

    Returns
    -------
    dff : pandas DataFrame
        Input dataframe with additions, specifying whether the data filepaths and filetypes appear to be correct.
    """
    # define the response datafileformat, machinetype, and assaytype from string (e.g. txt|versamax|12_doseconc)
    response_dataformat = dff.loc[fn, "response dataformat"]
    eccpy, resp_datafileformat, resp_machinetype, resp_assaytype = response_dataformat.split("|")
    # define the datafile path
    data_file_path = dff.loc[fn,"data_file_path"]

    # start by assuming the response datafile is not okay, mark as False
    dff.loc[fn, "resp_datafile_ok"] = False
    resp_datafile_contains_ident_txtstr, resp_datafile_ok = False, False

    if os.path.exists(data_file_path):
        if resp_datafileformat == "txt":
            if data_file_path[-4:] != ".txt":
                raise ValueError("Data file is not a text file. {} does not end in .txt. Please check settings file.".format(data_file_path))
            if resp_machinetype == "versamax" and resp_assaytype in ["8dose12sample", "12dose8sample", "24dose4sample"]:
                # look for the full response_dataformat string in the text file
                with open(data_file_path,'r') as f:
                    for line in f:
                        # check if the response_dataformat text is in the line
                        if resp_assaytype in line:
                            # mark the response datafile as okay
                            resp_datafile_ok = True
                            dff.loc[fn, "resp_datafile_ok"] = resp_datafile_ok
                            # string is found, break the loop
                            break

            elif resp_machinetype == "versamax" and resp_assaytype in ["8_ampconc", "12_ampconc", "24_ampconc"]:
                ident_txtstr = "Langosch_Lab_VERSAmax_export_file"
                # go through all lines in the text file, looking for the relevant textstring identifiers
                with open(data_file_path,'r') as f:
                    for line in f:
                        # check if the identifying text is in the line
                        if ident_txtstr in line:
                            resp_datafile_contains_ident_txtstr = True
                        # check if the doseconc textstring is in the line.
                        if resp_assaytype in line:
                            # mark the response datafile as okay
                            resp_datafile_ok = True
                            dff.loc[fn, "resp_datafile_ok"] = resp_datafile_ok
                            # string is found, break the loop
                            break

                if resp_datafile_contains_ident_txtstr == False:
                    raise tools.DatafileError("The input datafile does not contain the expected identifier"
                                        " phrase.\n'%s'\n Suggest checking the filepath and file integrity. "
                                        "For VersaMAX datafiles, double-check that when the data was"
                                        " exported, that both check-boxes were selected, so that"
                                        " 'all sections on plate' were exported.\n"
                                        "Relevant file:\n%s" % (ident_txtstr, data_file_path))
            else:
                raise tools.DatafileError("{} is not recognised as a valid format. "
                                    "Check settings excel file".format(response_dataformat))


            if resp_datafile_ok == False:
                raise tools.DatafileError("Response dataformat was listed as {a}, but the datafile "
                                    "does not seem to contain the resp_assaytype textstring ({b}). "
                                    "Suggest checking the settings "
                                    "excel file and raw text files.".format(a=response_dataformat,
                                                                            b=resp_assaytype))

        elif resp_datafileformat == "xlsx":
            if resp_machinetype == "generic":
                if resp_assaytype in ["vertical","horizontal"]:
                    try:
                        df_identifier = pd.read_excel(data_file_path, sheetname="template_identifier")
                        if dff.loc[fn, "response dataformat"] in df_identifier.columns:
                            dff.loc[fn, "resp_datafile_ok"] = True
                        else:
                            first_column_in_df = df_identifier.columns[0]
                            raise tools.DatafileError("response dataformat was listed as {a}, but the datafile contains the "
                                                "following string :{b}\nSuggest checking the settings excel file and "
                                                "response datafile templates.".format(a=response_dataformat,
                                                                                      b=first_column_in_df))
                    except:
                        raise tools.DatafileError("cannot open {a} file, {b}".format(a=resp_datafileformat,b=data_file_path))
            elif resp_machinetype == "fluostar":
                try:
                    df_resp = pd.read_excel(data_file_path, sheetname='A600', skiprows = 2)
                    if "A" in df_resp.iloc[0,0] and "X" in df_resp.iloc[0,1]:
                        dff.loc[fn, "resp_datafile_ok"] = True
                    else:
                        raise LookupError("input file with data does not seem to match the xlsx template for the "
                                          "12-well platereader :\n%s" %df_resp)
                except:
                    raise LookupError("input file with data does not seem to match the xlsx template for the "
                                     "12-well platereader.")
            else:
                raise tools.DatafileError("machine type {} not recognised".format(resp_machinetype))
        else:
            raise tools.DatafileError("datafile type {} not recognised".format(resp_datafileformat))
    else:
        raise IOError('File {0} does not exist'.format(data_file_path))

    return dff

def read_versamax_txt_datafile(fn, dff, resp_assaytype):
    """ Parses the VersaMax 96-well input file containing response data.

    Parameters
    ----------
    fn : int
        File number in list of files to analyse.
    dff : pandas DataFrame
        Dataframe for Files. Contains all the paths for input and output files.
        Created from the "files" tab of the settings excel file.
    resp_assaytype : string
        The response assay type. E.g. for "12dose8sample" for microplate data, or "vertical" for generic excel data.

    Returns
    -------
    df_resp_orig : pandas DataFrame
        Dataframe containing the the original response data.
    df_resp_all : String
        String, to be replaced with a dataframe later.
    df_dose_orig : pandas DataFrame
        Dataframe containing the dose concentrations. Needs to be standardised so that the index and columns
        match the response dataframe.

    Note
    ----

    SoftMax® Pro is designed for use with Molecular Devices’ Microplate Readers (Emax™, Vmax®, UVmax™, ThermoMax™,
    and VersaMax™), Micro-plate Spectrophotometers (SpectraMax®), and Microplate Spectroﬂuorometers (SpectraMax Gemini
    ﬂuorescence reader family and FlexStation instruments).

    A specific .pda template needs to be used for the VersaMax file, containing the relevant text string that
    indicates the assay type.
    """
    dose_conc_excel_path = os.path.join(dff.loc[fn, "input file directory"], dff.loc[fn, "dose conc file"])

    # identify where the data starts and ends in the csv file.
    identifier_start_of_data = "Data_for_export"
    #identifier_end_of_data = "Instrument type: Instrument serial number:"
    with open(dff.loc[fn,"data_file_path"],'r') as f:
        linenumber = 0
        # Create the initial local variable to keep the code inspection algorithm happy.
        linenumber_header = 0
        n_first_line_with_data = None
        for line in f:
            # the data starts one line after the "data_for_export" text header
            if identifier_start_of_data in line:
                linenumber_header = linenumber + 1
                n_first_line_with_data = linenumber + 2
            if linenumber == n_first_line_with_data:
                first_line_with_data = line
            # the data ends 3 lines before the "Instrument type:" text
            # not currently used in script, which assumes that there are always
            # 2 lines with text after the end of the data
            #if identifier_end_of_data in line:
            #    linenumber_end_data = linenumber - 3
            linenumber += 1
        # check if the data is saved using a German or English decimal point
        if "," in first_line_with_data:
            dec = ","
        elif "." in first_line_with_data:
            dec = "."
        else:
            raise TypeError("VersaMax exported text file seems to have an unknown decimal format. Try re-exporting data.")
        # open the VERSAmax microplate reader exported csv file as a new dataframe, df_resp_orig (MICROPLATE READER DATA)
        df_resp_orig = pd.read_csv(dff.loc[fn,"data_file_path"], sep='\t', skiprows=linenumber_header, decimal=dec)
        #remove the last two rows, which should contain "~End" and "Instrument type: Instrument serial number:"
        df_resp_orig = df_resp_orig.dropna()
    # create new dataframe "df_dose_orig" to contain the dose concentrations
    df_dose_orig = pd.read_excel(dose_conc_excel_path, sheetname = resp_assaytype)
    df_resp_all = "temporary object, will be overwritten by standardisation of df_resp_orig, 12-well or 96well"
    return df_resp_orig, df_resp_all, df_dose_orig

