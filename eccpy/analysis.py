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
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv
import os
import sys
import ast
import eccpy.settings as eccpysettings
import eccpy.tools as tools

# show wide pandas dataframes when using print function
pd.set_option('display.expand_frame_repr', False)

def run_analysis(settings_excel_file, **kwargs):
    """ Analysis of EC50 data from multiple experiments.

    Processes the datafiles in settings excel file marked as "TRUE" for "run analysis."
    Collects and analyses the output from the "run_curvefit" program.

    Parameters
    ----------
    settings_excel_file : filepath
        Settings file containing the list of datafiles for analysis, and also chosen parameters

    Keyword Arguments
    -----------------
    savefig_with_longnames : boolean
        If True, output figures will be created with both long and short sample names.

    Saved Files
    -----------
    EC50_barchart : Barchart of EC50 values from all experiments.
        Only data is included that is judged as good quality (i.e. not labelled as "data_needs_checking").
        Four barcharts are created.
        1) Original data, long sample names
        2) Original data, short sample names
        3) Adjusted data (e.g. fixed upper limit), long sample names
        4) Adjusted data (e.g. fixed upper limit), short sample names
    EC50_datapoints: Scattergram, sample names on x-axis, EC50 on y-axis.
        Effectively the same data as the EC50_barchart, except that the datapoints from each
        experiment are displayed individually, as in a scattergram. A legend with colours for each experiment
        is included. Very useful to determine if the EC50 values from one particular day were uniformly higher or
        lower than the EC50 values calculated on the other days.
        As with the barchart, four variations are created.
        1) Original data, long sample names
        2) Original data, short sample names
        3) Adjusted data (e.g. fixed upper limit), long sample names
        4) Adjusted data (e.g. fixed upper limit), short sample names

    Returns
    -------
    df_allp: pandas DataFrame
        For developers only. Dataframe for all individual EC50 datapoints.
    dff : pandas DataFrame
        For developers only. Dataframe for Files. Contains all the paths for input and output files.
        Created from the "files" tab of the settings excel file.

    Notes
    -------
    The output of the analysis is saved in the following folder:
    ORIGINAL_SUBFOLDER_WITH_SETTINGS_EXCEL_FILE/analysed/todays_date/

    Running this script will overwrite any previous files with the same name (i.e., analysed on the same day)

    The "run_curvefit" program needs to be run to create the relevent output files, before the analysis program
    can be started. Shifting the location of the output files will result in an error.

    The output files begin as follows:
    todays_date_analysed_
    todays_date is represented as YEAR|MONTH|DAY, e.g. 20151215.

    Figures are usually better displayed if a dictionary of long to short sample names is created. This can be saved
    in the settings_excel_file.

    Currently, the analysis automatically runs for original and adjusted datasets (e.g. fixed upper limit dataset).
    However summary graphs are created separately for each dataset.
    """
    print("\nStarting run_analysis program\n")
    # create output folder for analysed data, define basename
    outpath, basename = eccpysettings.setup_output_folder(settings_excel_file, "analysed")
    analysed_data_basename = os.path.join(outpath, basename)
    # add the relevant paths to the data files to the dataframe for files (dff)
    settings, dff, df_samplenames = eccpysettings.read_settings_file(settings_excel_file)
    # set the long sample name as the index
    df_samplenames.set_index("long name", inplace = True)
    # create t20 colour list
    t20 = tools.setup_t20_colour_list()
    # extract list of adjusted datasets for analysis
    datasets = ast.literal_eval(settings["adjust.datasets"])
    """
    COLLECT THE EC50 VALUES FROM ALL THE OUTPUT FILES
    """
    # fix the suffixes denoting the datasets (_orig for original data, _ful for fixed upper limit)
    # if any of the files are labelled True for "run analysis"
    if True in list(dff.loc[:, "run analysis"]):
        # create an empty dataframe to hold the average EC50 values
        dfc = pd.DataFrame()
        # create another dataframe to hold all the boolean data
        dfcb = pd.DataFrame()
        # create dataframe for all individual EC50 datapoints (df_allp), including samples repeated in a single experiment
        df_allp = pd.DataFrame()
        print("Analysing data from multiple experiments. List of experiments analysed :")
        # iterate through only the files labelled "True" for "run analysis", and join all output dataframes together
        for fn in dff.loc[dff["run analysis"] == True].index:
            # define the response data file
            data_file = dff.loc[fn, "response data file"]
            print(data_file)
            # if it is a real file, open
            if os.path.isfile(dff.loc[fn,"ofd_EC50_eval_excel"]):
                # open  as a new pandas dataframe
                df_eval_values = pd.read_excel(dff.loc[fn,"ofd_EC50_eval_excel"], sheetname="v_" + data_file[:20])
                # add the sample name to the dataframe, so it can be identified later
                df_eval_values["file"] = dff.loc[fn,"ofd_EC50_eval_excel"]
                # convert the sample_name column to a string datatype (forcing np.nan to be "nan")
                df_eval_values["sample_name"] = df_eval_values["sample_name"].astype(str)
                # drop any rows that contain "nan" as the sample name
                df_eval_values = df_eval_values.loc[df_eval_values["sample_name"] != "nan"]
                # join the dataframe with all previously joined dataframes
                dfc = pd.concat([dfc,df_eval_values], axis=0)
                # open the tab of the summary excel file that contains all the boolean values
                df_eval_bool = pd.read_excel(dff.loc[fn,"ofd_EC50_eval_excel"], sheetname="b_" + data_file[:20])
                # join the dataframe with all previously joined dataframes
                dfcb = pd.concat([dfcb,df_eval_bool], axis=0)
                # set the sample_name as the index
                df_eval_values = df_eval_values.set_index("sample_name")
                # iterate through datasets (save in the same df, by adding suffix to the column name)
                for d in datasets:
                    # change the dataset name (e.g. "_orig" to "") to an empty string if there is only one dataset for analysis
                    d_name = "" if len(datasets) == 1 else d
                    # define the column name in the dataframe
                    col_EC50 = "EC50" + d
                    # select data which "seems okay" according the the automatic data analysis
                    df_eval_values_OK = df_eval_values.loc[df_eval_values["data_seems_okay" + d] == True]
                    # create a list of unique sample names
                    unique_names = list(df_eval_values_OK.index.unique())
                    # create a new dataframe, called df_eval_uniq, which has a single row for each unique sample
                    df_eval_uniq = pd.DataFrame().astype(object)
                    for sn in unique_names:
                        # select only the data for that sample
                        df_sel = df_eval_values_OK.loc[sn,:]
                        # if there is only one sample, the selected data will form a series
                        if isinstance(df_sel, pd.Series):
                            # add the n, the EC50, and the std
                            df_eval_uniq.loc[sn,data_file + d] = df_sel["EC50{}".format(d)]
                        # if the name is not unique, the selected data will form a dataframe.
                        elif isinstance(df_sel, pd.DataFrame):
                            # transfer the EC50 values as a stringlist
                            df_eval_uniq.loc[sn,data_file + d] = str(["%0.2f"%l for l in df_sel[col_EC50]])
                        else:
                            raise TypeError("expected a series or dataframe.")
                    # add the dataframe containing _orig and other dataset columns for that exp to the final df with all data
                    df_allp = pd.concat([df_allp,df_eval_uniq], axis=1)
            else:
                print("File not found! {}".format(dff.loc[fn,"ofd_EC50_eval_excel"]))

        print("\nPercentage data okay:")
        if dfc.empty:
            raise ValueError("No data collected for analysis! Double-check that run_curvefit program has been"
                             " carried out for all the relevant files")

        for d in datasets:
            vc = dfc["data_seems_okay{}".format(d)].value_counts()
            if True in vc:
                n_data_okay = vc[True]
            else:
                n_data_okay = 0
            if False in vc:
                n_data_not_okay = vc[False]
            else:
                n_data_not_okay = 0
            perc_data_okay = n_data_okay / (n_data_okay + n_data_not_okay)*100
            print("{b:0.0f}% ({a} dataset)".format(a=d[1:], b=perc_data_okay))
        # select only the data labeled as "data_seems_okay"
        # save the current index as the sample letter
        dfc["sLet"] = dfc.index
        # convert the index to the sample name
        dfc.index = dfc.sample_name
        # create a new dataframe, called dfm, which has a single row for each unique sample, and contains MEAN values
        dfm = pd.DataFrame()
        for d in datasets:
             # select only the data labeled as "data_seems_okay"
            series_data_seems_okay = dfc["data_seems_okay{}".format(d)] == True
            dfc_ok = dfc.loc[series_data_seems_okay]
            # create list of unique sample names, where data is available
            list_unique_sample_names = list(dfc_ok.sample_name.dropna().unique())
            for sn in list_unique_sample_names:
                # select data for that one sample, resulting in either a series, or dataframe (indicating mult. exper.)
                data_1_sample_name = dfc_ok.loc[sn,:]
                # if there is only one sample, df_sel will be a series.
                if isinstance(data_1_sample_name, pd.Series):
                    # add the n, the EC50, and the std
                    dfm.loc[sn,"n{}".format(d)] = 1
                    dfm.loc[sn,"mean{}".format(d)] = data_1_sample_name["EC50{}".format(d)]
                    dfm.loc[sn,"std{}".format(d)] = 0
                    dfm.loc[sn,"SEM{}".format(d)] = 0
                # if there are multiple datapoints with the same name, df_sel will be a DataFrame
                elif isinstance(data_1_sample_name, pd.DataFrame):
                    # add the n, the mean EC50 of the samples, and the std
                    dfm.loc[sn,"n{}".format(d)] = data_1_sample_name["EC50{}".format(d)].shape[0]
                    dfm.loc[sn,"mean{}".format(d)] = data_1_sample_name["EC50{}".format(d)].mean()
                    dfm.loc[sn,"std{}".format(d)] = data_1_sample_name["EC50{}".format(d)].std()
                    dfm.loc[sn,"SEM{}".format(d)] = data_1_sample_name["EC50{}".format(d)].sem()

        # convert longnames in index to a new column
        dfm["longname"] = dfm.index
        df_allp["longname"] = df_allp.index
        # create a dictionary from the long and short sample names
        samplenames_dict = dict(zip(df_samplenames.index, df_samplenames["short name"]))
        # create a new column with the relevant shortnames for each sample, based on the samplenames_dict
        dfm["shortname"] = dfm["longname"].apply(lambda x : samplenames_dict[x] if x in list(samplenames_dict.keys()) else x)
        df_allp["shortname"] = df_allp["longname"].apply(lambda x : samplenames_dict[x] if x in list(samplenames_dict.keys()) else x)

        ##########################################################################################
        #          reorder the samples according to desired order in the settings file           #
        ##########################################################################################
        # sort original index, in case
        dfm.sort_index(inplace=True)
        df_allp.sort_index(inplace=True)
        # create a dictionary of the sample order
        sampleorder_dict = dict(zip(df_samplenames.index, df_samplenames["order in figure"]))
        # create a new column with the preferred sample order, if available in settings file.
        dfm["sampleorder"] = dfm.longname.apply(lambda x : sampleorder_dict[x] if x in list(sampleorder_dict.keys()) else np.nan)
        df_allp["sampleorder"] = df_allp.longname.apply(lambda x : sampleorder_dict[x] if x in list(sampleorder_dict.keys()) else np.nan)
        # sort by sample order
        dfm.sort_values(by="sampleorder", inplace=True)
        df_allp.sort_values(by="sampleorder", inplace=True)

        # save the dataframe with all mean data from all experiments to a csv
        # dfm.to_csv(analysed_data_basename + "_EC50_mean" + ".csv", sep=",", quoting=csv.QUOTE_NONNUMERIC)
        df_allp.to_csv(analysed_data_basename + "_EC50_indiv_exp" + ".csv", sep=",", quoting=csv.QUOTE_NONNUMERIC)
        # save both dataframes (mean data and indiv datapoints) from all experiments to excel
        writer = pd.ExcelWriter(analysed_data_basename + ".xlsx")#engine='xlsx_nonnorm'
        # dfm.to_excel(writer, sheet_name = "EC50_mean")
        df_allp.to_excel(writer, sheet_name="EC50_indiv_exp_nonnorm")

        # sort the columns
        df_allp.sort_index(axis=1, inplace=True)
        # replace index of df_allp with a simple range of integers
        df_allp.index = range(df_allp.shape[0])
        # give the index a name, indicating sample numbers
        df_allp.index.name = "sSnum"

        # set the fontsize
        if df_allp.shape[0] < 10:
            fontsize = 10
        else:
            fontsize = 6
        plt.rcParams['font.size'] = fontsize
        # set matplotlib legend parameters
        plt.rcParams['legend.numpoints'] = 3

        # iterate through the datasets (e.g. _orig, _ful)
        for d in datasets:
            col_mean = "mean" + d
            col_std = "std" + d
            col_SEM = "SEM" + d
            # create a subset to plot that contain data
            df_for_barchart = dfm.dropna(subset = [col_mean]).copy()
            # # use the dictionary of long to short names to obtain the short name for the datapoints
            # df_for_barchart["longname"] = df_for_barchart.index
            # df_for_barchart["shortname"] = df_for_barchart.longname.apply(lambda x : samplenames_dict[x] if x in list(samplenames_dict.keys()) else x)

            # setup normalisation by checking if a standard is listed in the settings file
            conduct_normalisation, list_norm_datasets = setup_normalisation(df_samplenames, df_allp)

            # identify the df_allp columns associated with that dataset
            col_contains_d = list(pd.Series(df_allp.columns).apply(lambda x: d in x))
            # select only data for that dataset (e.g. only orig data)
            sel_df_allp = df_allp.loc[:,col_contains_d]
            # create separate figures for the long names, and the short names
            for norm_dataset in list_norm_datasets:
                sys.stdout.write(".")
                sys.stdout.flush()
                # BARCHART WITH MEAN REPLACED BY BARCHART CREATED WITH SCATTER DATA, ALLOWING NORMALISATION
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                #                                                                                                     #
                #                                       Analysis Barcharts                                            #
                #               by dataset (ful/orig), by name (long,short), by normalisation (orig/norm)             #
                #                                                                                                     #
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                #            _________
                #           |XXXXXXXXX|          Full Canvas
                #           |XXXXXXXXX|
                #           |XXXXXXXXX|
                #           |XXXXXXXXX|

                # create a new figure object for the barchart
                #bar_fig, bar_ax = plt.subplots()

                # if norm_dataset == "_nonnorm":
                #     # define EC50 data (y-axis, height of bars)
                #     bar_y_meanEC50 = df_for_barchart[col_mean]
                #     # define y-error bars
                #     bar_yerr = df_for_barchart[col_SEM]
                # elif norm_dataset == "_norm":
                #     bar_samplenames_selected = df_samplenames.loc[df_for_barchart.longname, :]
                #     bar_y_meanEC50_nonnorm = df_for_barchart[col_mean]
                #     bar_yerr_nonnorm = df_for_barchart[col_SEM]
                #     # find the sample longname that is marked "True" as a standard
                #     standard_name = bar_samplenames_selected[bar_samplenames_selected["standard for normalisation?"] == True].index[0]
                #     standard_EC50 = df_for_barchart.loc[standard_name,col_mean]
                #     bar_y_meanEC50 = bar_y_meanEC50_nonnorm / standard_EC50 * 100
                #     bar_yerr = bar_yerr / standard_EC50 * 100

                # # define the indices of the boxes on the x-axis
                # bar_x_n_boxes = df_for_barchart.shape[0]
                # box_indices = range(bar_x_n_boxes)
                # # define error bar parameters
                # bar_error_kw = dict(ecolor='k', lw=1, capsize=2, capthick=1)
                # # plot the data as a barchart
                # barcontainer = bar_ax.bar(box_indices, bar_y_meanEC50, yerr=bar_yerr, align="center",
                #                           error_kw=bar_error_kw, color = '#1b9e77')
                # # set the xticks
                # bar_ax.set_xticks(box_indices)
                #
                # # set the limits of the x-axis
                # bar_ax.set_xlim([-1, bar_x_n_boxes])
                # # set the limit of the y-axis
                # bar_ax.set_ylim(0)
                # # set the y-axis title
                # # bar_ax.set_ylabel("EC50 (ug/ml)")
                # # set the ylabel extension string
                # if norm_dataset == "_nonnorm":
                #     dose_units  = settings["x-axis (dose) units"]
                # elif norm_dataset == "_norm":
                #     samplenames_selected = df_samplenames.loc[df_allp.longname, :]
                #     # find the sample longname that is marked "True" as a standard
                #     standard_name = samplenames_selected[samplenames_selected["standard for normalisation?"] == True].index[0]
                #     # check if the standard for normalisation has a short name
                #     if standard_name in samplenames_selected.index:
                #         standard_name_short = samplenames_selected.loc[standard_name, "short name"]
                #     else:
                #         # use the long name
                #         standard_name_short = standard_name
                #     # create a string to extend the y-axis label
                #     dose_units = "% {}".format(standard_name_short)
                # else:
                #     raise TypeError("dataset for standardisation{} is not recognised".format(norm_dataset))

                # ylabel_str = "{a}{b}, {c} ({d})".format(a=settings["calculation_type"],
                #                       b=str(settings["percentage_response"]),
                #                       c=settings["x-axis (dose) label"],
                #                       d=dose_units)
                # bar_ax.set_ylabel(ylabel_str)

                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                #                                                                                                     #
                #                  Analysis Scattergram, individual datapoints for multiple experiments               #
                #               by dataset (ful/orig), by name (long,short), by normalisation (orig/norm)             #
                #                                                                                                     #
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                #            _________
                #           |XXXXXXXXX|          Full Canvas
                #           |XXXXXXXXX|
                #           |XXXXXXXXX|
                #           |XXXXXXXXX|

                # create a scattergram figure object
                scat_fig, scat_ax = plt.subplots()

                # create an lists to hold the individual sample numbers (x-values) from all experiments
                xvalues_all_exp = []
                # create an lists to hold the individual datapoints  (y-values) from all experiments
                yvalues_all_exp = []

                # create a string for the dose units (e.g. "mg/L", or "% positive control")
                if norm_dataset == "_nonnorm":
                    dose_units  = settings["x-axis (dose) units"]
                elif norm_dataset == "_norm":
                    samplenames_selected = df_samplenames.loc[df_allp.longname, :]
                    # find the sample longname that is marked "True" as a standard
                    standard_name = samplenames_selected[samplenames_selected["standard for normalisation?"] == True].index[0]
                    # check if the standard for normalisation has a short name
                    if standard_name in samplenames_selected.index:
                        standard_name_short = samplenames_selected.loc[standard_name, "short name"]
                    else:
                        # use the long name
                        standard_name_short = standard_name
                    # create a string to extend the y-axis label
                    dose_units = "% {}".format(standard_name_short)
                else:
                    raise TypeError("dataset for standardisation{} is not recognised".format(norm_dataset))
                ylabel_str = "{a}{b}, {c} ({d})".format(a=settings["calculation_type"],
                                      b=str(settings["percentage_response"]),
                                      c=settings["x-axis (dose) label"],
                                      d=dose_units)

                # iterate through the experiment number (exp_nr) and columns (c) in the dataframe
                # Each column represents a single experiment
                for exp_nr, c in enumerate(sel_df_allp.columns):
                    data_colour = t20[exp_nr]
                    sys.stdout.write(".")
                    sys.stdout.flush()
                    # convert any stringlists to lists
                    df_allp.loc[:,c] = df_allp.loc[:,c].apply(lambda x: ast.literal_eval(x) if isinstance(x,str) else x)
                    # convert any lists of strings to lists of floats
                    df_allp.loc[:,c] = df_allp.loc[:,c].apply(lambda x: [float(s) for s in x] if isinstance(x,list) else x)
                    # convert any lists of floats to lists of numpy arrays
                    df_allp.loc[:,c] = df_allp.loc[:,c].apply(lambda x: np.array(x) if isinstance(x,list) else x)
                    
                    ####################################################################################################
                    #        NORMALISATION TO A STANDARD FOR EACH EXPERIMENT                                           #
                    #         a) find the standard sample number                                                       #
                    #         b) find the standard                                                                     #
                    #         c) divide all datapoints by the standard (or mean of standard, if duplicated on that day)#
                    #         d) multiply by 100 to give a percentage                                                  #
                    ####################################################################################################

                    if norm_dataset == "_nonnorm":
                        series_EC50_data = df_allp.loc[:,c].dropna()
                    elif norm_dataset == "_norm":
                        # find the sample longname that is marked "True" as a standard
                        standard_name = samplenames_selected[samplenames_selected["standard for normalisation?"] == True].index[0]
                        # find the sSNum (index number) of the standard
                        standard_sSnum = df_allp[df_allp.longname == standard_name].index[0]
                        # calculated the average EC50 for the standard
                        standard_EC50 = np.mean(df_allp.loc[standard_sSnum, c])
                        # calculate the normalised EC50 as a percentage (EC50/standard*100)
                        series_EC50_data = df_allp.loc[:,c]/standard_EC50*100
                        series_EC50_data.dropna(inplace=True)
                    ########################################################################################
                    #                                                                                      #
                    #    Add values to a "scatterplot" which has the sample number on the x-axis,          #
                    #    and the EC50 on the y-axis.                                                       #
                    #    For duplicate samples in that experiment, it's a bit tricky!!                     #
                    #      a)  find data with duplicates [S1_EC50_1,S1_EC50_2]                             #
                    #      b)  iterate through each sample with duplicates                                 #
                    #      c)  create matching list of sample numbers [S1,S1],[S1_EC50_1,S1_EC50_2]        #
                    #      d)  add to list of x-values and y-values from the samples with duplicates       #
                    #            x=[S1,S1,S2,S2,S5,S5],                                                    #
                    #            y=[S1_EC50_1,S1_EC50_2,S2_EC50_1,S2_EC50_2,S5_EC50_1,S5_EC50_2]           #
                    #      e)  append(extend) to the list of samples                                       #
                    #           x=[S2,S3,S4,S1,S1,S2,S2,S5,S5],                                            #
                    #            y=[S2_EC50_1,S3_EC50_1,S4_EC50_1,S1_EC50_1,S1_EC50_2,                     #
                    #               S2_EC50_1,S2_EC50_2,S5_EC50_1,S5_EC50_2]                               #
                    #      f) plot as a scattergram, so multiple datapoints are on each sample             #
                    #      g) append (extend) the sample names and datapoints to a list                    #
                    #         containing data from all experiments combined                                #
                    #      h) move and repeat for the next experiment                                      #
                    #      i) after iteration through all experiments, use collected list to               #
                    #         calculate the mean values (whether raw or normalised) for barchart           #
                    #                                                                                      #
                    ########################################################################################
                    
                    # set the transparency
                    alpha = 0.5
                    # create a series of bools, describing if the data is a float
                    is_float_series = series_EC50_data.apply(lambda x: isinstance(x,float))
                    # use boolean series to select the float data
                    values_single_sample = series_EC50_data.loc[is_float_series]
                    # create list of x-values (dose) and y-values (response) from the data with single samples
                    xvalues_from_single = list(values_single_sample.index)
                    yvalues_from_single = list(values_single_sample)
                    # create a series of bools, describing if the data is a list
                    is_array_series = series_EC50_data.apply(lambda x: isinstance(x,np.ndarray))
                    # if some of the datapoints are lists (multiple sample replicates in a single experiment)
                    if True in list(is_array_series):
                        # use boolean series to select the data with lists
                        values_mult_sample = series_EC50_data.loc[is_array_series]
                        # the multiple EC50 values belong to a single sample number. Identify list of sample numbers.
                        index_ssNum_with_lists = values_mult_sample.index.tolist()
                        # create a list of xvalues for the scattergram, for multiple samples per day
                        list_xvalues_from_mult = []
                        # create a list of yvalues for the scattergram, for multiple samples per day
                        list_yvalues_from_mult = []
                        # iterate through each selected sample number, associated with multiple EC50 values
                        for sSnum in index_ssNum_with_lists:
                            list_EC50_values = list(values_mult_sample[sSnum])
                            # create corresponding x-axis sample numbers (array of sSnum of length len(values_mult_sample))
                            list_index_values = list(np.ones(len(list_EC50_values)).astype(np.int64)*sSnum)
                            # add x and y values to the lists for all samples
                            list_xvalues_from_mult.extend(list_index_values)
                            list_yvalues_from_mult.extend(list_EC50_values)
                        # add the x and y-values from the multiple samples per experiment to the single sample data
                        xvalues_all_single_exp = xvalues_from_single + list_xvalues_from_mult
                        yvalues_all_single_exp = yvalues_from_single + list_yvalues_from_mult
                    else:
                        # there is no data from multiple samples in an experiment, therefore all data comes from single
                        xvalues_all_single_exp = xvalues_from_single
                        yvalues_all_single_exp = yvalues_from_single

                    # plot the float data as a scattergram (x-axis is the range, resembling a bar or line chart)
                    scat_ax.scatter(xvalues_all_single_exp, yvalues_all_single_exp, color=data_colour, s=40,
                                    alpha=alpha, label=c[:-5])

                    # add the values from that experiment to a list
                    xvalues_all_exp.extend(xvalues_all_single_exp)
                    yvalues_all_exp.extend(yvalues_all_single_exp)

                #######################################################
                #    format_and_save_analysis_scatterplot             #
                #######################################################

                settings_name = os.path.split(settings_excel_file)[1]

                # set the xticks and labels to match the index of df_allp
                scat_ax.set_xticks(np.arange(df_allp.shape[0]))

                # set the grid to go in between the sample names, as minor xticks
                scat_ax.set_xticks(np.arange(df_allp.shape[0])+0.5, minor=True)
                scat_ax.grid(which='minor', alpha=0.9)
                # set the x axis limits
                scat_ax.set_xlim(-0.5, df_allp.shape[0]-0.5)

                # set the y-axis title
                scat_ax.set_ylabel(ylabel_str)

                if "savefig_with_longnames" in kwargs.keys():
                    if kwargs["savefig_with_longnames"] == True:
                        list_nametypes = ["longname", "shortname"]
                    else:
                        list_nametypes = ["shortname"]
                else:
                    list_nametypes = ["shortname"]

                for nametype in list_nametypes:
                    # # define name on the x-axis
                    # bar_x_names = df_for_barchart[nametype]
                    # # set the labels of the x-bar_axis
                    # bar_ax.set_xticklabels(bar_x_names, rotation=90)
                    #
                    # # ax.annotate(s="%s%s" % (nametype,d), xy=(0.015,0.93), fontsize=af, xycoords=xyc)
                    # bar_ax.set_title("analysed data ({e} experiments),  {a}{b},  {c}".format(a=nametype,b=d,
                    #                                                   c=os.path.split(settings_excel_file)[1],
                    #                                                   e=dff.loc[dff["run analysis"] == True].shape[0]))
                    # # automatically tighten the layout and save figure
                    # bar_fig.tight_layout()
                    # # save the figure
                    # bar_fig.savefig(analysed_data_basename + "_bar_" +  nametype + d  + norm_dataset + '.png',
                    #                 format='png', dpi=150)

                    # add legend
                    if nametype == "longname":
                        scat_lgd = scat_ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1,
                                                  scatterpoints=1,numpoints=1, borderaxespad=1,
                                                  fontsize=fontsize)# mode="expand",
                        # scat_lgd = scat_ax.legend(loc='lower center', bbox_to_anchor=(0., 1.02, 1., .102), ncol=2,
                        #                           scatterpoints=1,numpoints=1, borderaxespad=2,
                        #                           fontsize=fontsize-1)# mode="expand",
                    elif nametype == "shortname":
                        scat_lgd = scat_ax.legend(loc='upper center', bbox_to_anchor=(0.5,-0.15), ncol=2,
                                                  scatterpoints=1,numpoints=1, fontsize=fontsize)
                    # define the x-axis labels (long or short)
                    scat_xticklabels = list(df_allp[nametype])
                    scat_ax.set_xticklabels(scat_xticklabels, rotation=90)
                    scat_ax.set_title("analysed data ({e} experiments), {b}".format(b=d_name,e=sel_df_allp.shape[1]))
                    # scat_ax.set_title("analysed data ({e} experiments),  "
                    #                   "{a}{b}{f},  {c}".format(a=nametype, b=d_name, f=norm_dataset,
                    #                                            c=settings_name, e=sel_df_allp.shape[1]))
                    # automatically tighten the layout and save figure
                    # scat_fig.tight_layout()
                    # scat_fig.savefig(analysed_data_basename + "_datapoints_" + nametype + d + norm_dataset + '.png',
                    #                  format='png', dpi=300, bbox_extra_artists=(scat_lgd,), bbox_inches='tight')
                    scat_fig.savefig(analysed_data_basename + "_datapoints_" + nametype + d_name + norm_dataset + '.png',
                                     format='png', dpi=300,bbox_extra_artists=(scat_lgd,), bbox_inches='tight')

                # create figure object for the barchart with normalised data (collected from the scatter data)
                barnorm_fig, barnorm_ax = plt.subplots()
                series_all_exp_redundant = pd.Series(yvalues_all_exp, index = xvalues_all_exp)
                # convert index to integer format
                series_all_exp_redundant.index = series_all_exp_redundant.index.astype(np.int64)
                # from the redundant index, create a unique, sorted index of sample numbers
                series_all_exp_redundant_index_uniq = series_all_exp_redundant.index.unique()
                series_all_exp_redundant_index_uniq.sort()
                df_all_exp_nonredundant = pd.DataFrame(index=series_all_exp_redundant_index_uniq)
                for sSnum_b in series_all_exp_redundant_index_uniq:
                    data_for_that_sSnum = series_all_exp_redundant.loc[sSnum_b]
                    if isinstance(data_for_that_sSnum, pd.Series):
                        n_samples = data_for_that_sSnum.shape[0]
                    else:
                        n_samples = 1
                    df_all_exp_nonredundant.loc[sSnum_b, "n{}{}".format(d, norm_dataset)] = n_samples
                    df_all_exp_nonredundant.loc[sSnum_b, "mean{}{}".format(d,norm_dataset)] = data_for_that_sSnum.mean()
                    df_all_exp_nonredundant.loc[sSnum_b, "std{}{}".format(d,norm_dataset)] = data_for_that_sSnum.std()
                    df_all_exp_nonredundant.loc[sSnum_b, "SEM{}{}".format(d,norm_dataset)] = pd.Series(data_for_that_sSnum).sem()

                # extract the original longnames from dfm
                longname_series = dfm.longname
                longname_series.index = range(len(longname_series))
                # add the longnames to df_all_exp_nonredundant (index should be sample number)
                df_all_exp_nonredundant["longname"] = longname_series
                # extract the shortnames from the dictionary in the settings file, if available
                df_all_exp_nonredundant["shortname"] = df_all_exp_nonredundant["longname"].apply(lambda x : samplenames_dict[x] if x in list(samplenames_dict.keys()) else x)

                # count the number of samples (number of boxes in barchart)
                bar_x_n_boxes = df_all_exp_nonredundant.shape[0]
                # set the range of number of boxes as the indices
                box_indices = range(bar_x_n_boxes)
                # define error bar parameters
                bar_error_kw = dict(ecolor='k', lw=1, capsize=2, capthick=1)

                # plot the data as a barchart
                barcontainer_norm = barnorm_ax.bar(box_indices,
                                                   df_all_exp_nonredundant["mean{}{}".format(d,norm_dataset)],
                                                   yerr=df_all_exp_nonredundant["SEM{}{}".format(d,norm_dataset)],
                                                   align="center",
                                                   error_kw=bar_error_kw, color = '#1b9e77')
                # set the xticks
                barnorm_ax.set_xticks(box_indices)

                # set the limits of the x-axis
                barnorm_ax.set_xlim([-1, bar_x_n_boxes])
                # set the limit of the y-axis
                barnorm_ax.set_ylim(0)
                # set the y-axis title
                # bar_ax.set_ylabel("EC50 (ug/ml)")
                # set the ylabel extension string
                if norm_dataset == "_nonnorm":
                    dose_units  = settings["x-axis (dose) units"]
                elif norm_dataset == "_norm":
                    samplenames_selected = df_samplenames.loc[df_allp.longname, :]
                    # find the sample longname that is marked "True" as a standard
                    standard_name = samplenames_selected[samplenames_selected["standard for normalisation?"] == True].index[0]
                    # check if the standard for normalisation has a short name
                    if standard_name in samplenames_selected.index:
                        standard_name_short = samplenames_selected.loc[standard_name, "short name"]
                    else:
                        # use the long name
                        standard_name_short = standard_name
                    # create a string to extend the y-axis label
                    dose_units = "% {}".format(standard_name_short)
                else:
                    raise TypeError("dataset for standardisation{} is not recognised".format(norm_dataset))

                ylabel_str = "{a}{b}, {c} ({d})".format(a=settings["calculation_type"],
                                      b=str(settings["percentage_response"]),
                                      c=settings["x-axis (dose) label"],
                                      d=dose_units)
                barnorm_ax.set_ylabel(ylabel_str)

                for nametype in list_nametypes:
                    # define name on the x-axis
                    barnorm_x_names = df_for_barchart[nametype]
                    # set the labels of the x-bar_axis
                    barnorm_ax.set_xticklabels(barnorm_x_names, rotation=90)

                    # ax.annotate(s="%s%s" % (nametype,d), xy=(0.015,0.93), fontsize=af, xycoords=xyc)
                    barnorm_ax.set_title("analysed data ({e} experiments),  {b}".format(b=d_name,
                                                             e=dff.loc[dff["run analysis"] == True].shape[0]))
                    # automatically tighten the layout and save figure
                    barnorm_fig.tight_layout()
                    # save the figure
                    barnorm_fig.savefig(analysed_data_basename + "_bar_" + nametype + d_name + norm_dataset + '.png',
                                        format='png', dpi=150)
                plt.close('all')

                # move sample names to index
                df_all_exp_nonredundant["sample_number"] = df_all_exp_nonredundant.index
                df_all_exp_nonredundant.set_index("longname", inplace=True)
                # save to excel
                df_all_exp_nonredundant.to_excel(writer, sheet_name = "EC50_mean{}{}".format(d_name,norm_dataset))
                df_all_exp_nonredundant.to_csv(analysed_data_basename + "_EC50_mean{}{}".format(d_name,norm_dataset) + ".csv",
                                               sep=",", quoting=csv.QUOTE_NONNUMERIC)
        writer.save()
        writer.close()


        #
        # '''
        # AllExperimentsPlots_01-04: Comparison EC50 different experiments. (longnames & shortnames, original and ful)
        # '''
        # # plt.rcParams['legend.numpoints'] = 3
        # # # sort index
        # # df_allp.sort_index(inplace=True)
        # # # copy the index with the full sample names to to a new column
        # # df_allp["longname"] = df_allp.index
        # # # create a column with the shortened sample names
        # # df_allp["shortname"] = df_allp.longname.apply(lambda x : samplenames_dict[x] if x in list(samplenames_dict.keys()) else x)
        # # reorder the index according to desired order in the settings file
        # # # create a dictionary of the sample order
        # # sampleorder_dict = dict(zip(df_samplenames.index, df_samplenames["order in figure"]))
        # # # create a new column with the preferred sample order
        # # df_allp["sampleorder"] = df_allp.longname.apply(lambda x : sampleorder_dict[x] if x in list(sampleorder_dict.keys()) else x)
        # # df_allp.sort_values(by="sampleorder", inplace=True)
        #
        # # # replace index with a simple range of integers
        # # df_allp.index = range(df_allp.shape[0])
        # # # give the index a name
        # # df_allp.index.name = "sSnum"
        # # # sort the columns
        # # df_allp.sort_index(axis=1, inplace=True)
        # for d in datasets:
        #     # identify the columns associated with that dataset
        #     col_contains_d = list(pd.Series(df_allp.columns).apply(lambda x: d in x))
        #     # select only data for that dataset (e.g. only orig data)
        #     sel_df_allp = df_allp.loc[:,col_contains_d]
        #     # set up a boolean to determine if any data is correctly normalised to the standard
        #     at_least_one_dataset_is_normalised = False
        #     # create separate figures for the long names, and the short names
        #     for name in ["longname", "shortname"]:
        #         conduct_normalisation, list_norm_datasets, \
        #         at_least_one_dataset_is_normalised = setup_normalisation(df_samplenames,
        #                                                                  df_allp,
        #                                                                  at_least_one_dataset_is_normalised)
        #
        #         samplenames_selected = df_samplenames.loc[df_allp.longname, :]
        #
        #         # if True in list(samplenames_selected["standard for normalisation?"]):
        #         #     n_standards_labelled_as_True = samplenames_selected["standard for normalisation?"].value_counts()[True]
        #         # if n_standards_labelled_as_True == 0:
        #         #     conduct_normalisation = False
        #         # elif n_standards_labelled_as_True == 1:
        #         #     conduct_normalisation = True
        #         #     at_least_one_dataset_is_normalised = True
        #         # elif n_standards_labelled_as_True > 1:
        #         #     conduct_normalisation = False
        #         #     raise ValueError("multiple samples are labelled as standards for normalisation. "
        #         #                      "Please check the samplenames tab of your settings file.")
        #         #
        #         # if "long name" in df_samplenames.columns:
        #         #     df_samplenames.set_index("long name", inplace = True)
        #         #
        #         # samplenames_selected = df_samplenames.loc[df_allp.longname, :]
        #         # if True in list(samplenames_selected["standard for normalisation?"]):
        #         #     n_standards_labelled_as_True = samplenames_selected["standard for normalisation?"].value_counts()[True]
        #         # if n_standards_labelled_as_True == 0:
        #         #     conduct_normalisation = False
        #         # elif n_standards_labelled_as_True == 1:
        #         #     conduct_normalisation = True
        #         # elif n_standards_labelled_as_True > 1:
        #         #     conduct_normalisation = False
        #         #     raise ValueError("multiple samples are labelled as standards for normalisation. "
        #         #                      "Please check the samplenames tab of your settings file.")
        #         #
        #         # if conduct_normalisation:
        #         #     if at_least_one_dataset_is_normalised:
        #         #         list_norm_datasets = ["_nonnorm", "_norm"]
        #         #     else:
        #         #         list_norm_datasets = ["_nonnorm"]
        #         # else:
        #         #     list_norm_datasets = ["_nonnorm"]
        #
        #         for norm_dataset in list_norm_datasets:
        #
        #             # create a new figure with a single plot
        #             plt.close("all")
        #             scat_fig, scat_ax = plt.subplots()
        #             # fig_norm, ax_norm = plt.subplots()
        #
        #             # iterate through the columns in the dataframe, each representing a single experiment
        #             for exp_nr, c in enumerate(sel_df_allp.columns):
        #                 data_colour = t20[exp_nr]
        #                 sys.stdout.write(".")
        #                 sys.stdout.flush()
        #                 # convert any stringlists lists
        #                 df_allp.loc[:,c] = df_allp.loc[:,c].apply(lambda x: eval(x) if isinstance(x,str) else x)
        #                 # convert any lists of strings to lists of floats
        #                 df_allp.loc[:,c] = df_allp.loc[:,c].apply(lambda x: [float(s) for s in x] if isinstance(x,list) else x)
        #                 # convert any lists of floats to lists of numpy arrays
        #                 df_allp.loc[:,c] = df_allp.loc[:,c].apply(lambda x: np.array(x) if isinstance(x,list) else x)
        #
        #                 if norm_dataset == "_nonnorm":
        #                     series_EC50_data = df_allp.loc[:,c].dropna()
        #                 elif norm_dataset == "_norm":
        #                     # find the sample longname that is marked "True" as a standard
        #                     standard_name = samplenames_selected[samplenames_selected["standard for normalisation?"] == True].index[0]
        #                     # find the sSNum (index number) of the standard
        #                     standard_sSnum = df_allp[df_allp.longname == standard_name].index[0]
        #                     # calculated the average EC50 for the standard
        #                     standard_EC50 = np.mean(df_allp.loc[standard_sSnum, c])
        #                     # calculate the normalised EC50 as a percentage (EC50/standard*100)
        #                     series_EC50_data = df_allp.loc[:,c]/standard_EC50*100
        #                     series_EC50_data.dropna(inplace=True)
        #
        #                 #######################################################
        #                 #    add_floats_arrays_to_scatterplot
        #                 #######################################################
        #                 # set the transparency
        #                 alpha = 0.5
        #                 # create a series of bools, describing if the data is a float
        #                 is_float_series = series_EC50_data.apply(lambda x: isinstance(x,float))
        #                 # use boolean series to select the float data
        #                 values_single_sample = series_EC50_data.loc[is_float_series]
        #                 # create a series of bools, describing if the data is a list
        #                 is_array_series = series_EC50_data.apply(lambda x: isinstance(x,np.ndarray))
        #                 # if some of the datapoints are lists (multiple sample replicates in a single experiment)
        #                 if True in list(is_array_series):
        #                     # use boolean series to select the data with lists
        #                     values_mult_sample = series_EC50_data.loc[is_array_series]
        #                     # the multiple EC50 values belong to a single sample number. Identify list of sample numbers.
        #                     index_ssNum_with_lists = values_mult_sample.index.tolist()
        #                     # create a list of xvalues for the scattergram, for multiple samples per day
        #                     list_xvalues_from_mult = []
        #                     # create a list of yvalues for the scattergram, for multiple samples per day
        #                     list_yvalues_from_mult = []
        #                     # iterate through each selected sample number, associated with multiple EC50 values
        #                     for sSnum in index_ssNum_with_lists:
        #                         list_EC50_values = list(values_mult_sample[sSnum])
        #                         # create corresponding x-axis sample numbers (array of sSnum of length len(values_mult_sample))
        #                         list_index_values = list(np.ones(len(list_EC50_values))*sSnum)
        #                         # add x and y values to the lists for all samples
        #                         list_xvalues_from_mult.extend(list_index_values)
        #                         list_yvalues_from_mult.extend(list_EC50_values)
        #
        #                 xvalues_from_single = list(values_single_sample.index)
        #                 yvalues_from_single = list(values_single_sample)
        #
        #
        #                 xvalues_all_single_exp = xvalues_from_single + list_xvalues_from_mult
        #                 yvalues_all_single_exp = yvalues_from_single + list_yvalues_from_mult
        #
        #                 # plot the float data as a scattergram (x-axis is the range, resembling a bar or line chart)
        #                 scat_ax.scatter(xvalues_all_single_exp, yvalues_all_single_exp, color=data_colour, s=40, alpha=alpha, label=c[:-5])
        #
        #             #######################################################
        #             #    format_and_save_analysis_scatterplot
        #             #######################################################
        #             # define the x-axis labels (long or short)
        #             xticklabels = list(df_allp[name])
        #             settings_name = os.path.split(settings_excel_file)[1]
        #
        #             # set the transparency
        #             alpha = 0.5
        #             # set the xticks and labels to match the index of df_allp
        #             scat_ax.set_xticks(np.arange(df_allp.shape[0]))
        #             scat_ax.set_xticklabels(xticklabels, rotation=90)
        #             # set the grid to go in between the sample names, as minor xticks
        #             scat_ax.set_xticks(np.arange(df_allp.shape[0])+0.5, minor=True)
        #             scat_ax.grid(which='minor', alpha=0.9)
        #             # set the x axis limits
        #             scat_ax.set_xlim(-0.5, df_allp.shape[0])
        #
        #             # set the y-axis title
        #             scat_ax.set_ylabel("{a}{b}, {c} ({d}){e}".format(a=settings["calculation_type"],
        #                                   b=str(settings["percentage_response"]),
        #                                   c=settings["x-axis (dose) label"],
        #                                   d=settings["x-axis (dose) units"],
        #                                   e=dose_units))
        #             # add legend and title
        #             scat_ax.set_title("analysed data ({e} experiments),  "
        #                          "{a}{b}{f},  {c}".format(a=name,b=d,f=norm_dataset,
        #                                                   c=settings_name,e=sel_df_allp.shape[1]))
        #             lgd = scat_ax.legend(loc='upper center', bbox_to_anchor=(0.5,-0.1), ncol=2, scatterpoints=1, numpoints =1)
        #             # automatically tighten the layout and save figure
        #             scat_fig.tight_layout()
        #             scat_fig.savefig(analysed_data_basename + "_datapoints_" + name + d_name + norm_dataset + '.png', format='png', dpi=300, bbox_extra_artists=(lgd,), bbox_inches='tight')
        #
        #         plt.close("all")

        print('\nCollection and analysis of data from multiple experiments is finished.\n'
              'Output files are saved in the following directory:\n{}'.format(outpath))

        return df_allp, dff
    else:
        print("\nNo files are selected for analysis! Double-check TRUE/FALSE columns in settings file.")
        return "no files selected", "no files selected"

def setup_normalisation(df_samplenames, df_allp):
    """ Sets up the normalisation of data for that experiment relative to a control.

    Checks to see that a single sample is labelled TRUE for "standard for normalisation?" in the settings file.

    Parameters
    ----------
    df_samplenames : pandas DataFrame
        Dataframe containing the samplenames and column identifying the sample for normalisation.
    df_allp : pandas Dataframe
        Dataframe for all individual EC50 datapoints from all experiments marked as "TRUE" in the run_analysis
        column of the settings file. Includes individual values for samples found twice in a single experiment.
        index : unique sample numbers (sSnum)
        columns: experiment names (e.g. "Mon_23.05", "Tue_24.05", "Wed_25.05")
        values : EC50 values, or lists of EC50 values as strings

    Returns
    -------
    conduct_normalisation : boolean
        If True, all necessary requirements for normalisation OF THAT DATASET are found, and it should proceed.
    list_norm_datasets : list
        Suffix for the datasets, e.g. ["_nonnorm", "_norm"] if conduct_normalisation is True
        If conduct_normalisation is False, will contain a suffix for a single, non-normalised dataset, ["_nonnorm"]
    """

    # if any of the sample names with valid EC50 data extracted from output files are listed in the settings file
    dataset_contains_samples_in_settings_file = bool(set(df_allp.longname).intersection(set(df_samplenames.index)))
    if dataset_contains_samples_in_settings_file:
        # select the relevant sample names in the settings file
        samplenames_selected = df_samplenames.loc[df_allp.longname, :]
        # if any of the samples are labelled as a standard
        if True in list(samplenames_selected["standard for normalisation?"]):
            # determine number of samples labelled True as the "standard for normalisation?"
            n_standards_labelled_as_True = samplenames_selected["standard for normalisation?"].value_counts()[True]
            # if a single sample is labelled as the standard
            if n_standards_labelled_as_True == 1:
                # this particular experiment can be normalised to a standard
                conduct_normalisation = True
                # prepare suffixes for the filenames that will include both non-normalised and normalised data
                list_norm_datasets = ["_nonnorm", "_norm"]
            elif n_standards_labelled_as_True > 1:
                raise ValueError("Multiple samples are labelled as standards for normalisation. "
                                 "Only a single sample can be labelled as a standard."
                                 "Please check the samplenames tab of your settings file.")
        else:
            conduct_normalisation = False
            list_norm_datasets = ["_nonnorm"]
    else:
        conduct_normalisation = False
        list_norm_datasets = ["_nonnorm"]
        print("Note: None of the samples with valid EC50 values was found in the samples tab of the settings file. "
              "Sample ordering and normalisation to a standard will not be conducted.")
    return conduct_normalisation, list_norm_datasets

def compare_rawdata(settings_excel_file, sample_names, **kwargs):
    """ Compare raw dose-response curves between selected samples, for selected experiments.

    Processes the datafiles in settings excel file marked as "TRUE" for "run analysis."
    Collects output from the "run_curvefit" program, but only for the selected samples.
    Recreates the fitted curves from the four-parameter Hill equation, with the previously calculated hill_constants.

    The output of the compare_rawdata is saved in the same folder as the "run_analysis":
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
    datasets = ast.literal_eval(settings["adjust.datasets"])

    # create boolean
    at_least_one_sample_found_in_selected_datafiles = False

    if not isinstance(sample_names[0], tuple):
        # the list of sample names contains only strings, and therefore only a single raw analysis is performed
        # convert to a list containing only one tuple, with the sample names for comparison
        sample_names = [tuple(sample_names)]

    for s_num, sample_tuple in enumerate(sample_names):
        if True in list(dff.loc[:, "run analysis"]):
            n_files_to_analyse = dff.loc[dff["run analysis"] == True].shape[0]
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
                for fn in dff.loc[dff["run analysis"] == True].index:
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
                ax.set_xlabel(settings["x-axis (dose) label"])
                ax.set_ylabel(settings["y-axis (response) label"],rotation='vertical')
                # save the figure in png format
                figpath = os.path.join(outpath, "{b}_{n}{d}.png".format(b=basename,n=fig_name,d=d_name))
                fig.savefig(figpath, format = "png", dpi = 150)
                plt.close("all")
    print('\nComparison of raw data is finished.\nOutput files are saved in the following directory:\n{}'.format(outpath))