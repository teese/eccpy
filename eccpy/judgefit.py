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
import ast
import eccpy.tools as tools

def judge_fit(dfe, sLet, df_settings):
    ''' Analyses the fit of the experimental data to a sigmoidal curve.

    Parameters
    ----------
    dfe : pandas dataframe
        Dataframe for evaluation of a single curve.
        The index is a long list of parameter names, e.g. EC50, hill_constants, etc.
        There are three columns. A) sLet (the relevant sample letter),
                                 B) sLet_okay (column for labeling if data seems okay)
                                 C) sLet_colour (column for defining text colour in output figures (red for low-quality)
    sLet : string
        Sample letter.
    df_settings : pandas DataFrame
        User-defined settings derived from the settings excel sheet.
        Includes all the relevant user-defined thresholds for the judge_fit scripts.

    Returns
    -------
    dfe : pandas dataframe
        Original dataframe, with added annotations in sLet_okay and sLet_colour columns concerning data quality.

    Notes
    -----
    How are low-quality EC50 values identified?
    1) Does the hillslope parameter from the fitted curve approach 0?
         Indicates an exponential rather than a sigmoid curve.
    2) Does the curve begin or end with a negative value on the y-axis?
        Indicates an exponential curve rather than a sharp Z-shaped curve
    3) Is the r_squared value below the user-defined threshold?
        Indicates a poor fit to data.
    4) Are the number of datapoints before (or after) the EC50 below a user-defined threshold?
        Suggests that the dose range is not optimal, curve may not be sigmoidal.
    5) Does the standard deviation of the datapoints to the left or right of the EC50 exceed a user-defined value?
        Suggests a poor fit to data, or the presence of an outlier.
    6) Does the Slope At X-axis Extremes (SAXE) exceed a user-defined threshold?
        High slopes at the initial or final datapoints suggest non-sigmoidal curves, usually due to an
        inappropriate range of doses used in that experiment.
    7) Does the dose stepsize at the EC50 exceed a user-defined threshold?
        The higher the dose steps, the less accurate the final EC50 value.

    What happens when a low-quality EC50 value or curve is identified?
    1) Each of the filters has a neighbouring "_okay" or "_colour" parameter
        If the EC50 value seems to be okay, for that particular filter
            "_okay" will be labelled True
            "_colour" will be labelled as "k" (black)
        If the EC50 value seems to be of low quality
            "_okay" will be labelled False
            "_colour" will be labelled as "r" (red)
    2) There will be a final parameter for each EC50 value, "data_seems_okay"
        if there are no False values in "_okay" column for all tested filters :
            data_seems_okay = True
            the EC50 will be coloured black in the output graph with the curve, and in the final barchart
            with all data for that day/experiment
        else if there is at least one False value in the "_okay" column :
            data_seems_okay = False
            The EC50 will be coloured red in the output graph with the curve, and in the final barchart
                with all data for that day/experiment.
            The EC50 value WILL BE IGNORED by the "analysis" scripts that compare results for different days/experiments.

    '''
    # setup cutoffs for judging data quality
    # number datapoints neighbouring the EC50 that are excluded from the highdose and lowdose data selection
    # set higher if you use a large number of dose concentrations
    n_neighb = df_settings.loc["n_neighb","B"]
    # maximum standard deviation of the response datapoints at high dose concentration
    max_std_resp_highdose_dp = df_settings.loc["max_std_resp_highdose_dp","B"]
    # maximum standard deviation of the response datapoints at low dose concentration
    max_std_resp_lowdose_dp = df_settings.loc["max_std_resp_lowdose_dp","B"]
    min_flat_lowdose_dp = df_settings.loc["min_flat_lowdose_dp","B"]
    min_flat_highdose_dp = df_settings.loc["min_flat_highdose_dp","B"]

    # minimum rsquared of the fit from sigmoidal curve to the data
    min_rsquared = df_settings.loc["min_rsquared","B"]
    # minimum acceptable dase concentration stepsizes. Smaller stepsizes give more accurate EC50 values!
    max_acceptable_doseconc_stepsize_at_EC50 = df_settings.loc["max_acceptable_doseconc_stepsize_at_EC50","B"]
    max_recommended_doseconc_stepsize_at_EC50 = df_settings.loc["max_recommended_doseconc_stepsize_at_EC50","B"]
    # minimum hillslope of the fit from sigmoidal curve to the data (below 1, tends not to be sigmoidal)
    weak_hillslope_range = ast.literal_eval(df_settings.loc["weak_hillslope_range","B"])
    # minimum value for the end of the curve, on the y-axis (below -1, tends not to be sigmoidal)
    min_curve_lowresp = df_settings.loc["min_curve_lowresp","B"]

    # create a list that contains the database suffixes (_orig for original, _ful for fixed upper limit)
    # datasets = ["_orig", "_ful"]
    datasets = ast.literal_eval(df_settings.loc["adjust.datasets", "B"])
    for d in datasets:
        x = np.array(dfe.loc["x{}".format(d), sLet])
        y = np.array(dfe.loc["y{}".format(d), sLet])
        # identify the datapoints at high dose concentrations
        dfe.loc["indices_highdose_datapoints{}".format(d),sLet] = np.where(x > dfe.loc["EC50{}".format(d), sLet])[0]
        # remove the datapoint closest to EC50
        dfe.loc["indices_highdose_datapoints_excl_nearest_EC50{}".format(d),sLet] = dfe.loc["indices_highdose_datapoints{}".format(d),sLet][n_neighb:]
        # slice using the indices to yield the OD600 values for the highdose datapoints
        dfe.loc["response_highdose_datapoints{}".format(d),sLet] = y[dfe.loc["indices_highdose_datapoints_excl_nearest_EC50{}".format(d),sLet]]
        # count the number of highdose datapoint
        dfe.loc["n_highdose_datapoints{}".format(d),sLet] = len(dfe.loc["response_highdose_datapoints{}".format(d),sLet])

        # identify the lowdose datapoints, count and measure standard deviation
        # identify the lowdose datapoints (x < EC50)
        dfe.loc["indices_lowdose_datapoints{}".format(d),sLet] = np.where(x < dfe.loc["EC50{}".format(d), sLet])[0]
        # exclude datapoint closest to the EC50
        dfe.loc["indices_lowdose_datapoints_excl_nearest_EC50{}".format(d),sLet] = dfe.loc["indices_lowdose_datapoints{}".format(d),sLet][:-n_neighb]
        # use index to select the y-axis (response) data
        dfe.loc["response_lowdose_datapoints{}".format(d),sLet] = y[dfe.loc["indices_lowdose_datapoints_excl_nearest_EC50{}".format(d),sLet]]
        # count the datapoints
        dfe.loc["n_lowdose_datapoints{}".format(d),sLet] = len(dfe.loc["response_lowdose_datapoints{}".format(d),sLet])

        # indices_highdose_datapoints_excl_nearest_EC50_orig = indices_highdose_datapoints_orig[1:]
        # response_highdose_datapoints_orig = y_orig[indices_highdose_datapoints_excl_nearest_EC50_orig]
        # # count the number of highdose datapoints
        # dfe.loc["n_highdose_datapoints_orig",sLet] = len(response_highdose_datapoints_orig)

        # # identify the ful datapoints at high dose concentrations, ignoring datapoint closest to EC50
        # indices_highdose_datapoints_ful = np.where(x_orig > EC50_ful)[0]
        # indices_highdose_datapoints_excl_nearest_EC50_ful = indices_highdose_datapoints_ful[1:]
        # response_highdose_datapoints_ful = y_orig[indices_highdose_datapoints_excl_nearest_EC50_ful]
        # # count the number of highdose datapoints
        # dfe.loc["n_highdose_datapoints_ful",sLet] = len(response_highdose_datapoints_ful)

        #######################################################################################################
        #                                                                                                     #
        #        Are the number of datapoints before (or after) the EC50 below a user-defined threshold?      #
        #                                                                                                     #
        #######################################################################################################

        # judge whether the data contains enough high and lowdose datapoints
        if dfe.loc["n_highdose_datapoints{}".format(d),sLet] >= min_flat_highdose_dp:
            dfe.loc["n_highdose_datapoints{}".format(d),"%s_okay" % sLet] = True
            dfe.loc["n_highdose_datapoints{}".format(d),"%s_colour" % sLet] = 'k'
        else:
            dfe.loc["n_highdose_datapoints{}".format(d),"%s_okay" % sLet] = False
            dfe.loc["n_highdose_datapoints{}".format(d),"%s_colour" % sLet] = 'r'

        # evaluate as "okay" if number of highdose or lowdose datapoints is more than two
        if dfe.loc["n_lowdose_datapoints{}".format(d),sLet] >= min_flat_lowdose_dp:
            dfe.loc["n_lowdose_datapoints{}".format(d),"%s_okay" % sLet] = True
            dfe.loc["n_lowdose_datapoints{}".format(d),"%s_colour" % sLet] = 'k'
        else:
            dfe.loc["n_lowdose_datapoints{}".format(d),"%s_okay" % sLet] = False
            dfe.loc["n_lowdose_datapoints{}".format(d),"%s_colour" % sLet] = 'r'

        #############################################################################################################
        #                                                                                                           #
        #  Does the standard deviation of datapoints to the left or right of the EC50 exceed a user-defined value?  #
        #                                                                                                           #
        #############################################################################################################

        # judge whether the standard deviation of the high and lowdose datapoints is acceptable
        if dfe.loc["n_highdose_datapoints{}".format(d),sLet] > 1:
            # calculate std of highdose datapoints
            dfe.loc["std_resp_highdose_datapoints{}".format(d),sLet] = np.std(dfe.loc["response_highdose_datapoints{}".format(d),sLet])
            # evaluate as "okay" if std of highdose datapoints is less than a cutoff value (max_std_resp_highdose_dp)
            if dfe.loc["std_resp_highdose_datapoints{}".format(d),sLet] < max_std_resp_highdose_dp:
                dfe.loc["std_resp_highdose_datapoints{}".format(d),"%s_okay" % sLet] = True
                dfe.loc["std_resp_highdose_datapoints{}".format(d),"%s_colour" % sLet] = 'k'
            else:
                dfe.loc["std_resp_highdose_datapoints{}".format(d),"%s_okay" % sLet] = False
                dfe.loc["std_resp_highdose_datapoints{}".format(d),"%s_colour" % sLet] = 'r'
        else:
            # there is either insufficient lowresponse or highresponse datapoints(insuff_lowresp_dp, or insuff_highresp_dp).
            # Replace std with 0, and colour black on the figure.
            dfe.loc["std_resp_highdose_datapoints{}".format(d),sLet] = 0
            dfe.loc["std_resp_highdose_datapoints{}".format(d),"%s_colour" % sLet] = 'k'

        if dfe.loc["n_lowdose_datapoints{}".format(d),sLet] > 1:
            # calculate std of lowdose datapoints
            dfe.loc["std_resp_lowdose_datapoints{}".format(d),sLet] = np.std(dfe.loc["response_lowdose_datapoints{}".format(d),sLet])
            # evaluate as "okay" if std of lowdose datapoints is less than a cutoff value
            if dfe.loc["std_resp_lowdose_datapoints{}".format(d),sLet] < max_std_resp_lowdose_dp:
                dfe.loc["std_resp_lowdose_datapoints{}".format(d),"%s_okay" % sLet] = True
                dfe.loc["std_resp_lowdose_datapoints{}".format(d),"%s_colour" % sLet] = 'k'
            else:
                dfe.loc["std_resp_lowdose_datapoints{}".format(d),"%s_okay" % sLet] = False
                dfe.loc["std_resp_lowdose_datapoints{}".format(d),"%s_colour" % sLet] = 'r'
        else:
            # there is either insufficient lowresponse or highresponse datapoints(insuff_lowresp_dp, or insuff_highresp_dp).
            # Replace std with 0, and colour black.
            dfe.loc["std_resp_lowdose_datapoints{}".format(d),sLet] = 0
            dfe.loc["std_resp_lowdose_datapoints{}".format(d),"%s_colour" % sLet] = 'k'

        #############################################################################################################
        #                                                                                                           #
        #                    Does the dose stepsize at the EC50 exceed a user-defined threshold?                    #
        #                                                                                                           #
        #############################################################################################################

        # identify the tested dose concentration below the EC50
        indices_lowdose_datapoints = np.where(x < dfe.loc["EC50{}".format(d),sLet])[0]
        if indices_lowdose_datapoints.size != 0:
            doseconc_before_EC50 = x[indices_lowdose_datapoints[-1]]
            # identify the tested dose concentration after the EC50
            doseconc_after_EC50 = x[dfe.loc["indices_highdose_datapoints{}".format(d),sLet][0]]
            # add values to output dataframe, so that the plot can be annotated
            dfe.loc["doseconc_steps_at_EC50{}".format(d),sLet] = (doseconc_before_EC50, doseconc_after_EC50)
            # calculate the stepsize at the EC50. Smaller is better!
            doseconc_stepsize_at_EC50 = doseconc_after_EC50 - doseconc_before_EC50
            dfe.loc["doseconc_stepsize_at_EC50{}".format(d),sLet] = doseconc_stepsize_at_EC50
            # evaluate as "okay" if the stepsize at the EC50 is smaller than the min acceptable value
            if doseconc_stepsize_at_EC50 <= max_acceptable_doseconc_stepsize_at_EC50:
                dfe.loc["doseconc_stepsize_at_EC50{}".format(d),"%s_okay" % sLet] = True
                # if the stepsize is small, colour to dark red as a warning that the doseconc should be optimised
                if doseconc_stepsize_at_EC50 <= max_recommended_doseconc_stepsize_at_EC50:
                    dfe.loc["doseconc_stepsize_at_EC50{}".format(d),"%s_colour" % sLet] = 'k'
                else:
                    # colour dark red
                    dfe.loc["doseconc_stepsize_at_EC50{}".format(d),"%s_colour" % sLet] = '#990033'
            else:
                # the stepsize is extremely high, and the data therefore has little value. doseconc needs to be optimised.
                dfe.loc["doseconc_stepsize_at_EC50{}".format(d),"%s_okay" % sLet] = False
                dfe.loc["doseconc_stepsize_at_EC50{}".format(d),"%s_colour" % sLet] = 'r'
        else:
            # there is either insufficient lowresponse or highresponse datapoints(insuff_lowresp_dp, or insuff_highresp_dp).
            # Stepsize can't be calculated. Replace with 0, and colour grey.
            dfe.loc["doseconc_stepsize_at_EC50{}".format(d),"%s_okay" % sLet] = False
            dfe.loc["doseconc_steps_at_EC50{}".format(d),sLet] = (0,0)
            dfe.loc["doseconc_stepsize_at_EC50{}".format(d),sLet] = 0
            dfe.loc["doseconc_stepsize_at_EC50{}".format(d),"%s_colour" % sLet] = '0.5'

        #############################################################################################################
        #                                                                                                           #
        #                          Is the r_squared value below the user-defined threshold?                         #
        #                                                                                                           #
        #############################################################################################################

        # evaluate rsquared_orig as okay if above a fixed limit (min_rsquared)
        if dfe.loc["rsquared{}".format(d),sLet] > min_rsquared:
            dfe.loc["rsquared{}".format(d),"%s_okay" % sLet] = True
            dfe.loc["rsquared{}".format(d),"%s_colour" % sLet] = 'k'
        else:
            dfe.loc["rsquared{}".format(d),"%s_okay" % sLet] = False
            dfe.loc["rsquared{}".format(d),"%s_colour" % sLet] = 'r'

        #############################################################################################################
        #                                                                                                           #
        #                      Does the hillslope parameter from the fitted curve approach 0?                       #
        #                                                                                                           #
        #############################################################################################################

        # evaluate slope as okay if it is outside the weak_hillslope_range, which is usually nonsigmoidal
        hillslope = dfe.loc["hillslope{}".format(d), sLet]
        if weak_hillslope_range[0] < hillslope < weak_hillslope_range[1]:
            # if it's outside the range, label as "data_needs_checking"
            dfe.loc["hillslope{}".format(d),"%s_okay" % sLet] = False
            dfe.loc["hillslope{}".format(d),"%s_colour" % sLet] = 'r'
        else:
            # if it's outside the range, label as okay
            dfe.loc["hillslope{}".format(d),"%s_okay" % sLet] = True
            dfe.loc["hillslope{}".format(d),"%s_colour" % sLet] = 'k'

        #############################################################################################################
        #                                                                                                           #
        #                   Does the curve begin or end with a negative value on the y-axis?                        #
        #                                                                                                           #
        #############################################################################################################

        if dfe.loc["curve_min_norm{}".format(d),sLet] > min_curve_lowresp:
            dfe.loc["curve_min_norm{}".format(d),"%s_okay" % sLet] = True
            dfe.loc["curve_min_norm{}".format(d),"%s_colour" % sLet] = 'k'
        else:
            dfe.loc["curve_min_norm{}".format(d),"%s_okay" % sLet] = False
            dfe.loc["curve_min_norm{}".format(d),"%s_colour" % sLet] = 'r'

        ##############################################################################################################
        #                                                                                                            #
        #                  Does the Slope At X-axis Extremes (SAXE) exceed a user-defined threshold?                 #
        #                  1) calculate average stepsize                                                             #
        #                  2) define points left and right of the xaxis extremes                                     #
        #                  3) calculate slope between points                                                         #
        #                  4) determine if slope is shallow enough to indicate an S- or Z-shaped curve               #
        #                                                                                                            #
        ##############################################################################################################

        hill_constants = dfe.loc["hill_constants{}".format(d),sLet]
        xnorm = dfe.loc["xnorm{}".format(d), sLet]
        # find the average stepsize by creating two truncated arrays from xnorm, and subtracting
        xnorm_left = np.array(list(xnorm)[:-1])
        xnorm_right = np.array(list(xnorm)[1:])
        xnorm_stepsizes = xnorm_right - xnorm_left
        xnorm_stepsize_mean = xnorm_stepsizes.mean()
        # define the width surrounding the datapoint for the slope measurement
        # calculated as the mean stepsize multiplied by a user value (0.001 to 1.0)
        width_lowdose_slope = xnorm_stepsize_mean/2 * df_settings.loc["width_lowdose_slope","B"]
        width_highdose_slope = xnorm_stepsize_mean/2 * df_settings.loc["width_highdose_slope","B"]
        # define the min and max of normalised datapoints (will simply be 0 and 1 for normalised data)
        xnorm_min, xnorm_max = xnorm.min(), xnorm.max()
        # define SAXE lowdose/highdose x-axis datapoints (to the left and right of the original datapoints)
        saxe_lowdose_x_dp_left = xnorm_min - width_lowdose_slope
        # if it is negative (results in nan in the sigmoidal function), replace with 0
        saxe_lowdose_x_dp_left = saxe_lowdose_x_dp_left if saxe_lowdose_x_dp_left > 0 else 0
        saxe_lowdose_x_dp_right = xnorm_min + width_lowdose_slope
        saxe_highdose_x_dp_left = xnorm_max - width_highdose_slope
        saxe_highdose_x_dp_right = xnorm_max + width_highdose_slope
        # calculate the y-values on the curve, for the x-values surrounding the min and max datapoints
        saxe_lowdose_y_dp_left = tools.hill_eq(hill_constants, saxe_lowdose_x_dp_left)
        saxe_lowdose_y_dp_right = tools.hill_eq(hill_constants, saxe_lowdose_x_dp_right)
        saxe_highdose_y_dp_left = tools.hill_eq(hill_constants, saxe_highdose_x_dp_left)
        saxe_highdose_y_dp_right = tools.hill_eq(hill_constants, saxe_highdose_x_dp_right)
        # calculate the linear slope (y2 - y1)/(x2 - x1) between the chosen datapoints from the curve
        saxe_lowdose_1 = (saxe_lowdose_y_dp_right - saxe_lowdose_y_dp_left)/(saxe_lowdose_x_dp_right - saxe_lowdose_x_dp_left)
        saxe_highdose_1 = (saxe_highdose_y_dp_right - saxe_highdose_y_dp_left)/(saxe_highdose_x_dp_right - saxe_highdose_x_dp_left)
        # convert slopes to positive numbers
        saxe_lowdose = abs(saxe_lowdose_1)
        saxe_highdose = abs(saxe_highdose_1)
        # add to output dataframe, dfe
        dfe.loc["saxe_lowdose{}".format(d), sLet] = saxe_lowdose
        dfe.loc["saxe_highdose{}".format(d), sLet] = saxe_highdose
        dfe.loc["saxe_lowdose_values{}".format(d), sLet] = [[saxe_lowdose_x_dp_left, saxe_lowdose_x_dp_right],
                                                             [saxe_lowdose_y_dp_left, saxe_lowdose_y_dp_right]]
        dfe.loc["saxe_highdose_values{}".format(d), sLet] = [[saxe_highdose_x_dp_left, saxe_highdose_x_dp_right],
                                                            [saxe_highdose_y_dp_left, saxe_highdose_y_dp_right]]

        # obtain the max allowed values for the slopes
        max_lowdose_slope = df_settings.loc["max_lowdose_slope","B"]
        max_highdose_slope = df_settings.loc["max_highdose_slope","B"]

        # check that the calculated slopes of the curve do not exceed the saxe_max_slope
        if saxe_lowdose < max_lowdose_slope:
            dfe.loc["saxe_lowdose{}".format(d),"%s_okay" % sLet] = True
            dfe.loc["saxe_lowdose{}".format(d),"%s_colour" % sLet] = 'k'
        else:
            dfe.loc["saxe_lowdose{}".format(d),"%s_okay" % sLet] = False
            dfe.loc["saxe_lowdose{}".format(d),"%s_colour" % sLet] = 'r'

        if saxe_highdose < max_highdose_slope:
            dfe.loc["saxe_highdose{}".format(d),"%s_okay" % sLet] = True
            dfe.loc["saxe_highdose{}".format(d),"%s_colour" % sLet] = 'k'
        else:
            dfe.loc["saxe_highdose{}".format(d),"%s_okay" % sLet] = False
            dfe.loc["saxe_highdose{}".format(d),"%s_colour" % sLet] = 'r'

    # #######################################################################################################
    # #                                                                                                     #
    # #                          Calculate the slope at the EC50 (not currently in use)                     #
    # #                                                                                                     #
    # #######################################################################################################
    #
    # # calculate the slope surrounding the EC50
    # x1 = dfe.loc["EC50_norm_bq{}".format(d),sLet]
    # x2 = x1 + 0.01
    # y1 = tools.hill_eq(hill_constants, x1)
    # y2 = tools.hill_eq(hill_constants, x2)
    # middle_slope = (y2 - y1)/(x2 - x1)
    # # print("middle_slope:", middle_slope)

    return dfe