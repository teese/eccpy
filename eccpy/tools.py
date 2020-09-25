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
import os
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
import csv

import eccpy


def aaa(df_or_series):
    """ Function for use in debugging.
    Saves pandas Series or Dataframes to a user-defined csv file.
    """
     # convert any series to dataframe
    if isinstance(df_or_series, pd.Series):
        df_or_series = df_or_series.to_frame()
    csv_out = r"D:\data\000_aaa_temp_df_out.csv"
    df_or_series.to_csv(csv_out, sep=",", quoting=csv.QUOTE_NONNUMERIC)

def hill_eq(hill_constants, x):
    """ Four parameter sigmoidal Hill equation.

    y = upper + (lower-upper)/(1+(x/EC50)**-hillslope)

    Parameters
    ----------
    hill_constants : tuple
        Tuple of the four parameters : upper, lower, EC50 and hillslope
    x : float
        x-value for use in the equation

    Constants
    ---------
    upper : float
        Renamed from "bottom" in some sources.
        Will approach the minimum response (perhaps zero) in the dose-response S-curve.
        Will approach the maximum response in an inverse LD50 curve.
        Is not currently used as a parameter to judge whether curve is sigmoidal.
    lower : float
        Renamed from "top" in some sources.
        Will approach the maximum response in the dose-response S-curve.
        Will approach zero in an inverse LD50 curve.
        Is not currently used as a parameter to judge whether curve is sigmoidal.
    EC50 : float
        Does not always accurately reflect the EC50 from the data.
        Is not currently used as a parameter to judge whether curve is sigmoidal.
        Is not currently used as the calculated EC50 value.
    hillslope : float
        A high hillslope generally indicates a steep curve, an accurate EC50 calculation, and a strong sigmoidal shape.
        Note that the hillslope can be strongly negative.
        Hillslope values approaching zero (e.g. -1 > hill_slope > +1) generally indicate an exponential, rather than
        a sigmoidal curve. For this reason, hillslope is currently used as a parameter to judge whether the curve is
        sigmoidal, using the parameters chosen in the settings_excel_file and implemented in the "judge_fit" program.

    Notes
    -----
    The hill equation is used for the EC50 calculation, but there are several variants of the formula.

    Dose-Response Equations
    variant 1: http://www.graphpad.com/guides/prism/6/curve-fitting/index.htm?reg_classic_dr_variable.htm
    variant 2: http://en.wikipedia.org/wiki/EC50
    variant 3: http://pydoc.net/Python/cgptoolbox/0.1.2/cgp.sigmoidmodels.doseresponse/

    Sigmoid equations:
    variant 4: y = c * 1.0 / (1.0 + ((k/x)**g))
    variant 5: y = c / (1 + np.exp(-k*(x-x0))) + y0
        (http://stackoverflow.com/questions/7588371/scipy-leastsq-goodness-of-fit-estimator)

    This equation uses variant 2. See Wikipedia for a more detailed explanation.

    In theory, any formula that accurately fits a sigmoidal curve is appropriate for the purpose. The Hill equation
    is preferred, because under some circumstances the EC50 value is one of the fitted parameters in the hill_constants
    tuple. However this value is extremely unreliable, and regularly gives EC50 values that are out of the range of the
    datapoints. A more accurate EC50 calculation is achieved through root finding using the brent equation.

    """
    upper, lower, EC50, hillslope = hill_constants
    y = upper + (lower-upper)/(1+(x/EC50)**-hillslope)
    return y

def residuals(constants, function, x, y):
    """
    Function used to optimise the fit of the curve to the data.
    It calculates the distance between y-value from real data and y-value from the function (sigmoid/sine/etc).
    """
    return y - function(constants, x)

def hill_eq_brentq(xvalues_for_curve, hill_constants, y_value_curve_center):
    """ Residual function for the four parameter sigmoidal Hill equation.

    For further detail on the Hill equation, see the relevant docstring for hill_eq.

    y = hill_eq(x) - y_value_curve_center

    Parameters
    ----------
    xvalues_for_curve : array
        Numpy array of >250 x-values between the lowest and highest dose.
    hill_constants : tuple
        Tuple of the four parameters : upper, lower, EC50 and hillslope
    y_value_curve_center : array
        Represents the y-value of the curve. The difference between this and the real y-value from experimental
        data (e.g. y_orig_norm) is used to optimise the fit of the curve to the data.

    Returns
    -------
    y - y_value_curve_center : array
        Array that approaches zero, when an optimum fit is found.
    """
    upper, lower, EC50, hillslope = hill_constants
    y = upper + (lower-upper)/(1+(xvalues_for_curve/EC50)**-hillslope)
    return y - y_value_curve_center

def normalise_0_1(arraylike):
    """ Normalise an array to values between 0 and 1.

    The following linear formula is used.
    norm_array = (orig_array - array_min)/(array_max - array_min)

    The use of this simple linear formula allows the normalised data to be "denormalised" later, so long as
    the min and max values of the original array are known.

    Parameters
    ----------
    arraylike : array
        Numpy array (or other arraylike) dataset of floats or ints to be normalised.

    Returns
    -------
    normalised : array
        Array of floats, containing the normalised datapoints.
    array_min : float
        Minimum value of the original data. Necessary in order to "denormalise" the data later, back to the effective
        original values.
    array_max : float
        Maximum value of the original data. Necessary in order to "denormalise" the data later, back to the effective
        original values.

    Usage
    -----
    normalised_array, min_, max_ = normalise_0_1(original_array)
    # or, if denormalisation is not necessary
    normalised_array = normalise_0_1(original_array)[0]
    # for further usage examples, see the docstring for denormalise_0_1
    """
    array_min = np.min(arraylike)
    array_max = np.max(arraylike)
    normalised = (arraylike - array_min)/(array_max - array_min)
    # convert to float
    normalised = np.array(normalised).astype(float)
    return normalised, array_min, array_max

def denormalise_0_1(value_or_array, array_min, array_max):
    """ Denormalise a value or array to orig values.

    For use after normalisation between 0 and 1 with the normalise_0_1 function.

    The normalisation formula (normalise_0_1):
        norm_array = (orig_array - array_min)/(array_max - array_min)

    The denormalisation formula (denormalise_0_1):
        denormalised_array = norm_array*(array_max - array_min) + array_min

    Parameters
    ----------
    value_or_array : int, float or arraylike
        Int or float to be denormalised.
        Numpy array (or other arraylike) of data (float, int, etc) to be denormalised.

    Returns
    -------
    normalised : float, or numpy array
        Array of floats, containing the normalised datapoints.
    array_min : float
        Minimum value of the original data. Necessary in order to "denormalise" the data later, back to the effective
        original values.
    array_max : float
        Maximum value of the original data. Necessary in order to "denormalise" the data later, back to the effective
        original values.

    Usage
    -----
    from eccpy.tools import normalise_0_1, denormalise_0_1
    import numpy as np
    original_array = np.linspace(10,130,10)
    original_array[2], original_array[4] = 3, 140
    print(original_array)
    # normalise original array
    normalised_array, min_, max_ = normalise_0_1(original_array)
    print(normalised_array)
    # do stuff to normalised array (e.g., multiply by 0.5)
    normalised_array_halved = normalised_array * 0.5
    # denormalise values to match the equivalents in the original array.
    # Note that the min value (3) was normalised to zero, and was therefore not affected by multiplication.
    normalised_array_halved_denorm = denormalise_0_1(normalised_array_halved, min_, max_)
    print(normalised_array_halved_denorm)
    # now calculate average values, and check that they match
    norm_array_mean = np.mean(normalised_array)
    norm_array_mean_denormalised = denormalise_0_1(norm_array_mean, min_, max_)
    orig_array_mean = np.mean(original_array)
    # print the two mean values. They should be equal.
    print(norm_array_mean_denormalised)
    print(orig_array_mean)
    """
    if isinstance(value_or_array, list):
        raise ValueError('this function accepts arraylike data, not a list. '
                         'Please check data or convert list to numpy array')
    elif isinstance(value_or_array, float):
        #print("found a float")
        denormalised = value_or_array*(array_max - array_min) + array_min
    elif isinstance(value_or_array, np.ndarray):
        #print("found an array")
        denormalised = value_or_array*(array_max - array_min) + array_min
    elif isinstance(value_or_array, pd.Series):
        #print("found a series")
        denormalised = value_or_array*(array_max - array_min) + array_min
    else:
        print("Unknown datatype. denormalise_0_1 has been given an input that does not appear to be "
              "an int, float, np.ndarray or pandas Series\n"
              "Attempting to process as if it is arraylike.....")
    return denormalised


def normalise_between_2_values(arraylike, min_value, max_value, invert=False):
    """Normalises an array of data between two desired values.

    Any values below min_value will be converted to 0.
    Any values above max_value will be converted to 1.
    Optionally, the normalised array can be inverted, so that the original highest
    values are 0, and the original lowest values are now 1.

    Parameters
    ----------
    arraylike : np.ndarray
        Arraylike original data (numpy array or pandas Series)
    min_value : float
        Desired minimum value for normalisation
    max_value : float
        Desired max value for normalisation
    invert : bool
        If True, normalised data will be inverted (former highest value = 0)

    Returns
    -------
    normalised : np.ndarray
        Normalised array of data

    Usage
    -----
    from eccpy.tools import normalise_between_2_values
    # for array
    orig_array = np.array(range(0, 15))
    norm_array = normalise_between_2_values(orig_array, 3, 10)
    # for pandas Dataframe
    df["norm_data"] = normalise_between_2_values(df["orig_data"], 3, 10)
    """
    # normalise array between min and max values
    normalised = (arraylike - min_value)/(max_value - min_value)
    # replace anything above 1 with 1
    normalised[normalised > 1] = 1
    # replace anything below 0 with 0
    normalised[normalised < 0] = 0
    # if desired, invert the normalised values
    if invert:
        normalised = abs(normalised - 1)
    return normalised


def setup_t20_colour_list():
    """ Setup a list of colours for the figures.

    Returns
    -------
    t20 : list
        List of RGB colour tuples, normalised between 0 and 1 as used by python colours.
    """
    """ Setup colours for the figures
    return : t20, a list of colour tuples
    """
    colour_lists_dict = {}
    #define colour lists
    colour_lists = {
                    'tableau20' : [
                                 (31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
                                 (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
                                 (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
                                 (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
                                 (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)
                                    ],
                    'tableau20blind' : [
                                         (0, 107, 164), (255, 128, 14), (171, 171, 171), (89, 89, 89),
                                         (95, 158, 209), (200, 82, 0), (137, 137, 137), (163, 200, 236),
                                         (255, 188, 121), (207, 207, 207)
                                          ]
                    }
    #normalise the colours for the colour lists
    for rgb_list in colour_lists:
        colour_array = np.array(colour_lists[rgb_list])/255.
        colour_array_tup = tuple(map(tuple,colour_array))
        colour_lists[rgb_list] = colour_array_tup
        #add normalised colours to output dictionary
        colour_lists_dict[rgb_list] = colour_lists[rgb_list]
    # extract the tableau20 lists, join together
    t20 = list(colour_lists_dict["tableau20"] + colour_lists_dict["tableau20blind"])
    # extend the list in case someone need a long list of colours (after 30 colours, they will be redundant!)
    t20 = t20 + t20 + t20 + t20 + t20
    return t20

def reindex_df_so_selected_cols_are_first(df, selected_cols, raise_error=True):
    """ Reindex DataFrame so that selected columns are displayed first (on the left)

    Parameters
    ----------
    df : pandas DataFrame
        Pandas 2D DataFrame with unique columns
    selected_cols : list
        List of strings representing the column names to place first

    Returns
    -------
    df : pandas DataFrame
        Original DataFrame with altered column order
    """

    # convert columns to list
    col_list_orig = list(df.columns)
    # remove the list_cols_to_place_first from the original columns
    for col in selected_cols:
        if col in col_list_orig:
            col_list_orig.remove(col)
        else:
            if raise_error == True:
                raise ValueError("\n\nError, reindex_df_so_selected_cols_are_first, '%s' not in columns" % col)
            else:
                print("\n\nError, reindex_df_so_selected_cols_are_first, '%s' not in columns" % col)
    # join to create desired list of columns, and reindex the dataframe
    col_list_final = selected_cols + col_list_orig
    return df.reindex(columns = col_list_final)

def convert_listlike_cols_to_str(df, list_cols, convert_nan=False, nanreplacement = "[]"):
    """ Convert listlike values in pandas DataFrame to stringlists.

    Writing a DataFrame to excel and csv raises errors due to presence of listlike or arraylike data.
    Here, all listlike and arraylike data is converted to a "stringlist"

    Parameters
    ----------
    df : pandas DataFrame
        Pandas DataFrame where some values are listlike or arraylike
    list_cols : list
        List of columns that contain listlike/arraylike should be converted to stringlists
    convert_nan : bool
        Whether np.nan values will be converted
    nanreplacement : string
        String to replace np.nan with. Default is a string representing an empty list. "[]"

    Returns
    ----------
    df : pandas DataFrame
        DataFrame with listlike and arraylike converted to stringlists

    Note:
    ----------
    # Convert individual stringlist back to a numpy array as follows
    np.array(ast.literal_eval(stringlist))
    # convert a column of stringlists back to arrays as follows
    df.loc[:,"x"] = df.loc[:,"x"].apply(lambda x : np.array(ast.literal_eval(x)))
    """
    for col in list_cols:
        if col in df:
            # convert each np.array to a list, and then a stringlist
            series = df[col].dropna()
            # check if the series is empty (no data)
            if series.empty == False:
                example_data1 = series[0]
                # if the datatype is a numpy array or pandas series, convert to a list
                if "ndarray" in str(type(example_data1)) or "Series" in str(type(example_data1)):
                    df[col] = series.apply(lambda x: list(x))
                example_data2 = df[col].dropna()[0]
                # check that the first nonnan datapoint is now a list
                if "list" in str(type(example_data2)):
                    df[col] = df[col].dropna().apply(lambda x: str(x))
                else:
                    raise TypeError("datatype for col {a}, ({b}) is not listlike".format(a=col, b=str(type(example_data2))))
                # if desired, convert np.nan to empty stringlists
                if convert_nan == True:
                    df[col] = df[col].fillna(nanreplacement)
            else:
                #do nothing. There is no listlike data in the column, because all values are empty
                pass
        else:
            raise KeyError("The column {} is not in the dataframe".format(col))
    return df

class DatafileError(Exception):
    """
    Custom error, to be used when there is a problem with the datafiles
    """
    pass

class DataMismatchError(Exception):
    """
    Custom error, to be used when there is a mismatch in the data somehow, such as an unequal number of
    dose and response concentrations, or the Contains_Data columns/rows do not match up.
    """
    pass

def create_dict_organising_subplots(n_plots_per_fig,n_rows):
    '''
    Function to help organise the creation of figures that contain multiple plots.

    Parameters
    ----------
    n_plots_per_fig : int
        Number of plots per figure in total
    n_rows : int
        Number of rows of plots per figure

    Notes
    -----
    For example, 15 histograms printed in figures with 8 histograms per figure/page.
    Returns a dict that gives a tuple for each plot/graph.
    newfig, savefig, fig_nr, plot_nr_in_fig, r, c
    r and c are used to index pyplot subplots as follows
    fig, axarr = plt.subplots(2,2)
    plotcontainer = axarr[r,c].plot(x, y)
    '''
    dict_organising_subplots = {}
    #figure number
    fig_nr = 0
    #plot number in figure
    plot_nr_in_fig = 0
    #row number in figure
    r = 0
    #column number in figure
    c = 0
    #whether the figure needs to be saved
    savefig = False
    #whether a new figure needs to be created
    newfig = True

    for plotnr in range(1, 500):
        #add current counters to dict
        dict_organising_subplots[plotnr] = (newfig, savefig, fig_nr, plot_nr_in_fig, r, c)
        plot_nr_in_fig += 1
        r += 1
        newfig = False
        savefig = False
        #if plot_nr_in_fig is the last one before the new figure, then savefig = True
        if plot_nr_in_fig % (n_plots_per_fig - 1) == 0 and plot_nr_in_fig != 0:
            savefig = True
        #if plot_nr_in_fig is in a multiple of n_rows, then the plot goes to the second column
        if plot_nr_in_fig % n_rows == 0 and plot_nr_in_fig != 0:
            c += 1
            r = 0
        #if the plotnr is in a multiple of n_plots_per_fig, then a new figure needs to created, and everything else reset
        if plotnr % n_plots_per_fig == 0 and plotnr != 0:
            #go to second figure
            fig_nr += 1
            #reset values
            plot_nr_in_fig = 0
            r = 0
            c = 0
            newfig = True
    return dict_organising_subplots

def convert_truelike_to_bool(input_item, convert_float=False, convert_nontrue=True):
    """Converts true-like values ("true", 1, True", "WAHR", etc) to python boolean True.

    Parameters
    ----------
    input_item : string or int
        Item to be converted to bool (e.g. "true", 1, "WAHR" or the equivalent in several languagues)
    convert_float: bool
        Convert floats to bool.
        If True, "1.0" will be converted to True
    convert_nontrue : bool
        If True, the output for input_item not recognised as "True" will be False.
        If True, the output for input_item not recognised as "True" will be the original input_item.

    Returns
    -------
    return_value : True, or input_item
        If input_item is True-like, returns python bool True. Otherwise, returns the input_item.

    Usage
    -----
    # convert a single value or string
    convert_truelike_to_bool("true")
    # convert a column in a pandas DataFrame
    df["column_name"] = df["column_name"].apply(convert_truelike_to_bool)
    """
    list_True_items = [True, 'True', "true","TRUE","T","t",'wahr', 'WAHR', 'prawdziwy', 'verdadeiro', 'sann', 'istinit',
                       'veritable', 'Pravda', 'sandt', 'vrai', 'igaz', 'veru', 'verdadero', 'sant', 'gwir', 'PRAWDZIWY',
                       'VERDADEIRO', 'SANN', 'ISTINIT', 'VERITABLE', 'PRAVDA', 'SANDT', 'VRAI', 'IGAZ', 'VERU',
                       'VERDADERO', 'SANT', 'GWIR', 'bloody oath', 'BLOODY OATH', 'nu', 'NU','damn right','DAMN RIGHT']

    # if you want to accept 1.0 as a true value, add it to the list
    if convert_float:
        list_True_items += [1.0,"1.0"]
    # check if the user input string is in the list_True_items
    input_item_is_true = input_item in list_True_items
    # if you want to convert non-True values to "False", then nontrue_return_value = False
    if convert_nontrue:
        nontrue_return_value = False
    else:
        # otherwise, for strings not in the True list, the original string will be returned
        nontrue_return_value = input_item
    # return True if the input item is in the list. If not, return either False, or the original input_item
    return_value = input_item_is_true if input_item_is_true == True else nontrue_return_value

    return return_value

def convert_nonelike_to_none(input_item):
    """Converts None-like values ("none", "NULL", None, etc) to the uniform string "None".

    Note, the output is NOT the python None, but a string.

    Parameters
    ----------
    input_item : string or int
        Item to be converted to None (e.g. "none", "NULL" or the equivalent in several languagues)

    Returns
    -------
    return_value : string
        If input_item is None-like, returns python string "None". Otherwise, returns the input_item.

    Usage
    -------
    # convert a single value or string
    convert_nonelike_to_none("none")
    # convert a column in a pandas DataFrame
    df["column_name"] = df["column_name"].apply(convert_nonelike_to_none)
    """
    list_None_items = [None, "none","NONE","null","NULL",'Nijedna', 'Cap', 'Niti', 'Ingen', 'Geen', 'Aucun',
                       'Keine', 'Okenn', 'Egyik', 'Tidak', 'Nessuno', 'Hakuna', 'pagh', 'Neviens', 'Tiada', 'L-eda',
                       'Mix', 'Ingen', 'Ninguno', 'Brak', 'Nenhum', 'Nici', 'Niko', 'Nobena', 'Ninguno', 'Ingen',
                       'Dim','NIJEDNA', 'CAP', 'NITI', 'INGEN', 'GEEN', 'AUCUN', 'KEINE', 'OKENN', 'EGYIK', 'TIDAK',
                       'NESSUNO', 'HAKUNA', 'PAGH', 'NEVIENS', 'TIADA', 'L-EDA', 'MIX', 'INGEN', 'NINGUNO', 'BRAK',
                       'NENHUM', 'NICI', 'NIKO', 'NOBENA', 'NINGUNO', 'INGEN', 'DIM']
    # determine if input_item is in the list
    input_item_is_None = input_item in list_None_items
    # define the return_value as either the string "None" or the original input item
    return_value = "None" if input_item_is_None == True else input_item

    return return_value

def format_cols_2digit(df, skip_last_col=True):
    """Formats a dataframes columns so that numbers are always two-digits (padded with 0)

    Parameters
    ----------
    df : pandas DataFrame
        Input DataFrame.
    skip_last_col : bool
        A special case, where the very last column contains text, rather than a number, and should be excluded.

    Returns
    -------
    reformatted_cols : list
        The DataFrame columns, reformatted with padded 0 to make 2 digits.

    """
    if skip_last_col:
        # reformat the columns to be padded stringnumbers. (exclude the last "Contains_Data")
        reformatted_cols = ["%02d" % col for col in df.columns[:-1]]
        # add last column back to list
        reformatted_cols.append(df.columns[-1])
    else:
        # reformat the columns to be padded stringnumbers. (exclude the last "Contains_Data")
        reformatted_cols = ["%02d" % col for col in df.columns]

    return reformatted_cols


def assert_df_contains_no_nan_values(df: pd.DataFrame):

    df_contains_nan = df.isnull().any().any()

    if df_contains_nan:
        positions_with_nan: List[str] = []

        for column_name in df.columns:
            for index_name in df.index:
                value = df.at[index_name, column_name]
                if pd.isnull(value):
                    position_string = f"(row='{index_name}', column='{column_name}')"
                    positions_with_nan.append(position_string)

        all_positions_str = ", and ".join(positions_with_nan)

        raise ValueError(f"\n\nThe 'files' tab of the excel settings appears to contain empty cells at {all_positions_str}. "
                         "Please delete any partially filled rows. "
                         "Fill empty cells with the text 'None' if necessary.")


def get_eccpy_module_path()-> Path:
    eccpy_module_path = Path(os.path.abspath(eccpy.__file__)).parents[1]
    #eccpy_module_path = Path(__file__).parents[1]
    return eccpy_module_path