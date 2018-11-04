# This will aim to produce a figure from the SP output for the camp data
# the data will be split up by the three species and then by the three location types
# these are mangrove, reef_one and reef_two. I will call the reefs reef_far and reef_near as one is much
# closer than the other to the mangrove site. I have recently produced some figures for teh TARA work which would
# have a similar layout to what we want for this. Ideally I would like to code this so that it can generically be
# used to work with an input that is the excel file, and the seq and type outputs and then split up according
# to the parameters provided. This can then be reused for additional studies.

# read in the excel and generate a df from it that will be an info df used for the figure
# Create a pandas df from the data_sheet if it was provided

import pandas as pd
import sys
import random
import matplotlib as mpl
mpl.use('TkAgg')

from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
# from matplotlib.pyplot import *
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
import numpy as np

def create_info_df(data_sheet_path, param_one, param_two):
    # param one and two are the parameters (column names) that we want to be dividing the samples by
    info_df = pd.read_excel(io=data_sheet_path, header=1, index_col=0)

    # get list of the param_one categories
    param_one_list = list(set(info_df.loc[:, param_one].values.tolist()))

    # get list of the param_two categories
    param_two_list = list(set(info_df.loc[:, param_two].values.tolist()))

    return info_df, param_one_list, param_two_list


def create_parameter_matrix_fig(data_sheet_path, param_one, param_two):


    # this will return the info_df and the list of categories for param one and two
    info_df, param_one_list, param_two_list = create_info_df(
        data_sheet_path = data_sheet_path, param_one = param_one, param_two = param_two
    )

    # here we define the number of columns and rows we want in each of the legends (one div, one type)
    # number of cols for the plot will be the number of param_one categories
    max_n_leg_cols_div = 5
    # number of rows will the the number of param_two categories plus one for the legend
    max_n_leg_rows_div = 8
    # number of legend cells
    num_leg_cells_div = max_n_leg_cols_div * max_n_leg_rows_div

    # number of cols for the plot will be the number of param_one categories
    max_n_leg_cols_type = 4
    # number of rows will the the number of param_two categories plus one for the legend
    max_n_leg_rows_type = 8
    # number of legend cells
    num_leg_cells_type = max_n_leg_cols_type * max_n_leg_rows_type

    # read in the SymPortal relative seq abundance output and create df for plotting
    path_to_tab_delim_count_DIV = '/Users/humebc/Google_Drive/projects/camp_gardener/sp_outputs_041118/38_DBV_011118_2018-11-04_02-36-01.981594.DIVs.relative.txt'
    smp_id_to_smp_name_dict, smp_name_to_smp_id_dict, sp_output_df_div, colour_dict_div, ordered_list_of_seqs = process_div_df(path_to_tab_delim_count_DIV=path_to_tab_delim_count_DIV, num_leg_cells=num_leg_cells_div)

    # read in the SymPortal relative type abundance output and create df for plotting
    path_to_tab_delim_count_type = '/Users/humebc/Google_Drive/projects/camp_gardener/sp_outputs_041118/38_DBV_011118_2018-11-04_02-36-01.981594.profiles.relative.txt'
    colour_dict_type, sp_output_df_type, sorted_type_prof_names_by_local_abund,  = process_type_df(
        path_to_tab_delim_count_type = path_to_tab_delim_count_type, num_leg_cells=num_leg_cells_type)


    ordered_sample_list = sp_output_df_type.index.values.tolist()


    sp_output_df_div = sp_output_df_div[ordered_list_of_seqs]

    # SETUP AXES
    num_ax_cols = len(param_one_list)
    num_ax_rows = len(param_two_list) + 1
    ax_list, leg_axes = setup_axes(num_ax_cols, num_ax_rows)


    # Here we have the main data axes popualated, the blank spacer and the legend axes
    # later we can make it so that the axes have a constant height. same with the spacer

    plot_data_axes(param_one_list, param_two_list, ax_list, colour_dict_div, colour_dict_type, info_df, ordered_sample_list,
                   smp_id_to_smp_name_dict,
                   smp_name_to_smp_id_dict, sp_output_df_div, sp_output_df_type, param_one, param_two)

    #
    # PLOT DIV LEGEND
    plot_div_legend(colour_dict_div, leg_axes, max_n_leg_cols_div, max_n_leg_rows_div, num_leg_cells_div, ordered_list_of_seqs)

    # PLOT TYPE LEGEND
    plot_type_legend(colour_dict_type, leg_axes, max_n_leg_cols_type, max_n_leg_rows_type, num_leg_cells_type,
                     sorted_type_prof_names_by_local_abund)
    apples = 'pears'
    #
    # # add the labels and text here so that we don't have to debug through all of the plotting each time
    # add_labels(ax_list, leg_axes, extra_ax_list)

def plot_type_legend(colour_dict_type, leg_axes, max_n_cols_type, max_n_rows_type, num_leg_cells_type,
                     sorted_type_prof_names_by_local_abund, string_cut_off=10):
    # Since the matplotlib legends are pretty rubbish when made automatically, I vote that we make our own axes
    # all in favour... Ok.
    # Let's plot the boxes and text that are going to make up the legend in another subplot that we will put underneath
    # the one we currenty have. So.. we will add a subplot when we initially create the figure. We will make the axis
    # 100 by 100 just to make our coordinate easy to work with. We can get rid of all of the axes lines and ticks
    # The type names are generally quite long so we will cut the type legends down to 4 x 8
    # we should start plotting in the top left working right and then down
    # until we have completed 100 sequences.
    # Y axis coordinates
    # we will allow a buffer of 0.5 of the legend box's height between each legend box.
    # as such the coordinates of each y will be in increments of 100 / (1.5 * num rows)
    # the depth of the Rectangle for the legend box will be 2/3 * the above.
    y_coord_increments = 100 / (max_n_rows_type)
    leg_box_depth = 2 / 3 * y_coord_increments
    # X axis coordinates
    # for the x axis we will work in sets of three columns were the first col will be for the box
    # and the second and third cols will be for the text
    # as such the x coordinates will be in increments of 100 / (3 * numcols) starting with 0
    # the width of the legend Rectangle will be the above number * 1/6 (I am making this smaller for the types).
    x_coord_increments = 100 / max_n_cols_type
    leg_box_width = x_coord_increments / 6
    # go column by column
    # we can now calculate the actual number of columns and rows we are going to need.
    if len(sorted_type_prof_names_by_local_abund) < num_leg_cells_type:
        if len(sorted_type_prof_names_by_local_abund) % max_n_cols_type != 0:
            n_rows_type = int(len(sorted_type_prof_names_by_local_abund) / max_n_cols_type) + 1
        else:
            n_rows_type = int(len(sorted_type_prof_names_by_local_abund) / max_n_cols_type)
        last_row_len = len(sorted_type_prof_names_by_local_abund) % max_n_cols_type
    else:
        n_rows_type = max_n_rows_type
        last_row_len = max_n_cols_type
    its2_profile_count = 0
    # Once we know the number of rows, we can also adjust the y axis limits
    leg_axes[1].set_xlim(0, 100)
    # axarr[-1].set_ylim(0, 100)
    leg_axes[1].set_ylim(0, ((n_rows_type - 1) * y_coord_increments) + leg_box_depth)
    leg_axes[1].invert_yaxis()
    # If there are more sequences than there are rows x cols then we need to make sure that we are only going
    # to plot the first row x cols number of sequences.
    sys.stdout.write(
        '\nGenerating figure legend for {} most common sequences\n'.format(str(max_n_rows_type * max_n_cols_type)))

    for row_increment in range(min(n_rows_type, max_n_rows_type)):
        # if not in the last row then do a full set of columns
        if row_increment + 1 != n_rows_type:
            for col_increment in range(max_n_cols_type):
                # add the legend Rectangle
                leg_box_x = col_increment * x_coord_increments
                leg_box_y = row_increment * y_coord_increments
                leg_axes[1].add_patch(Rectangle((leg_box_x, leg_box_y),
                                                width=leg_box_width, height=leg_box_depth,
                                                color=colour_dict_type[
                                                    sorted_type_prof_names_by_local_abund[its2_profile_count]]))

                # add the text
                text_x = leg_box_x + leg_box_width + (0.2 * leg_box_width)
                text_y = leg_box_y + (0.5 * leg_box_depth)
                # lets limit the name to 15 characters and '...'
                if len(sorted_type_prof_names_by_local_abund[its2_profile_count]) > string_cut_off:
                    text_for_legend = sorted_type_prof_names_by_local_abund[its2_profile_count][
                                      :string_cut_off] + '...'
                else:
                    text_for_legend = sorted_type_prof_names_by_local_abund[its2_profile_count]
                leg_axes[1].text(text_x, text_y, text_for_legend,
                                 verticalalignment='center',
                                 fontsize=8)

                # increase the sequence count
                its2_profile_count += 1
        # else just do up to the number of last_row_cols
        else:
            for col_increment in range(last_row_len):
                # add the legend Rectangle
                leg_box_x = col_increment * x_coord_increments
                leg_box_y = row_increment * y_coord_increments
                leg_axes[1].add_patch(Rectangle((leg_box_x, leg_box_y),
                                                width=leg_box_width, height=leg_box_depth,
                                                color=colour_dict_type[
                                                    sorted_type_prof_names_by_local_abund[its2_profile_count]]))

                # add the text
                text_x = leg_box_x + leg_box_width + (0.2 * leg_box_width)
                text_y = leg_box_y + (0.5 * leg_box_depth)
                # lets limit the name to 15 characters and '...'
                if len(sorted_type_prof_names_by_local_abund[its2_profile_count]) > string_cut_off:
                    text_for_legend = sorted_type_prof_names_by_local_abund[its2_profile_count][
                                      :string_cut_off] + '...'
                else:
                    text_for_legend = sorted_type_prof_names_by_local_abund[its2_profile_count]
                leg_axes[1].text(text_x, text_y, text_for_legend,
                                 verticalalignment='center',
                                 fontsize=8)

                # Increase the sequences count
                its2_profile_count += 1
    remove_axes_but_allow_labels(leg_axes[1])

def plot_div_legend(colour_dict_div, leg_axes, max_n_cols_div, max_n_rows_div, num_leg_cells_div, ordered_list_of_seqs):
    # Since the matplotlib legends are pretty rubbish when made automatically, I vote that we make our own axes
    # all in favour... Ok.
    # Let's plot the boxes and text that are going to make up the legend in another subplot that we will put underneath
    # the one we currenty have. So.. we will add a subplot when we initially create the figure. We will make the axis
    # 100 by 100 just to make our coordinate easy to work with. We can get rid of all of the axes lines and ticks
    # lets aim to plot a 10 by 10 legend max
    # we should start plotting in the top left working right and then down
    # until we have completed 100 sequences.
    # Y axis coordinates
    # we will allow a buffer of 0.5 of the legend box's height between each legend box.
    # as such the coordinates of each y will be in increments of 100 / (1.5 * num rows)
    # the depth of the Rectangle for the legend box will be 2/3 * the above.
    y_coord_increments = 100 / (max_n_rows_div)
    leg_box_depth = 2 / 3 * y_coord_increments
    # X axis coordinates
    # for the x axis we will work in sets of three columns were the first col will be for the box
    # and the second and third cols will be for the text
    # as such the x coordinates will be in increments of 100 / (3 * numcols) starting with 0
    # the width of the legend Rectangle will be the above number * 1/3.
    x_coord_increments = 100 / max_n_cols_div
    leg_box_width = x_coord_increments / 3
    # go column by column
    # we can now calculate the actual number of columns and rows we are going to need.
    if len(ordered_list_of_seqs) < num_leg_cells_div:
        if len(ordered_list_of_seqs) % max_n_cols_div != 0:
            n_rows_div = int(len(ordered_list_of_seqs) / max_n_cols_div) + 1
        else:
            n_rows_div = int(len(ordered_list_of_seqs) / max_n_cols_div)
        last_row_len = len(ordered_list_of_seqs) % max_n_cols_div
    else:
        n_rows_div = max_n_rows_div
        last_row_len = max_n_cols_div
    sequence_count = 0
    # Once we know the number of rows, we can also adjust the y axis limits
    leg_axes[0].set_xlim(0, 100)
    # axarr[-1].set_ylim(0, 100)
    leg_axes[0].set_ylim(0, ((n_rows_div - 1) * y_coord_increments) + leg_box_depth)
    leg_axes[0].invert_yaxis()
    # If there are more sequences than there are rows x cols then we need to make sure that we are only going
    # to plot the first row x cols number of sequences.
    sys.stdout.write(
        '\nGenerating figure legend for {} most common sequences\n'.format(str(max_n_rows_div * max_n_cols_div)))
    for row_increment in range(min(n_rows_div, max_n_rows_div)):
        # if not in the last row then do a full set of columns
        if row_increment + 1 != n_rows_div:
            for col_increment in range(max_n_cols_div):
                # add the legend Rectangle
                leg_box_x = col_increment * x_coord_increments
                leg_box_y = row_increment * y_coord_increments
                leg_axes[0].add_patch(Rectangle((leg_box_x, leg_box_y),
                                                width=leg_box_width, height=leg_box_depth,
                                                color=colour_dict_div[ordered_list_of_seqs[sequence_count]]))

                # add the text
                text_x = leg_box_x + leg_box_width + (0.2 * leg_box_width)
                text_y = leg_box_y + (0.5 * leg_box_depth)
                leg_axes[0].text(text_x, text_y, ordered_list_of_seqs[sequence_count], verticalalignment='center',
                                 fontsize=8)

                # increase the sequence count
                sequence_count += 1
        # else just do up to the number of last_row_cols
        else:
            for col_increment in range(last_row_len):
                # add the legend Rectangle
                leg_box_x = col_increment * x_coord_increments
                leg_box_y = row_increment * y_coord_increments
                leg_axes[0].add_patch(Rectangle((leg_box_x, leg_box_y),
                                                width=leg_box_width, height=leg_box_depth,
                                                color=colour_dict_div[ordered_list_of_seqs[sequence_count]]))

                # add the text
                text_x = leg_box_x + leg_box_width + (0.2 * leg_box_width)
                text_y = leg_box_y + (0.5 * leg_box_depth)
                leg_axes[0].text(text_x, text_y, ordered_list_of_seqs[sequence_count], verticalalignment='center',
                                 fontsize=8)

                # Increase the sequences count
                sequence_count += 1
    remove_axes_but_allow_labels(leg_axes[0])

def setup_axes(max_n_cols, max_n_rows):
    # https://matplotlib.org/users/gridspec.html
    fig = plt.figure(figsize=(14, 10))
    # the bottom row will be for the legend
    # the second to last will just be invisible to give a space between the legend and the other plots
    # we also want to include a gridspec plot after each of the main three. These will hold the csw and surface
    # samples
    # we will put in an invisible spacer so as an extra row
    gs = plt.GridSpec(max_n_rows + 1, max_n_cols)
    ax_list = []
    # first make an ax for each of the main data plots
    for row in range(max_n_rows - 1):
        for col in range(max_n_cols):
            ax = plt.Subplot(fig, gs[row, col])
            ax_list.append(ax)
            fig.add_subplot(ax)
    # now make the ax for the space
    blank_ax = plt.subplot(gs[max_n_rows - 1, :])
    # make the axes invisible for the space ax
    remove_axes_but_allow_labels(blank_ax)
    # now split up the final row to put the legend in. One for DIVs and one for TYPEs
    temp_grid_spec_subplot_leg = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[-1, :])
    leg_axes = []
    for i in range(2):
        ax = plt.Subplot(fig, temp_grid_spec_subplot_leg[i])
        leg_axes.append(ax)
        fig.add_subplot(ax)
    return ax_list, leg_axes


def process_div_df(path_to_tab_delim_count_DIV, num_leg_cells):
    sp_output_df = pd.read_csv(path_to_tab_delim_count_DIV, sep='\t', lineterminator='\n', header=0, index_col=0)

    # In order to be able to drop the DIV row at the end and the meta information rows, we should
    # drop all rows that are after the DIV column. We will pass in an index value to the .drop
    # that is called here. To do this we need to work out which index we are working with
    index_values_as_list = sp_output_df.index.values.tolist()
    for i in range(-1, -(len(index_values_as_list)), -1):
        if index_values_as_list[i].startswith('DIV'):
            # then this is the index (in negative notation) that we need to cut from
            meta_index_to_cut_from = i
            break
    sp_output_df = sp_output_df.iloc[:meta_index_to_cut_from]

    # create sample id to sample name dict
    smp_id_to_smp_name_dict = {ID: '_'.join(nm.split('_')[:3]) for ID, nm in
                               zip(sp_output_df.index.values.tolist(), sp_output_df['sample_name'].values.tolist())}
    smp_name_to_smp_id_dict = {'_'.join(nm.split('_')[:3]): ID for ID, nm in
                               zip(sp_output_df.index.values.tolist(), sp_output_df['sample_name'].values.tolist())}

    # now lets drop the QC columns from the SP output df and also drop the clade summation columns
    # we will be left with just clumns for each one of the sequences found in the samples
    sp_output_df.drop(columns=['sample_name', 'noName Clade A', 'noName Clade B', 'noName Clade C', 'noName Clade D',
                               'noName Clade E', 'noName Clade F', 'noName Clade G', 'noName Clade H',
                               'noName Clade I', 'raw_contigs', 'post_qc_absolute_seqs', 'post_qc_unique_seqs',
                               'post_taxa_id_absolute_symbiodinium_seqs', 'post_taxa_id_unique_symbiodinium_seqs',
                               'post_taxa_id_absolute_non_symbiodinium_seqs',
                               'post_taxa_id_unique_non_symbiodinium_seqs',
                               'size_screening_violation_absolute', 'size_screening_violation_unique',
                               'post_med_absolute', 'post_med_unique'
                               ], inplace=True)
    sp_output_df = sp_output_df.astype('float')

    colour_palette_div = get_colour_list()
    grey_palette_div = ['#D0CFD4', '#89888D', '#4A4A4C', '#8A8C82', '#D4D5D0', '#53544F']
    # get a list of the sequences in order of their abundance and use this list to create the colour dict
    # the abundances can be got by simply summing up the columns making sure to ommit the last columns
    abundance_dict = {}
    for col in list(sp_output_df):
        abundance_dict[col] = sum(sp_output_df[col])
    # get the names of the sequences sorted according to their totalled abundance
    ordered_list_of_seqs = [x[0] for x in sorted(abundance_dict.items(), key=lambda x: x[1], reverse=True)]
    # create the colour dictionary that will be used for plotting by assigning a colour from the colour_palette
    # to the most abundant seqs first and after that cycle through the grey_pallette assigning colours
    # If we aer only going to have a legend that is cols x rows as shown below, then we should only use
    # that many colours in the plotting.


    colour_dict_div = {}
    for i in range(len(ordered_list_of_seqs)):
        if i < num_leg_cells:
            colour_dict_div[ordered_list_of_seqs[i]] = colour_palette_div[i]
        else:
            grey_index = i % len(grey_palette_div)
            colour_dict_div[ordered_list_of_seqs[i]] = grey_palette_div[grey_index]
    return smp_id_to_smp_name_dict, smp_name_to_smp_id_dict, sp_output_df, colour_dict_div, ordered_list_of_seqs

def process_type_df(path_to_tab_delim_count_type, num_leg_cells):
    sp_output_df_type = pd.read_csv(path_to_tab_delim_count_type, sep='\t', lineterminator='\n',
                                    skiprows=[0, 1, 2, 3, 5],
                                    header=None)
    # get a list of tups that are the seq names and the abundances zipped together
    type_profile_to_abund_tup_list = [(name, int(abund)) for name, abund in
                                      zip(sp_output_df_type.iloc[1][2:].values.tolist(),
                                          sp_output_df_type.iloc[0][2:].values.tolist())]
    # convert the names that are numbers into int strings rather than float strings.
    int_temp_list = []
    for name_abund_tup in type_profile_to_abund_tup_list:
        try:
            int_temp_list.append((str(int(name_abund_tup[0])), int(name_abund_tup[1])))
        except:
            int_temp_list.append((name_abund_tup[0], int(name_abund_tup[1])))
    type_profile_to_abund_tup_list = int_temp_list
    # need to drop the rows that contain the sequence accession and species descriptions
    for i, row_name in enumerate(sp_output_df_type.iloc[:, 0]):
        if 'Sequence accession' in row_name:
            # then we want to drop all rows from here until the end
            index_to_drop_from = i
            break
    sp_output_df_type = sp_output_df_type.iloc[:index_to_drop_from]
    # now drop the sample name columns
    sp_output_df_type.drop(columns=1, inplace=True)
    # make headers
    sp_output_df_type.columns = ['sample_id'] + [a[0] for a in type_profile_to_abund_tup_list]
    # now drop the local abund row and promote the its2_type_prof names to columns headers.
    sp_output_df_type.drop(index=[0, 1], inplace=True)
    sp_output_df_type = sp_output_df_type.set_index(keys='sample_id', drop=True).astype('float')
    # we should plot sample by sample and its2 type by its2 type in the order of the output
    # the problem with doing he convert_to_pastel is that the colours become very similar
    # colour_palette = convert_to_pastel(get_colour_list())
    # Rather, I will attempt to generate a quick set of colours that are pastel and have a minimum distance
    # rule for any colours that are generated from each other.
    # let's do this for 50 colours to start with and see how long it takes.
    # turns out it is very quick. Easily quick enough to do dynamically.
    # When working with pastel colours (i.e. mixing with 255,255,255 it is probably best to work with a smaller dist cutoff
    colour_palette_pas = ['#%02x%02x%02x' % rgb_tup for rgb_tup in
                          create_colour_list(mix_col=(255, 255, 255), sq_dist_cutoff=1000, num_cols=50,
                                             time_out_iterations=10000)]
    # # The below 3d scatter produces a 3d scatter plot to examine the spread of the colours created
    # from mpl_toolkits.mplot3d import Axes3D
    # colour_palette = create_colour_list(sq_dist_cutoff=5000)
    # hex_pal = ['#%02x%02x%02x' % rgb_tup for rgb_tup in colour_palette]
    # colcoords = [list(a) for a in zip(*colour_palette)]
    # print(colcoords)
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # ax.scatter(colcoords[0], colcoords[1], colcoords[2], c=hex_pal, marker='o')
    # colour_palette = get_colour_list()
    grey_palette_type = ['#D0CFD4', '#89888D', '#4A4A4C', '#8A8C82', '#D4D5D0', '#53544F']
    # we will use the col headers as the its2 type profile order for plotting but we
    # we should colour according to the abundance of the its2 type profiles
    # as we don't want to run out of colours by the time we get to profiles that are very abundant.
    # The sorted_type_prof_names_by_local_abund object has the names of the its2 type profile in order of abundance
    # we will use the index order as the order of samples to plot
    # create the colour dictionary that will be used for plotting by assigning a colour from the colour_palette
    # to the most abundant seqs first and after that cycle through the grey_pallette assigning colours
    sorted_type_prof_names_by_local_abund = [a[0] for a in
                                             sorted(type_profile_to_abund_tup_list, key=lambda x: x[1], reverse=True)]


    colour_dict_type = {}
    for i in range(len(sorted_type_prof_names_by_local_abund)):
        if i < num_leg_cells:
            colour_dict_type[sorted_type_prof_names_by_local_abund[i]] = colour_palette_pas[i]
        else:
            grey_index = i % len(grey_palette_type)
            colour_dict_type[sorted_type_prof_names_by_local_abund[i]] = grey_palette_type[grey_index]
    return colour_dict_type, sp_output_df_type, sorted_type_prof_names_by_local_abund

def get_colour_list():
    colour_list = ["#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059", "#FFDBE5",
                  "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87", "#5A0007", "#809693",
                  "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80", "#61615A", "#BA0900", "#6B7900",
                  "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100", "#DDEFFF", "#000035", "#7B4F4B", "#A1C299",
                  "#300018", "#0AA6D8", "#013349", "#00846F", "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744",
                  "#C0B9B2", "#C2FF99", "#001E09", "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68",
                  "#7A87A1", "#788D66", "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED",
                  "#886F4C", "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
                  "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00", "#7900D7",
                  "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700", "#549E79", "#FFF69F",
                  "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329", "#5B4534", "#FDE8DC", "#404E55",
                  "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C", "#83AB58", "#001C1E", "#D1F7CE", "#004B28",
                  "#C8D0F6", "#A3A489", "#806C66", "#222800", "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59",
                  "#8ADBB4", "#1E0200", "#5B4E51", "#C895C5", "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94",
                  "#7ED379", "#012C58", "#7A7BFF", "#D68E01", "#353339", "#78AFA1", "#FEB2C6", "#75797C", "#837393",
                  "#943A4D", "#B5F4FF", "#D2DCD5", "#9556BD", "#6A714A", "#001325", "#02525F", "#0AA3F7", "#E98176",
                  "#DBD5DD", "#5EBCD1", "#3D4F44", "#7E6405", "#02684E", "#962B75", "#8D8546", "#9695C5", "#E773CE",
                  "#D86A78", "#3E89BE", "#CA834E", "#518A87", "#5B113C", "#55813B", "#E704C4", "#00005F", "#A97399",
                  "#4B8160", "#59738A", "#FF5DA7", "#F7C9BF", "#643127", "#513A01", "#6B94AA", "#51A058", "#A45B02",
                  "#1D1702", "#E20027", "#E7AB63", "#4C6001", "#9C6966", "#64547B", "#97979E", "#006A66", "#391406",
                  "#F4D749", "#0045D2", "#006C31", "#DDB6D0", "#7C6571", "#9FB2A4", "#00D891", "#15A08A", "#BC65E9",
                  "#FFFFFE", "#C6DC99", "#203B3C", "#671190", "#6B3A64", "#F5E1FF", "#FFA0F2", "#CCAA35", "#374527",
                  "#8BB400", "#797868", "#C6005A", "#3B000A", "#C86240", "#29607C", "#402334", "#7D5A44", "#CCB87C",
                  "#B88183", "#AA5199", "#B5D6C3", "#A38469", "#9F94F0", "#A74571", "#B894A6", "#71BB8C", "#00B433",
                  "#789EC9", "#6D80BA", "#953F00", "#5EFF03", "#E4FFFC", "#1BE177", "#BCB1E5", "#76912F", "#003109",
                  "#0060CD", "#D20096", "#895563", "#29201D", "#5B3213", "#A76F42", "#89412E", "#1A3A2A", "#494B5A",
                  "#A88C85", "#F4ABAA", "#A3F3AB", "#00C6C8", "#EA8B66", "#958A9F", "#BDC9D2", "#9FA064", "#BE4700",
                  "#658188", "#83A485", "#453C23", "#47675D", "#3A3F00", "#061203", "#DFFB71", "#868E7E", "#98D058",
                  "#6C8F7D", "#D7BFC2", "#3C3E6E", "#D83D66", "#2F5D9B", "#6C5E46", "#D25B88", "#5B656C", "#00B57F",
                  "#545C46", "#866097", "#365D25", "#252F99", "#00CCFF", "#674E60", "#FC009C", "#92896B"]
    return colour_list

def create_colour_list(sq_dist_cutoff=None, mix_col=None, num_cols=50, time_out_iterations=10000, avoid_black_and_white=True):
    new_colours = []
    min_dist = []
    attempt = 0
    while len(new_colours) < num_cols:
        attempt += 1
        # Check to see if we have run out of iteration attempts to find a colour that fits into the colour space
        if attempt > time_out_iterations:
            sys.exit('Colour generation timed out. We have tried {} iterations of colour generation '
                     'and have not been able to find a colour that fits into your defined colour space.\n'
                     'Please lower the number of colours you are trying to find, '
                     'the minimum distance between them, or both.'.format(attempt))
        if mix_col:
            r = int((random.randint(0, 255) + mix_col[0]) /2)
            g = int((random.randint(0, 255) + mix_col[1]) /2)
            b = int((random.randint(0, 255) + mix_col[2]) /2)
        else:
            r = random.randint(0, 255)
            g = random.randint(0, 255)
            b = random.randint(0, 255)

        # now check to see whether the new colour is within a given distance
        # if the avoids are true also
        good_dist = True
        if sq_dist_cutoff:
            dist_list = []
            for i in range(len(new_colours)):
                distance = (new_colours[i][0] - r)**2 + (new_colours[i][1] - g)**2 + (new_colours[i][2] - b)**2
                dist_list.append(distance)
                if distance < sq_dist_cutoff:
                    good_dist = False
                    break
            # now check against black and white
            d_to_black = (r - 0)**2 + (g - 0)**2 + (b - 0)**2
            d_to_white = (r - 255)**2 + (g - 255)**2 + (b - 255)**2
            if avoid_black_and_white:
                if d_to_black < sq_dist_cutoff or d_to_white < sq_dist_cutoff:
                    good_dist = False
            if dist_list:
                min_dist.append(min(dist_list))
        if good_dist:
            new_colours.append((r,g,b))
            attempt = 0

    return new_colours

def remove_axes_but_allow_labels(ax, x_tick_label_list=None):
    ax.set_frame_on(False)
    if not x_tick_label_list:
        ax.set_xticks([])
    ax.set_yticks([])

def plot_data_axes(param_list_one, param_list_two, ax_list, colour_dict_div, colour_dict_type, info_df, ordered_sample_list, smp_id_to_smp_name_dict,
                   smp_name_to_smp_id_dict, sp_output_df_div, sp_output_df_type, param_one, param_two):
    ax_count = 0
    extra_ax_count = 0
    for param_two_cat in param_list_two:
        for param_one_cat in param_list_one:

            ax = ax_list[ax_count]
            patches_list = []
            ind = 0
            colour_list = []

            # for each set of param_one_cat and param_two_cat, we basically want to get a list of the samples
            # that meet the set criteria, we then want to plot samples according to the ordered_sample_list
            # order which will be in IDs. As such we will have to convert the sample_name in the info_df
            # to a sample ID using the smp_name_to_smp_id_dict.

            # get sample_names that fit the requirements
            sample_names_of_set = info_df.loc[
                (info_df[param_one] == param_one_cat) &
                (info_df[param_two] == param_two_cat)
            ].index.values.tolist()


            if not sample_names_of_set:
                remove_axes_but_allow_labels(ax)
                ax_count += 1
                continue
            # convert these to sample IDs
            # The sample names in symportal are actually the full file names version rather than
            # the shorter versions in the info_df. As such we should we will have to do a conversion here
            # TODO in future it should be that the sample names in the excel should match the sampele names in the
            # abundance df. for the time being we will manually make the link.
            smple_ids_of_set = []
            id_to_short_name_dict = {}
            for short_name in sample_names_of_set:
                full_name = '_'.join(info_df.loc[short_name]['fastq_fwd_file_name'].split('/')[-1].split('_')[:3])
                smpl_id = smp_name_to_smp_id_dict[full_name]
                id_to_short_name_dict[smpl_id] = short_name
                smple_ids_of_set.append(smpl_id)

            # full_sample_names = [
            #     '_'.join(info_df.loc[smp_name]['fastq_fwd_file_name'].split('/')[-1].split('_')[:3]) for smp_name in
            #     sample_names_of_set]
            # try:
            #     smple_ids_of_set = [smp_name_to_smp_id_dict[smp_name] for smp_name in full_sample_names]
            # except:
            #     apples = 'asdf'

            # now we want to plot in the order of the ordered_sample_list
            ordered_smple_ids_of_set = [smpl_id for smpl_id in ordered_sample_list if smpl_id in smple_ids_of_set]

            num_smp_in_this_subplot = len(ordered_smple_ids_of_set)
            x_tick_label_list = []

            for smple_id_to_plot in ordered_smple_ids_of_set:


                # General plotting
                sys.stdout.write('\rPlotting sample: {}'.format(smple_id_to_plot))
                x_tick_label_list.append(id_to_short_name_dict[smple_id_to_plot])
                # for each sample we will start at 0 for the y and then add the height of each bar to this

                # PLOT DIVs
                plot_div_over_type(colour_dict_div, colour_list, ind, patches_list, smple_id_to_plot,
                                   sp_output_df_div)

                # PLOT type
                plot_type_under_div(colour_dict_type, colour_list, ind, patches_list, smple_id_to_plot,
                                    sp_output_df_type)
                ind += 1


            paint_rect_to_axes_div_and_type(ax=ax, colour_list=colour_list,
                                            num_smp_in_this_subplot=num_smp_in_this_subplot,
                                            patches_list=patches_list,
                                            x_tick_label_list=x_tick_label_list,
                                            max_num_smpls_in_subplot=10)

            ax_count += 1



def plot_div_over_type(colour_dict_div, colour_list, ind, patches_list, smple_id_to_plot, sp_output_df_div):
    bottom_div = 0
    # for each sequence, create a rect patch
    # the rect will be 1 in width and centered about the ind value.
    for seq in list(sp_output_df_div):
        # class matplotlib.patches.Rectangle(xy, width, height, angle=0.0, **kwargs)
        rel_abund_div = sp_output_df_div.loc[smple_id_to_plot, seq]
        if rel_abund_div > 0:
            patches_list.append(Rectangle((ind - 0.5, bottom_div), 1, rel_abund_div, color=colour_dict_div[seq]))
            # axarr.add_patch(Rectangle((ind-0.5, bottom), 1, rel_abund, color=colour_dict[seq]))
            colour_list.append(colour_dict_div[seq])
            bottom_div += rel_abund_div

def plot_type_under_div(colour_dict_type, colour_list, ind, patches_list, smple_id_to_plot, sp_output_df_type):
    # the idea of the type is to put it as a reflection below the y=0 line
    # as such we should just want to make everything negative
    bottom_type = 0
    # for each sequence, create a rect patch
    # the rect will be 1 in width and centered about the ind value.
    # we want to plot the rects so that they add to 1. As such we want to divide
    # each value by the total for that sample.
    tot_for_sample = sp_output_df_type.loc[smple_id_to_plot].sum()
    for its2_profile in list(sp_output_df_type):
        rel_abund = sp_output_df_type.loc[smple_id_to_plot, its2_profile]
        if rel_abund > 0:
            depth = -0.2 * (rel_abund / tot_for_sample)
            patches_list.append(
                Rectangle((ind - 0.5, bottom_type), 1, depth,
                          color=colour_dict_type[its2_profile]))
            # axarr.add_patch(Rectangle((ind-0.5, bottom), 1, rel_abund, color=colour_dict[seq]))
            colour_list.append(colour_dict_type[its2_profile])
            bottom_type += depth

def paint_rect_to_axes_div_and_type(ax, colour_list, num_smp_in_this_subplot,  patches_list, coral_csw_x_val_tup_list=None, x_tick_label_list=None,  max_num_smpls_in_subplot=10):
    # We can try making a custom colour map
    # https://matplotlib.org/api/_as_gen/matplotlib.colors.ListedColormap.html
    this_cmap = ListedColormap(colour_list)
    # here we should have a list of Rectangle patches
    # now create the PatchCollection object from the patches_list
    patches_collection = PatchCollection(patches_list, cmap=this_cmap)
    patches_collection.set_array(np.arange(len(patches_list)))
    # if n_subplots is only 1 then we can refer directly to the axarr object
    # else we will need ot reference the correct set of axes with i
    # Add the pathces to the axes
    ax.add_collection(patches_collection)
    ax.autoscale_view()
    ax.figure.canvas.draw()
    # also format the axes.
    # make it so that the x axes is constant length
    ax.set_xlim(0 - 0.5, max_num_smpls_in_subplot - 0.5)
    ax.set_ylim(-0.2, 1)
    ax.set_xticks(range(num_smp_in_this_subplot))
    ax.set_xticklabels(x_tick_label_list, rotation='vertical', fontsize=6)

    remove_axes_but_allow_labels(ax, x_tick_label_list)

    # as well as getting rid of the top and right axis splines
    # I'd also like to restrict the bottom spine to where there are samples plotted but also
    # maintain the width of the samples
    # I think the easiest way to do this is to hack a bit by setting the x axis spines to invisible
    # and then drawing on a line at y = 0 between the smallest and largest ind (+- 0.5)
    # ax.spines['bottom'].set_visible(False)
    ax.add_line(Line2D((0 - 0.5, num_smp_in_this_subplot - 0.5), (0, 0), linewidth=2, color='black'))

    # once we have added the black line we should add the grey and brown line that will associate the
    # coral samples to the csw samples
    if coral_csw_x_val_tup_list:
        for coral_csw_x_value_tup in coral_csw_x_val_tup_list:
            ax.add_line(Line2D((coral_csw_x_value_tup[0], coral_csw_x_value_tup[1]), (0, 0), linewidth=2, color=coral_csw_x_value_tup[2]))

d_s_p = '/Users/humebc/Google_Drive/projects/camp_gardener/SymPortal_submission_input_Emma_Camp_NewCal.xlsx'




create_parameter_matrix_fig(data_sheet_path = d_s_p, param_one = 'host_species', param_two = 'site_name')



