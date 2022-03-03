#
# Python version of Surface Builder 24/7
#
# Jan 2022
# GeoData Institute
# University of Southampton
# on behalf of ONS

# Import core modules

import logging
import datetime

# Import additional modules (which will need to be installed)

import numpy as np
import pandas as pd


# A class for reading in and storing all Project parameters

class ProjectParams:

    def loadFromDict(self, dictionary):

        self.analysisarray = dictionary['analysisarray']
        self.buffer = dictionary['buffer']
        self.background = dictionary['background']
        self.timeseries = dictionary['timeseries']
        self.origin = dictionary['origin']
        self.destarray = dictionary['destarray']

    def loadFromFile(self, filename):

        try:
            with open(filename, 'r') as file_opened:
                text = file_opened.read()
                projparamslist = text.split()

                logging.debug('Full projectParamsList: ' + str(projparamslist))

                # set analysisparams to first list line, split into array by commas
                analysisParams = projparamslist[0].split(',')

                # turn this into an integer array assigned to analysisarray
                self.analysisarray = list(map(int, analysisParams))

                # assign buffer to second line (make integer)
                self.buffer = int(projparamslist[1])

                # assign background to third line
                self.background = projparamslist[2]

                # assign timeseries to fourth line
                self.timeseries = projparamslist[3]

                # assign origin to fifth line
                self.origin = projparamslist[4]

                # assign destarray to fifth line
                self.destarray = projparamslist[5].split(',')

        except IOError as e:
            logging.error(e)

    def outputProjectParams(self):

        logging.info('Loaded Project Params:')
        logging.info('  Analysis Array    (.projParams.analysisarray): ' + str(self.analysisarray))
        logging.info('  Buffer            (.projParams.buffer):        ' + str(self.buffer))
        logging.info('  Background        (.projParams.background):    ' + self.background)
        logging.info('  Timeseries        (.projParams.timeseries):    ' + self.timeseries)
        logging.info('  Origin            (.projParams.origin):        ' + self.origin)
        logging.info('  Destination Array (.projParams.destarray):     ' + str(self.destarray))

    def loadBackground(self, filename):

        header_rows = 6  # six rows for header information

        # store header information as a dictionary
        # including: ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value
        self.background_header = {}
        row_iterate = 1

        try:
            with open(filename, 'r') as file_opened:
                for line in file_opened:
                    if row_iterate <= header_rows:
                        header_values = line.split()  # split line using whitespace as separator
                        header_var = header_values[0]
                        header_val = header_values[1]

                        # add var:value to the header dictionary
                        #   n.b. some vals are float so needs to be converted
                        self.background_header[header_var] = int(float(header_val))
                        row_iterate = row_iterate + 1
                    else:
                        break

            # read background data array - up down flip applied due to ASCII grid top to bottom writing
            self.background_array = np.flipud(np.loadtxt(filename, skiprows=header_rows, dtype='float64'))

            # set some values from the header 

            self.background_bl_east = int(self.background_header['xllcorner'])
            self.background_bl_north = int(self.background_header['yllcorner'])
            self.background_rows = int(self.background_header['nrows'])
            self.background_cols = int(self.background_header['ncols'])
            self.background_csize = int(self.background_header['cellsize'])

            # background top right easting = bottom left easting + (cols*cellsize)
            self.background_tr_east = self.background_bl_east \
                                      + (self.background_cols * self.background_csize)

            # background top right northing = bottom left northing + (rows*cellsize)
            self.background_tr_north = self.background_bl_north \
                                       + (self.background_rows * self.background_csize)

            logging.info('Background Header:')
            logging.info('  BL Easting  (.projParams.background_bl_east):  ' + str(self.background_bl_east))
            logging.info('  BL Northing (.projParams.background_bl_north): ' + str(self.background_bl_north))
            logging.info('  TR Easting  (.projParams.background_tr_east):  ' + str(self.background_tr_east))
            logging.info('  TR Northing (.projParams.background_tr_north): ' + str(self.background_tr_north))
            logging.info('  Rows        (.projParams.background_rows):     ' + str(self.background_rows))
            logging.info('  Cols        (.projParams.background_cols):     ' + str(self.background_cols))
            logging.info('  Cellsize    (.projParams.background_csize):    ' + str(self.background_csize))

        except IOError as e:
            logging.error(e)

    def loadTimeSeries(self, filename):
        try:

            with open(filename, 'r') as file_opened:
                ts_data = pd.read_excel(filename)

        except IOError as e:
            logging.error(e)

        # create an empty dictionary for timeseries data, structure:
        #   { 'TS01':
        #         { 'InTravel': [ (time1, weight1), .... ],
        #           'OnSite':   [ (time1, weight1), .... ]  }
        # etc - dictionary -> dictionary -> array

        self.timeseries_data = {}
        msg = ''

        for (columnName, columnData) in ts_data.iteritems():  # loop through the columns

            if columnName[:7] != 'Unnamed':  # a new timeseries key
                ts_times = columnData.values
                ts_key = columnName.upper()  # store for the next set of values (weights), UPPER case
                self.timeseries_data[ts_key] = {}  # create a new dictionary for this key
                msg = msg + ts_key + '; '

            else:
                ts_weights = columnData.values

                for row in range(0, len(ts_times)):  # now go through every row

                    if ts_times[row] == 'InTravel':  # the start of an InTravel category
                        category = 'InTravel'
                        # self.timeseries_data[ts_key][category] = []  # empty array for the values
                        self.timeseries_data[ts_key][category] = {}  # empty dict for the values

                    elif ts_times[row] == 'OnSite':  # the start of an OnSite category
                        category = 'OnSite'
                        # self.timeseries_data[ts_key][category] = []  # empty array for the values
                        self.timeseries_data[ts_key][category] = {}  # empty dict for the values

                    else:
                        # if a valid set of values add it to the current ts dictionary array
                        if str(ts_times[row]) != 'nan' and str(ts_weights[row]) != 'nan':
                            logging.debug('TS: ' + category + ' ' + ts_key + ' '
                                          + str(ts_times[row]) + ' ' + str(ts_weights[row]))
                            # self.timeseries_data[ts_key][category].append((ts_times[row],ts_weights[row]))
                            timekey = datetime.time(ts_times[row].hour, ts_times[row].minute)
                            self.timeseries_data[ts_key][category][timekey] = ts_weights[row]

        logging.info('  TS Keys: ' + msg)
        logging.info('  TS Last key/category/time/weight: ' + ts_key + ' ' + category + ' '
                     + str(timekey) + ' ' + str(self.timeseries_data[ts_key][category][timekey]))

    def loadOrigin(self, filename):
        try:

            with open(filename, 'r') as file_opened:

                # TODO - read and process header info on first 22 lines
                # header information is on row 23, so start here
                originData = pd.read_csv(filename, header=[22])

                # get the output area codes from the first column
                self.origin_OA = originData.iloc[:, 0].to_list()

                logging.info('  Origin Output Areas  (.projParams.origin_OA):  ' + str(
                    self.origin_OA[0:5]) + ' ... Count: ' + str(len(self.origin_OA)))

                orig_east_list = originData.OSEAST.to_list()
                # round the coords as they are float
                self.origin_eastings = [round(elem) for elem in orig_east_list]

                orig_nrth_list = originData.OSNRTH.to_list()
                # round the coords as they are float
                self.origin_northings = [round(elem) for elem in orig_nrth_list]

                logging.info('  Origin Eastings  (.projParams.origin_eastings):  ' + str(
                    self.origin_eastings[0:5]) + ' ... Count: ' + str(len(self.origin_eastings)))
                logging.info('  Origin Northings (.projParams.origin_northings): ' + str(
                    self.origin_northings[0:5]) + ' ... Count: ' + str(len(self.origin_northings)))

                self.origin_eastings_min = min(self.origin_eastings)
                self.origin_eastings_max = max(self.origin_eastings)
                self.origin_northings_min = min(self.origin_northings)
                self.origin_northings_max = max(self.origin_northings)

                logging.info('  BL Origin point  (.projParams.origin_eastings_min / northings_min): ' + str(
                    self.origin_eastings_min) + ', ' + str(self.origin_northings_min))
                logging.info('  TR Origin point  (.projParams.origin_eastings_max / northings_max): ' + str(
                    self.origin_eastings_max) + ', ' + str(self.origin_northings_max))

                # calculate Cell X and Y references for each origin
                csize = self.background_csize
                self.origin_XY = []  # array of grid X,Y coords, corresponding with the Easting/Northing values

                for origin in range(0, len(self.origin_eastings)):
                    orig_E = self.origin_eastings[origin]
                    orig_N = self.origin_northings[origin]
                    orig_X = round((orig_E - self.background_bl_east) / csize)
                    orig_Y = round((orig_N - self.background_bl_north) / csize)
                    self.origin_XY.append((orig_X, orig_Y))

                logging.info('  Origin X,Y indexes (.projParams.origin_XY): '
                             + str(self.origin_XY[0:5]) + ' ... Count: ' + str(len(self.origin_XY)))

                # get the total population data
                pop_data = originData.iloc[:, 4]
                self.origin_pop_name = pop_data.name
                self.origin_pop_data = pop_data.to_list()

                logging.info('  Origin Population name (.projParams.origin_pop_name): ' + self.origin_pop_name)
                logging.info('  Origin Population data (.projParams.origin_pop_data): '
                             + str(self.origin_pop_data[0:5]) + ' ... Count: ' + str(len(self.origin_pop_data))
                             + ' Total: ' + str(sum(self.origin_pop_data)))

                # get the age breakdown data, add names to a list and data arrays to a dictionary
                #   we may well want to use a different (more efficient) structure here

                self.origin_subgroup_names = []
                self.origin_subgroups_data = {}  # we may not need to store these eventually
                self.origin_subgroups_pop = {}  # matching dictionary for actual populations
                subgroups_total = 0

                for (columnName, columnData) in originData.iloc[:, 5:].iteritems():  # loop through the columns
                    logging.debug('Origin Group Data: ' + columnName + ' '
                                  + str(len(columnData)) + '\n' + str(columnData[0:4]))

                    # if there's a number in the first row we'll consider it as valid
                    if str(columnData.values[0]) != 'nan':
                        self.origin_subgroup_names.append(columnName)
                        self.origin_subgroups_data[columnName] = columnData.to_list()

                        # calculate subgroup population (total pop * percent / 100 rounded to int)
                        self.origin_subgroups_pop[columnName] = [(a * b / 100.0) for a, b in
                                                                 zip(self.origin_subgroups_data[columnName],
                                                                     self.origin_pop_data)]
                        subgroups_total = subgroups_total + sum(self.origin_subgroups_pop[columnName])

                logging.info('  Origin Population subgroups (.projParams.origin_subgroup_names): '
                             + str(self.origin_subgroup_names))
                logging.info('  Origin Population first subgroup data (.projParams.origin_subgroups_data): '
                             + str(self.origin_subgroups_data[self.origin_subgroup_names[0]][0:5])
                             + ' ... Count: ' + str(len(self.origin_subgroups_data[self.origin_subgroup_names[0]])))

                logging.info('  Origin Population first subgroup pops (.projParams.origin_subgroups_pop): '
                             + str(self.origin_subgroups_pop[self.origin_subgroup_names[0]][0:5])
                             + ' ... Count: ' + str(len(self.origin_subgroups_pop[self.origin_subgroup_names[0]]))
                             + ' Total: ' + str(sum(self.origin_subgroups_pop[self.origin_subgroup_names[0]])))

                logging.info('  Origin Population all subgroups pops total: ' + str(round(subgroups_total, 2)))

        except IOError as e:
            logging.error(e)

    def loadDestFiles(self, pathname):

        self.destination_data = []  # array of each destination dictionary
        total_rows = 0

        for dest_file in self.destarray:

            filename = pathname + dest_file
            logging.info('\n  Reading: ' + filename + '\n')

            try:
                with open(filename, 'r') as file_opened:  # currently unused, but checks existence!

                    dest_data = {}  # empty dictionary to hold destination data

                    dest_header = pd.read_csv(filename, index_col=False, header=None, nrows=23)

                    dest_data['time_profile'] = dest_header.iloc[10, 1].upper()  # UPPER to always match timeseries

                    logging.info('    Time profile: ' + dest_data['time_profile'])

                    # create a Pandas dataframe and miss out the first 23 rows

                    dest_df = pd.read_csv(filename, header=[22])  # , index_col=False, header=None, skiprows=23)

                    # just the eastings, northings and WAD for now.
                    # dest_read = dest_df.iloc[:, [2,3,14]]

                    logging.info('    Rows: ' + str(len(dest_df)))
                    total_rows = total_rows + len(dest_df)

                    dest_data['OA'] = dest_df.iloc[:, 0].to_list()
                    dest_data['eastings'] = dest_df.iloc[:, 2].to_list()
                    dest_data['northings'] = dest_df.iloc[:, 3].to_list()

                    logging.info('    OAs: ' + str(dest_data['OA'][0:5]))
                    logging.info('    Eastings: ' + str(dest_data['eastings'][0:5]))
                    logging.info('    Northings: ' + str(dest_data['northings'][0:5]))

                    # calculate Cell X and Y references for each destination
                    csize = self.background_csize
                    dest_data[
                        'XY'] = []  # array of grid X,Y coords, corresponding with the Easting/Northing values

                    for dest in range(0, len(dest_data['eastings'])):
                        dest_E = dest_data['eastings'][dest]
                        dest_N = dest_data['northings'][dest]
                        dest_X = round((dest_E - self.background_bl_east) / csize)
                        dest_Y = round((dest_N - self.background_bl_north) / csize)
                        dest_data['XY'].append((dest_X, dest_Y))

                    logging.info('    X,Y indexes: '
                                 + str(dest_data['XY'][0:5]) + ' ... Count: ' + str(len(dest_data['XY'])))

                    pop_column = int(dest_header.iloc[7, 1]) - 1

                    pop_data = dest_df.iloc[:, pop_column]
                    dest_data['pop_name'] = pop_data.name
                    dest_data['pop_data_orig'] = pop_data.to_list()
                    dest_data['pop_data'] = list(map(round, dest_data['pop_data_orig']))

                    logging.info('    Pop name: ' + dest_data['pop_name'])
                    logging.info('    Pop data orig: ' + str(dest_data['pop_data_orig'][0:5])
                                 + ' Total: ' + str(sum(dest_data['pop_data'])))
                    logging.info('    Pop data (rounded): ' + str(dest_data['pop_data'][0:5])
                                 + ' Total: ' + str(sum(dest_data['pop_data'])))

                    pop_subgroups_count = int(dest_header.iloc[8, 1])

                    sg_start = pop_column + 1
                    sg_end = sg_start + pop_subgroups_count
                    pop_subgroups = dest_df.iloc[:, sg_start:sg_end]

                    # get the age breakdown data, add names to a list and data arrays to a dictionary
                    #   we may well want to use a different (more efficient) structure here

                    dest_data['subgroup_names'] = []
                    dest_data['subgroups_data'] = {}  # we may not need to store these eventually
                    dest_data['subgroups_pop'] = {}  # matching dictionary for actual populations
                    subgroups_total = 0

                    for (columnName, columnData) in pop_subgroups.iteritems():  # loop through the columns
                        logging.debug('Destination Group Data: ' + columnName + ' '
                                      + str(len(columnData)) + '\n' + str(columnData[0:4]))

                        # if there's a number in the first row we'll consider it as valid
                        if str(columnData.values[0]) != 'nan':
                            dest_data['subgroup_names'].append(columnName)
                            dest_data['subgroups_data'][columnName] = columnData.to_list()

                            # calculate subgroup population (total pop * percent / 100 rounded to int)
                            dest_data['subgroups_pop'][columnName] = [(a * b / 100.0) for a, b in
                                                                      zip(dest_data['subgroups_data'][columnName],
                                                                          dest_data['pop_data'])]
                            subgroups_total = subgroups_total + sum(dest_data['subgroups_pop'][columnName])

                    logging.info('    Population subgroups: ' + str(dest_data['subgroup_names']))
                    logging.info('    Population first subgroup data: '
                                 + str(dest_data['subgroups_data'][dest_data['subgroup_names'][0]][0:5])
                                 + ' ... Count: ' + str(
                        len(dest_data['subgroups_data'][dest_data['subgroup_names'][0]])))
                    logging.info('    Population first subgroup pop: '
                                 + str(dest_data['subgroups_pop'][dest_data['subgroup_names'][0]][0:5])
                                 + ' ... Count: ' + str(len(dest_data['subgroups_pop'][dest_data['subgroup_names'][0]]))
                                 + ' Total: ' + str(sum(dest_data['subgroups_pop'][dest_data['subgroup_names'][0]])))
                    logging.info('    Dest file Population all subgroups pops total: ' + str(round(subgroups_total, 2)))

                    wads = dest_df.iloc[:, 14].str.split('|')
                    dest_data['WAD'] = self.covertWADs(wads)

                    logging.info('    WADs: ' + str(dest_data['WAD'][0:5]))
                    logging.info('      First tuple/radius/percent: '
                                 + str(dest_data['WAD'][0][0]) + ' '
                                 + str(dest_data['WAD'][0][0][0]) + ' '
                                 + str(dest_data['WAD'][0][0][1]))

                    # append to the destination array
                    self.destination_data.append(dest_data)

            except IOError as e:
                logging.error(e)

        logging.info('\n  Total Destination datasets / rows (.projParams.destination_data): '
                     + str(len(self.destination_data)) + ' / ' + str(total_rows))

        # go across to numpy array at this point? 
        # dest_array = np.column_stack([dest_east, dest_north, dest_wad])

    def covertWADs(self, wad_list):
        # take a list of WADs as strings, return a new list of lists of tuples (now lists, so we can add to them):
        # [  [ [radius, percent], ... ], [ ... ] ... ]

        # get rid of the pipe separator
        # dest_df_wad = wad_list.str.split('|')

        # Empty list of all WADs for ALL destinations
        dest_wad = []

        for each_row in wad_list:  # each row of data

            wad_pairs_list = []
            for wad_pairs in each_row:  # each pair of radii/perc (or the final perc)

                if '>' in wad_pairs:
                    pair = list(wad_pairs.split('>'))
                    radius = int(pair[0])
                    perc = int(pair[1])
                else:
                    perc = int(wad_pairs)
                    radius = 0  # set the radius to zero if not specified

                wad_pair = [radius, perc]  # NOT a tuple, use array so we can append to it later
                wad_pairs_list.append(wad_pair)

            dest_wad.append(wad_pairs_list)

        return dest_wad
