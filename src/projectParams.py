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
import math

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

            # populate a background values array with non-zero cells (and filter values to within study area)
            #   tuples of (x, y, easting, northing, value)
            self.background_values = []
            cellcentre = round(self.background_csize / 2)
            bg_histogram = {}
            out_of_range = 0

            for (X,Y) in np.ndindex(self.background_array.shape):
                val = self.background_array[X,Y]
                if val not in bg_histogram.keys():
                    bg_histogram[val] = 1
                else:
                    bg_histogram[val] += 1
                if val > 0.001:  # threshold for inclusion
                    bg_E = self.background_bl_east + self.background_csize * X + cellcentre
                    bg_N = self.background_bl_north + self.background_csize * Y + cellcentre
                    if bg_E >= self.sarea_bl_east and bg_E <= self.sarea_tr_east \
                            and bg_N >= self.sarea_bl_north and bg_N <= self.sarea_tr_north:
                        self.background_values.append((X, Y, bg_E, bg_N, val))
                    else:
                        out_of_range += 1

            logging.info('  Non-zero background cells within Study Area: ' + str(len(self.background_values)) + '  Out of range: ' + str(out_of_range))

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

                header = pd.read_csv(filename, index_col=False, header=None, nrows=23)

                # identify locations of various data columns and their defaults
                col_ID = int(header.iloc[4, 1])
                col_E = int(header.iloc[5, 1])
                col_N = int(header.iloc[6, 1])
                col_Pop = int(header.iloc[7, 1])
                col_PopSubGroups = int(header.iloc[8, 1])

                lst = header.iloc[9, 1:col_PopSubGroups+1].to_list()
                col_PopSubGroupsDefaults = [int(elem) for elem in lst]

                lst = header.iloc[9, col_PopSubGroups + 1:2*col_PopSubGroups+1].to_list()
                col_PopSubGroupsCols = [int(elem) for elem in lst]

                lst = header.iloc[16, 1:col_PopSubGroups+1].to_list()
                col_MobSubGroupsDefaults = [int(elem) for elem in lst]

                lst = header.iloc[16, col_PopSubGroups + 1:2*col_PopSubGroups+1].to_list()
                col_MobSubGroupsCols = [int(elem) for elem in lst]

                col_LD = int(header.iloc[12, 2])
                col_LDdefault = int(header.iloc[12, 1])

                # create a Pandas dataframe and miss out the first 23 rows
                # header information is on row 23, so start here
                origin_df = pd.read_csv(filename, header=[22])

                # filter the data to within the Study area
                # TODO: use column numbers instead of names here
                bbquery = 'OSEAST >= {} & OSEAST <= {} & OSNRTH >= {} & OSNRTH <= {}'.format(
                    self.sarea_bl_east, self.sarea_tr_east, self.sarea_bl_north, self.sarea_tr_north)
                csvData = origin_df.query(bbquery)

                # all origin data to be stored inside origins dictionary
                self.origin_data = {}

                self.origin_data['ID'] = csvData.iloc[:,col_ID-1].to_list()

                self.origin_data['eastings'] = [round(elem) for elem in csvData.iloc[:,col_E-1].to_list()]
                self.origin_data['northings'] = [round(elem) for elem in csvData.iloc[:, col_N - 1].to_list()]

                # get the total population data
                pop_data = csvData.iloc[:, col_Pop-1]
                self.origin_data['pop_data'] = pop_data.to_list()
                self.origin_data['pop_name'] = pop_data.name

                pop_subgroups = csvData.iloc[:, min(col_PopSubGroupsCols)-1:max(col_PopSubGroupsCols)]

                # read in the population subgroups data, calculate the population (for each age band)
                self.origin_data['subgroup_names'] = []
                self.origin_data['subgroups_data'] = {}  # we may not need to store these eventually
                self.origin_data['subgroups_pop'] = {}  # matching dictionary for actual populations
                subgroups_total = 0

                idx = 0  # which population data subgroup this is
                for (columnName, columnData) in pop_subgroups.iteritems():  # loop through the columns
                    self.origin_data['subgroup_names'].append(columnName)
                    self.origin_data['subgroups_data'][columnName] = []
                    for pop in columnData.to_list():
                        if not pop or math.isnan(pop):
                            self.origin_data['subgroups_data'][columnName].append(col_PopSubGroupsDefaults[idx] / 100)
                        else:
                            self.origin_data['subgroups_data'][columnName].append(pop)

                    idx += 1

                    self.origin_data['subgroups_pop'][columnName] = [(a * b / 100.0) for a, b in zip(self.origin_data['subgroups_data'][columnName],pop_data)]
                    subgroups_total = subgroups_total + sum(self.origin_data['subgroups_pop'][columnName])

                mob_subgroups = csvData.iloc[:, min(col_MobSubGroupsCols) - 1:max(col_MobSubGroupsCols)]

                # read in the mobility data (by age band)
                self.origin_data['subgroups_mob'] = {}

                idx = 0   # which mob subgroup this is
                for (columnName, columnData) in mob_subgroups.iteritems():  # loop through the columns
                    self.origin_data['subgroups_mob'][columnName] = []
                    for mob in columnData.to_list():  # loop through mobility data groups, apply default value
                        if not mob or math.isnan(mob):
                            self.origin_data['subgroups_mob'][columnName].append(col_MobSubGroupsDefaults[idx])
                        else:
                            self.origin_data['subgroups_mob'][columnName].append(mob)
                    idx += 1

                # read in the local dispersion data

                self.origin_data['LD'] = []
                for val in csvData.iloc[:, col_LD - 1].to_list():
                    if not val or math.isnan(val):
                        self.origin_data['LD'].append(col_LDdefault)
                    else:
                        self.origin_data['LD'].append(val)

                logging.info('  Origin Eastings  (.projParams.origin_data[eastings]):  ' + str(
                    self.origin_data['eastings'][0:5]) + ' ... Count: ' + str(len(self.origin_data['eastings'])))
                logging.info('  Origin Northings (.projParams.origin_data[northings]): ' + str(
                    self.origin_data['northings'][0:5]) + ' ... Count: ' + str(len(self.origin_data['northings'])))
                # logging.info('  Out of Study Area range origins: ' + str(out_of_range))

                self.origin_eastings_min = min(self.origin_data['eastings'])
                self.origin_eastings_max = max(self.origin_data['eastings'])
                self.origin_northings_min = min(self.origin_data['northings'])
                self.origin_northings_max = max(self.origin_data['northings'])

                logging.info('  BL Origin point  (.projParams.origin_eastings_min / northings_min): ' + str(
                    self.origin_eastings_min) + ', ' + str(self.origin_northings_min))
                logging.info('  TR Origin point  (.projParams.origin_eastings_max / northings_max): ' + str(
                    self.origin_eastings_max) + ', ' + str(self.origin_northings_max))

                # calculate Cell X and Y references for each origin
                csize = self.background_csize
                self.origin_data['XY'] = []  # array of grid X,Y coords, corresponding with the Easting/Northing values

                for origin in range(0, len(self.origin_data['eastings'])):
                    orig_E = self.origin_data['eastings'][origin]
                    orig_N = self.origin_data['northings'][origin]
                    orig_X = round((orig_E - self.background_bl_east) / csize)
                    orig_Y = round((orig_N - self.background_bl_north) / csize)
                    self.origin_data['XY'].append((orig_X, orig_Y))

                logging.info('  Origin X,Y indexes (.projParams.origin_data[XY]): '
                             + str(self.origin_data['XY'][0:5]) + ' ... Count: ' + str(len(self.origin_data['XY'])))

                logging.info('  Origin Population name (.projParams.origin_data[pop_name]): ' + self.origin_data['pop_name'])
                logging.info('  Origin Population data (.projParams.origin_data[pop_data]): '
                             + str(self.origin_data['pop_data'][0:5]) + ' ... Count: ' + str(len(self.origin_data['pop_data']))
                             + ' Total: ' + str(sum(self.origin_data['pop_data'])))

                logging.info('  Origin Population subgroups (.projParams.origin_data[subgroup_names]): '
                             + str(self.origin_data['subgroup_names']))
                
                logging.info('  Origin Population first subgroup data (.projParams.origin_data[subgroups_data]): '
                             + str(self.origin_data['subgroups_data'][self.origin_data['subgroup_names'][0]][0:5])
                             + ' ... Count: ' + str(len(self.origin_data['subgroups_data'][self.origin_data['subgroup_names'][0]])))

                logging.info('  Origin Population first subgroup pops (.projParams.origin_data[subgroups_pop]): '
                             + str(self.origin_data['subgroups_pop'][self.origin_data['subgroup_names'][0]][0:5])
                             + ' ... Count: ' + str(len(self.origin_data['subgroups_pop'][self.origin_data['subgroup_names'][0]]))
                             + ' Total: ' + str(sum(self.origin_data['subgroups_pop'][self.origin_data['subgroup_names'][0]])))

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

                    header = pd.read_csv(filename, index_col=False, header=None, nrows=23)

                    # identify locations of various data columns and their defaults
                    col_ID = int(header.iloc[4, 1])
                    col_E = int(header.iloc[5, 1])
                    col_N = int(header.iloc[6, 1])
                    col_Pop = int(header.iloc[7, 1])
                    col_PopSubGroups = int(header.iloc[8, 1])

                    lst = header.iloc[9, 1:col_PopSubGroups + 1].to_list()
                    col_PopSubGroupsDefaults = [int(elem) for elem in lst]

                    lst = header.iloc[9, col_PopSubGroups + 1:2 * col_PopSubGroups + 1].to_list()
                    col_PopSubGroupsCols = [int(elem) for elem in lst]

                    col_TimeProfile = int(header.iloc[10, 2])
                    col_TimeProfileDefault = header.iloc[10, 1].upper()  # UPPER to always match timeseries

                    col_LD = int(header.iloc[12, 2])
                    col_LDdefault = int(header.iloc[12, 1])

                    col_WAD = int(header.iloc[13, 2])
                    col_WADdefault = header.iloc[13, 1]

                    lst = header.iloc[14, 1:col_PopSubGroups + 1].to_list()
                    col_MajorFlowCols = [int(elem) for elem in lst]

                    # create a Pandas dataframe and miss out the first 23 rows
                    dest_df = pd.read_csv(filename, header=[22])  # , index_col=False, header=None, skiprows=23)

                    # filter the data to within the Study area
                    # TODO: use column numbers instead of names here
                    bbquery = 'OSEAST >= {} & OSEAST <= {} & OSNRTH >= {} & OSNRTH <= {}'.format(
                        self.sarea_bl_east, self.sarea_tr_east, self.sarea_bl_north, self.sarea_tr_north)
                    csvData = dest_df.query(bbquery)

                    logging.info('    Rows within Study Area: ' + str(len(csvData)))
                    total_rows = total_rows + len(csvData)

                    dest_data['ID'] = csvData.iloc[:, col_ID - 1].to_list()

                    dest_data['eastings'] = [round(elem) for elem in csvData.iloc[:, col_E - 1].to_list()]
                    dest_data['northings'] = [round(elem) for elem in csvData.iloc[:, col_N - 1].to_list()]

                    logging.info('    IDs: ' + str(dest_data['ID'][0:5]))
                    logging.info('    Eastings: ' + str(dest_data['eastings'][0:5]))
                    logging.info('    Northings: ' + str(dest_data['northings'][0:5]))

                    # calculate Cell X and Y references for each destination
                    csize = self.background_csize
                    dest_data['XY'] = []  # array of grid X,Y coords, corresponding with the Easting/Northing values

                    for dest in range(0, len(dest_data['eastings'])):
                        dest_E = dest_data['eastings'][dest]
                        dest_N = dest_data['northings'][dest]
                        dest_X = round((dest_E - self.background_bl_east) / csize)
                        dest_Y = round((dest_N - self.background_bl_north) / csize)
                        dest_data['XY'].append((dest_X, dest_Y))

                    logging.info('    X,Y indexes: '
                                 + str(dest_data['XY'][0:5]) + ' ... Count: ' + str(len(dest_data['XY'])))

                    # get the total population data
                    pop_data = csvData.iloc[:, col_Pop - 1]
                    dest_data['pop_data'] = pop_data.to_list()
                    # dest_data['pop_data'] = list(map(round, dest_data['pop_data_orig']))
                    dest_data['pop_name'] = pop_data.name

                    logging.info('    Pop name: ' + dest_data['pop_name'])
                    logging.info('    Pop data: ' + str(dest_data['pop_data'][0:5])
                                 + '... Total: ' + str(sum(dest_data['pop_data'])))
                    #logging.info('    Pop data (rounded): ' + str(dest_data['pop_data'][0:5])
                    #             + ' Total: ' + str(sum(dest_data['pop_data'])))

                    pop_subgroups = csvData.iloc[:, min(col_PopSubGroupsCols) - 1:max(col_PopSubGroupsCols)]

                    # get the age breakdown data, add names to a list and data arrays to a dictionary
                    #   we may well want to use a different (more efficient) structure here

                    dest_data['subgroup_names'] = []
                    dest_data['subgroups_data'] = {}  # we may not need to store these eventually
                    dest_data['subgroups_pop'] = {}  # matching dictionary for actual populations
                    subgroups_total = 0

                    idx = 0  # which population data subgroup this is
                    for (columnName, columnData) in pop_subgroups.iteritems():  # loop through the columns
                        dest_data['subgroup_names'].append(columnName)
                        dest_data['subgroups_data'][columnName] = []
                        for pop in columnData.to_list():
                            if not pop or math.isnan(pop):
                                dest_data['subgroups_data'][columnName].append(
                                    col_PopSubGroupsDefaults[idx] / 100)
                            else:
                                dest_data['subgroups_data'][columnName].append(pop)

                        idx += 1

                        dest_data['subgroups_pop'][columnName] = [(a * b / 100.0) for a, b in zip(
                            dest_data['subgroups_data'][columnName], pop_data)]
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


                    # read in the time profile data

                    dest_data['time_profiles'] = []
                    for val in csvData.iloc[:, col_TimeProfile - 1].to_list():
                        if not val or str(val) == 'nan':
                            dest_data['time_profiles'].append(col_TimeProfileDefault)
                        else:
                            dest_data['time_profiles'].append(val)

                    # read in the local dispersion data

                    dest_data['LD'] = []
                    for val in csvData.iloc[:, col_LD - 1].to_list():
                        if not val or math.isnan(val):
                            dest_data['LD'].append(col_LDdefault)
                        else:
                            dest_data['LD'].append(val)

                    # read in the WAD data
                    wads = []
                    wad_default = col_WADdefault.split('|')
                    for val in csvData.iloc[:, col_WAD - 1].to_list():
                        if not val:
                            wads.append(wad_default)
                        else:
                            wads.append(val.split('|'))

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
        # amount and list are appended to hold totals and indexes of origins, background cells, etc. later
        # [  [ [radius, percent, amount, list], ... ], [ ... ] ... ]

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

                wad_pair = [ radius, perc, 0, [] ]  # NOT a tuple, use array so we can append to it later
                wad_pairs_list.append(wad_pair)

            dest_wad.append(wad_pairs_list)

        return dest_wad
