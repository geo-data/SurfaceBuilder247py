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

from locationIndex import LocationIndex


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

    def loadBackground(self, filename, threshold):

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

            # Clip the background grid to the size of the Study Area, which is all we will need

            start_x = int((self.sarea_bl_east - self.background_bl_east) / self.background_csize)
            start_y = int((self.sarea_bl_north - self.background_bl_north) / self.background_csize)
            end_x = int(start_x + (self.sarea_tr_east - self.sarea_bl_east) / self.background_csize)
            end_y = int(start_y + (self.sarea_tr_north - self.sarea_bl_north) / self.background_csize)

            # include an extra background grid row and column to capture locations along the very top and right
            end_x += 1
            end_y += 1

            # NOTE a numpy array is rows of columns [rownum - Y][colnum - X]
            self.background_array = self.background_array[start_y:end_y, start_x:end_x]

            # update the background header dictionary for when this is saved, and update the instance vars
            self.background_header['ncols'] = end_x - start_x
            self.background_header['nrows'] = end_y - start_y
            self.background_header['xllcorner'] = self.sarea_bl_east
            self.background_header['yllcorner'] = self.sarea_bl_north

            self.background_bl_east = int(self.background_header['xllcorner'])
            self.background_bl_north = int(self.background_header['yllcorner'])
            self.background_rows = int(self.background_header['nrows'])
            self.background_cols = int(self.background_header['ncols'])

            # background top right easting = bottom left easting + (cols*cellsize)
            self.background_tr_east = self.background_bl_east \
                                      + (self.background_cols * self.background_csize)

            # background top right northing = bottom left northing + (rows*cellsize)
            self.background_tr_north = self.background_bl_north \
                                       + (self.background_rows * self.background_csize)

            logging.info('Background array clipped to the dimensions of Analysis Area')

            logging.info('\nClipped Background Header:')
            logging.info('  BL Easting  (.projParams.background_bl_east):  ' + str(self.background_bl_east))
            logging.info('  BL Northing (.projParams.background_bl_north): ' + str(self.background_bl_north))
            logging.info('  TR Easting  (.projParams.background_tr_east):  ' + str(self.background_tr_east))
            logging.info('  TR Northing (.projParams.background_tr_north): ' + str(self.background_tr_north))
            logging.info('  Rows        (.projParams.background_rows):     ' + str(self.background_rows))
            logging.info('  Cols        (.projParams.background_cols):     ' + str(self.background_cols))

            # populate a background values dictionary of arrays with non-zero cells
            #   (and filter values to within study area)
            self.background_data = {}
            self.background_data['eastings'] = []
            self.background_data['northings'] = []
            self.background_data['XY'] = []
            self.background_data['bgval'] = []

            halfcell = int(self.aarea_csize / 2)
            # bg_histogram = {}
            out_of_range = 0

            # Y, X  is row, col
            for (Y, X) in np.ndindex(self.background_array.shape):
                val = self.background_array[Y, X]
                # not being used, but uncomment if we would like a record of how many of each unique value:
                #if val not in bg_histogram.keys():
                #    bg_histogram[val] = 1
                #else:
                #    bg_histogram[val] += 1

                # threshold for inclusion in list of relevant grid cells
                #   if we use 0 this list is potentially huge
                #   if, alternatively, we use 0.0001, much quicker and v v v slightly less accurate.
                if val > threshold:
                    bg_E = self.sarea_bl_east + self.aarea_csize * X + halfcell
                    bg_N = self.sarea_bl_north + self.aarea_csize * Y + halfcell
                    if bg_E >= self.sarea_bl_east and bg_E <= self.sarea_tr_east \
                            and bg_N >= self.sarea_bl_north and bg_N <= self.sarea_tr_north:
                        self.background_data['eastings'].append(bg_E)
                        self.background_data['northings'].append(bg_N)
                        self.background_data['XY'].append((X, Y))
                        self.background_data['bgval'].append(val)
                    else:
                        out_of_range += 1

            logging.info('  Non-zero background cells within Study Area: ' \
                         + str(len(self.background_data['eastings'])) + '  Out of range: ' + str(out_of_range))

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

                for row in range(len(ts_times)):  # now go through every row

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
                col_DataStart = int(header.iloc[3, 1])  # Start row of data
                col_DataRows = int(header.iloc[3, 3])   # Total number of data rows

                col_ID = int(header.iloc[4, 1])   # Unique ID
                col_E = int(header.iloc[5, 1])    # Easting
                col_N = int(header.iloc[6, 1])    # Northing
                col_Pop = int(header.iloc[7, 1])  # Population
                col_PopSubGroups = int(header.iloc[8, 1])  # Population subgroups count

                lst = header.iloc[9, 1:col_PopSubGroups+1].to_list()  # Default subgroup population percentages
                def_PopSubGroups = [int(elem) for elem in lst]

                lst = header.iloc[9, col_PopSubGroups + 1:2*col_PopSubGroups+1].to_list()  # Population subgroup cols
                col_PopSubGroupsCols = [int(elem) for elem in lst]

                lst = header.iloc[16, 1:col_PopSubGroups+1].to_list()  # Default mobility values
                def_MobSubGroups = [((100 - float(elem)) / 100) for elem in lst]

                lst = header.iloc[16, col_PopSubGroups + 1:2*col_PopSubGroups+1].to_list()  # Mobility subgroup cols
                col_MobSubGroupsCols = [int(elem) for elem in lst]

                col_LD = int(header.iloc[12, 2])  # Local dispersion col
                def_LD = int(header.iloc[12, 1])  # Local dispersion default

                # create a Pandas dataframe and miss out the specified header

                origin_df = pd.read_csv(filename, header=[col_DataStart-1], nrows=col_DataRows, low_memory=False)

                # filter the data to within the Study area
                # note that we include (>=) BL, and exclude (<) TR to avoid overflow
                #  may need to reduce further by csize / 2 to ensure this due to rounding
                name_E = origin_df.columns.values[col_E-1]    # Easting column name
                name_N = origin_df.columns.values[col_N - 1]  # Northing column name
                bbquery = '{} >= {} & {} <= {} & {} >= {} & {} <= {}'.format(
                    name_E, self.sarea_bl_east, name_E, self.sarea_tr_east,
                    name_N, self.sarea_bl_north, name_N, self.sarea_tr_north)

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
                            self.origin_data['subgroups_data'][columnName].append(def_PopSubGroups[idx] / 100)
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
                    # match the pop ageband array name, perhaps we should forcefully use the same index ageband name?
                    ageband = columnName.replace('Mob_', '')
                    self.origin_data['subgroups_mob'][ageband] = []
                    for mob in columnData.to_list():  # loop through mobility data groups, apply default value
                        if self.is_number(mob):
                            self.origin_data['subgroups_mob'][ageband].append((100 - int(mob)) / 100)
                        else:
                            self.origin_data['subgroups_mob'][ageband].append(def_MobSubGroups[idx])  # % immobile
                    idx += 1

                # read in the local dispersion data

                self.origin_data['LD'] = []
                for val in csvData.iloc[:, col_LD - 1].to_list():
                    if self.is_number(val):
                        self.origin_data['LD'].append(val)
                    else:
                        self.origin_data['LD'].append(def_LD)

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

                for origin in range(len(self.origin_data['eastings'])):
                    orig_E = self.origin_data['eastings'][origin]
                    orig_N = self.origin_data['northings'][origin]
                    orig_X = int((orig_E - self.background_bl_east) / csize)
                    orig_Y = int((orig_N - self.background_bl_north) / csize)
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

                # create an index for quick access to Origin locations
                # used in modelRun if area is large and Local Dispersion
                logging.info('\n  Populating Origin Location Index...')
                self.originLocationIndex = LocationIndex(self, self.origin_data)

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
                    dest_data['Filename'] = dest_file

                    header = pd.read_csv(filename, index_col=False, header=None, nrows=23)

                    # identify locations of various data columns and their defaults
                    col_DataStart = int(header.iloc[3, 1])  # Start row of data
                    col_DataRows = int(header.iloc[3, 3])   # Total number of data rows

                    col_ID = int(header.iloc[4, 1])   # Unique ID
                    col_E = int(header.iloc[5, 1])    # Easting
                    col_N = int(header.iloc[6, 1])    # Northing
                    col_Pop = int(header.iloc[7, 1])  # Population
                    col_PopSubGroups = int(header.iloc[8, 1])  # Population subgroups count

                    lst = header.iloc[9, 1:col_PopSubGroups + 1].to_list()  # Default subgroup population percentages
                    def_PopSubGroups = [int(elem) for elem in lst]

                    lst = header.iloc[9, col_PopSubGroups + 1:2 * col_PopSubGroups + 1].to_list()  # Pop subgroup cols
                    col_PopSubGroupsCols = [int(elem) for elem in lst]

                    col_TimeProfile = int(header.iloc[10, 2])  # Time profile column
                    def_TimeProfile = header.iloc[10, 1].upper()  # Default (UPPER to always match timeseries)

                    col_LD = int(header.iloc[12, 2])  # Local dispersion col
                    def_LD = int(header.iloc[12, 1])  # Default LD

                    col_WAD = int(header.iloc[13, 2])  # WAD col
                    def_WAD = header.iloc[13, 1]  # Default WAD

                    lst = header.iloc[14, 1:].to_list()  # Major Flows columns
                    col_MajorFlowCols = []
                    for elem in lst:
                        if self.is_number(elem):
                            col_MajorFlowCols.append(int(elem))

                    # create a Pandas dataframe and miss out the specified header

                    dest_df = pd.read_csv(filename, header=[col_DataStart-1], nrows=col_DataRows, low_memory=False)

                    # filter the data to within the Study area
                    # note that we include (>=) BL, and exclude (<) TR to avoid overflow
                    #  may need to reduce by csize / 2 to ensure this
                    name_E = dest_df.columns.values[col_E-1]    # Easting column name
                    name_N = dest_df.columns.values[col_N - 1]  # Northing column name
                    bbquery = '{} >= {} & {} <= {} & {} >= {} & {} <= {}'.format(
                        name_E, self.sarea_bl_east, name_E, self.sarea_tr_east,
                        name_N, self.sarea_bl_north, name_N, self.sarea_tr_north)

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

                    for dest in range(len(dest_data['eastings'])):
                        dest_E = dest_data['eastings'][dest]
                        dest_N = dest_data['northings'][dest]
                        dest_X = int((dest_E - self.background_bl_east) / csize)
                        dest_Y = int((dest_N - self.background_bl_north) / csize)
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
                            if self.is_number(pop):
                                dest_data['subgroups_data'][columnName].append(pop)
                            else:
                                dest_data['subgroups_data'][columnName].append(def_PopSubGroups[idx] / 100)

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
                            dest_data['time_profiles'].append(def_TimeProfile)
                        else:
                            dest_data['time_profiles'].append(val)

                    # read in the local dispersion data

                    dest_data['LD'] = []
                    for val in csvData.iloc[:, col_LD - 1].to_list():
                        if self.is_number(val):
                            dest_data['LD'].append(val)
                        else:
                            dest_data['LD'].append(def_LD)

                    # read in the WAD data
                    wads = []
                    wad_default = def_WAD.split('|')
                    for val in csvData.iloc[:, col_WAD - 1].to_list():
                        if not val:
                            wads.append(wad_default)
                        else:
                            wads.append(val.split('|'))

                    dest_data['WAD'] = self.convertWADs(wads)

                    logging.info('    WADs: ' + str(dest_data['WAD'][0:5]))
                    logging.info('      First tuple/radius/percent: '
                                 + str(dest_data['WAD'][0][0]) + ' '
                                 + str(dest_data['WAD'][0][0][0]) + ' '
                                 + str(dest_data['WAD'][0][0][1]))

                    # read in the Major Flows data
                    major_flows = csvData.iloc[:, min(col_MajorFlowCols) - 1:max(col_MajorFlowCols)]

                    # read in each major flow column, add to major_flows dictionary
                    #   if none provided, we could leave out of the dictionary entirely?
                    #   MF columns correspond to the Population ageband subgroups
                    dest_data['major_flows'] = {}
                    idx = 0
                    for (columnName, columnData) in major_flows.iteritems():  # loop through the columns
                        if idx >= len(dest_data['subgroup_names']):
                            break  # no more agebands, ignore the remaining MF columns
                        ageband = dest_data['subgroup_names'][idx]
                        dest_data['major_flows'][ageband] = []
                        mf_count = 0
                        for elem in columnData.to_list():
                            if elem and str(elem) != 'nan':
                                dest_data['major_flows'][ageband].append(self.convertMajorFlows(elem))
                                mf_count += 1
                            else:
                                dest_data['major_flows'][ageband].append(None)
                                #dest_data['major_flows'][columnName].append(self.convertMajorFlows('00MSNE>5|00MSMY>5|00MSNA>2'))
                        if mf_count > 0:
                            logging.info('    Major flow {} with {} values'.format(ageband, mf_count))
                        idx += 1


                    # append to the destination array
                    self.destination_data.append(dest_data)

            except IOError as e:
                logging.error(e)

        logging.info('\n  Total Destination datasets / rows (.projParams.destination_data): '
                     + str(len(self.destination_data)) + ' / ' + str(total_rows))

        # go across to numpy array at this point? 
        # dest_array = np.column_stack([dest_east, dest_north, dest_wad])

    def convertWADs(self, wad_list):
        # take a list of WADs as strings, return a new list of lists of tuples (now lists, so we can add to them):
        # amount and list are appended to hold totals and indexes of origins, background cells, etc. later
        # [  [ [radius, percent, amount, list], ... ], [ ... ] ... ]
        # to optimise the model run and grid creation, we will store radius squared

        # Empty list of all WADs for ALL destinations
        dest_wad = []

        for each_row in wad_list:  # each row of data

            wad_pairs_list = []
            for wad_pairs in each_row:  # each pair of radii/perc (or the final perc)

                if '>' in wad_pairs:
                    pair = list(wad_pairs.split('>'))
                    radius_sq = int(pair[0]) ** 2
                    perc = int(pair[1])
                else:
                    perc = int(wad_pairs)
                    # we don't have a radius, but we can calculate a radius that delivers a pop density
                    # consistent with what is inside the last radius:
                    #   new radius = square root ( last radius squared / percentage inside last radius)
                    #radius = math.sqrt(radius ** 2 / ((100 - perc) / 100))
                    radius_sq = radius_sq / ((100 - perc) / 100)

                wad_pair = [ radius_sq, perc, 0, [] ]  # NOT a tuple, use array so we can append to it later
                wad_pairs_list.append(wad_pair)

            dest_wad.append(wad_pairs_list)

        return dest_wad

    def convertMajorFlows(self, mfstr):
        # take a list of Major Flows (mfstr) as a string, return a list of tuples
        # [  [ (IDcode, percent), ... ], [ ... ] ... ]
        # Note we are working with a single mf list, not the whole lot of them as per WADs

        # Empty list of all WADs for ALL destinations
        mf_list = []

        for mf_pair in mfstr.split('|'):  # each mf

            if '>' in mf_pair:
                pair = list(mf_pair.split('>'))
                idcode = pair[0]
                perc = int(pair[1])
                mf_item = (idcode, perc)  # MF tuple
                mf_list.append(mf_item)

        return mf_list

    def is_number(self, s):
        # Returns True is string is a number
        try:
            f = float(s)
            if math.isnan(f):
                return False
            else:
                return True
        except ValueError:
            return False
