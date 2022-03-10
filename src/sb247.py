#
# Python version of Surface Builder 24/7
#
# Jan 2022
# GeoData Institute
# University of Southampton
# on behalf of ONS

# Developed to work with Python version 3.8
#   should work with most versions above 3

# Import core modules

import logging

from projectParams import ProjectParams
from modelRun import ModelRun

SPACER = '\n____ '
SPACER2 = ' ____________\n'


# The main SB247 Class, the container for all functionality and data

class SB247:

    def __init__(self, projDir):

        self.projParams = ProjectParams()
        self.modelRun = None

        if projDir.endswith('/'):
            self.projDir = projDir
        else:
            self.projDir = projDir + '/'

        logging.info(SPACER + 'SB247 Class instantiated, project directory: ' + self.projDir + SPACER2)

    def loadProjectParamsFromFile(self, filename):

        logging.info(SPACER + 'Loading Project Params from file: ' + self.projDir + filename + ' ...' + SPACER2)

        self.projParams.loadFromFile(self.projDir + filename)

        self.projParams.outputProjectParams()

    def loadProjectParamsFromDict(self, dictionary):

        logging.info(SPACER + 'Loading Project Params from Dictionary...' + SPACER2)

        self.projParams.loadFromDict(dictionary)

        self.projParams.outputProjectParams()

    def calcAreaCoords(self):

        logging.info(SPACER + 'Calculating the Analysis area and Study area coords...' + SPACER2)

        # Analysis area params that we use to calculate the top right of bbox
        # (bl = bottom left, tr = top right)

        self.projParams.aarea_bl_east = int(self.projParams.analysisarray[0])
        self.projParams.aarea_bl_north = int(self.projParams.analysisarray[1])
        self.projParams.aarea_rows = int(self.projParams.analysisarray[2])
        self.projParams.aarea_cols = int(self.projParams.analysisarray[3])
        self.projParams.aarea_csize = int(self.projParams.analysisarray[4])

        logging.info('Analysis Area:')
        logging.info('  BL Easting  (.projParams.aarea_bl_east):  ' + str(self.projParams.aarea_bl_east))
        logging.info('  BL Northing (.projParams.aarea_bl_north): ' + str(self.projParams.aarea_bl_north))
        logging.info('  Rows        (.projParams.aarea_rows):     ' + str(self.projParams.aarea_rows))
        logging.info('  Cols        (.projParams.aarea_cols):     ' + str(self.projParams.aarea_cols))
        logging.info('  Cellsize    (.projParams.aarea_csize):    ' + str(self.projParams.aarea_csize))

        # analysis area top right easting = bottom left easting + (cols * cellsize)
        self.projParams.aarea_tr_east = self.projParams.aarea_bl_east + (self.projParams.aarea_cols * self.projParams.aarea_csize)

        # analysis area top right northing = bottom left northing + (rows * cellsize)
        self.projParams.aarea_tr_north = self.projParams.aarea_bl_north + (self.projParams.aarea_rows * self.projParams.aarea_csize)

        logging.info('  TR Easting  (.projParams.aarea_tr_east):  ' + str(self.projParams.aarea_tr_east))
        logging.info('  TR Northing (.projParams.aarea_tr_north): ' + str(self.projParams.aarea_tr_north))
        logging.info('  Width: {}m  Height: {}m'.format(
            self.projParams.aarea_tr_east-self.projParams.aarea_bl_east,
            self.projParams.aarea_tr_north-self.projParams.aarea_bl_north))

        aarea_buffer = self.projParams.buffer

        # the bottom left easting of the study area is the analysis area bottom left minus buffer
        self.projParams.sarea_bl_east = self.projParams.aarea_bl_east - aarea_buffer
        self.projParams.sarea_bl_north = self.projParams.aarea_bl_north - aarea_buffer

        # the top right easting of the study area is the analysis area top right plus buffer
        self.projParams.sarea_tr_east = self.projParams.aarea_tr_east + aarea_buffer
        self.projParams.sarea_tr_north = self.projParams.aarea_tr_north + aarea_buffer

        logging.info('\nStudy Area:')
        logging.info('  BL Easting  (.projParams.sarea_bl_east):  ' + str(self.projParams.sarea_bl_east))
        logging.info('  BL Northing (.projParams.sarea_bl_north): ' + str(self.projParams.sarea_bl_north))
        logging.info('  TR Easting  (.projParams.sarea_tr_east):  ' + str(self.projParams.sarea_tr_east))
        logging.info('  TR Northing (.projParams.sarea_tr_north): ' + str(self.projParams.sarea_tr_north))
        logging.info('  Width: {}m  Height: {}m'.format(
            self.projParams.sarea_tr_east-self.projParams.sarea_bl_east,
            self.projParams.sarea_tr_north-self.projParams.sarea_bl_north))

    def geometryChecks(self):

        logging.info(SPACER + 'Coordinate Geometry checks...' + SPACER2)

        # Determine if study area is fully inside the background

        if (self.projParams.background_bl_east > self.projParams.sarea_bl_east
            or self.projParams.background_bl_north > self.projParams.sarea_bl_north) \
                and (self.projParams.background_tr_east < self.projParams.sarea_tr_east
                     or self.projParams.background_tr_north < self.projParams.sarea_tr_north):
            logging.error('  Background does not cover study area, please check the extent.')
        else:
            logging.info('  Background covers the study area.')

        # check to see if origin is wholly within the study area

        if (self.projParams.origin_eastings_min < self.projParams.sarea_bl_east
                or self.projParams.origin_northings_min < self.projParams.sarea_bl_north):

            logging.info('  One of the origin centroids is beyond the bottom left extent of the study area.')
        elif (self.projParams.origin_eastings_max > self.projParams.sarea_tr_east
              or self.projParams.origin_northings_max > self.projParams.sarea_tr_north):

            logging.info('  One of the origin centroids is beyond the top right extent of the study area.')
        else:

            logging.info('  Centroids are all within the study area.')

    def loadBackgroundFromFile(self, filename):

        logging.info(SPACER + 'Loading the Background file ' + self.projDir + filename + ' ...' + SPACER2)

        self.projParams.loadBackground(self.projDir + filename)

    def loadTimeSeriesFromFile(self, filename):

        logging.info(SPACER + 'Loading the Timeseries file ' + self.projDir + filename + ' ...' + SPACER2)

        self.projParams.loadTimeSeries(self.projDir + filename)

    def loadOriginFromFile(self, filename):

        logging.info(SPACER + 'Loading Origin data from file: ' + self.projDir + filename + ' ...' + SPACER2)

        self.projParams.loadOrigin(self.projDir + filename)

    def loadDestinationData(self, pathname):

        logging.info(SPACER + 'Loading Destination data...' + SPACER2)

        if pathname.endswith('/'):
            self.projParams.loadDestFiles(self.projDir + pathname)
        else:
            self.projParams.loadDestFiles(self.projDir + pathname + '/')

    def runSBModel(self, ageBand, runDate, runTime,
                 destination_sample_rate = 1,
                 origin_sample_rate = 1):

        logging.info(SPACER + 'Running the model...' + SPACER2)

        self.modelRun = ModelRun(ageBand, runDate, runTime,
                                 destination_sample_rate, origin_sample_rate)

        # sample_rate optional parameters are for testing model runs more quickly by sampling destinations / origins
        # not sensical for real data modelling

        # pass this object as a parameter for access to instance data

        self.modelRun.runModel(self)

    def createGridData(self):

        logging.info(SPACER + 'Creating grid data from model outputs...' + SPACER2)

        self.modelRun.createGridData(self)

    def saveGridData(self, file_prefix):

        logging.info(SPACER + 'Saving grid data to files...' + SPACER2)

        self.modelRun.saveGridData(self, file_prefix)