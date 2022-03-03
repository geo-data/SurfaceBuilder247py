#!/usr/bin/env python3
#
# Python version of Surface Builder 24/7
#
# Jan 2022
# GeoData Institute
# University of Southampton
# on behalf of ONS

# Sample main.py

import logging
import datetime

from sb247 import SB247

# Set up the destination, level and style of logging (INFO, WARNING, ERROR)
#
# logging.basicConfig(filename='log_sb2472py.txt', level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
# logging.basicConfig(level=logging.WARNING, format='%(message)s')

logging.basicConfig(level=logging.INFO, format='%(message)s')


def main():
    try:

        # instantiate an object of the SB247 class
        #   project directory parameter will prefix all file location references
        sb = SB247('./Data')

        # instance variables of the class are accessible as required
        print('Project directory is... ' + sb.projDir)

        # Project parameters can be loaded from a Dictionary or (VB compatible) file

        # Load Project parameters from a dictionary
        """
        proj_dict = {
            'analysisarray': [373000,160000,40,40,200],
            'buffer': 8000,
            'background': 'rastp1_monfri_00_06_2011.txt',
            'timeseries': 'TimeSeries.xls',
            'origin': 'Origin_Eng_OTT_2011.csv',
            'destarray': ['Dest_Eng_Accom_OTT_2011.csv',
                          'Dest_Eng_Agri+Fish_OTT_2011.csv',
                          'Dest_Eng_Healthcare_OTT_2011.csv',
                          'Dest_Eng_Mine+Transp_OTT_2011.csv',
                          'Dest_Eng_Public+Office_OTT_2011.csv',
                          'Dest_Eng_Retail+Arts_OTT_2011.csv',
                          'Dest_Eng_Service_OTT_2011.csv']            
        }

        sb.loadProjectParamsFromDict(proj_dict)
        """

        # Load Project parameters from a file

        sb.loadProjectParamsFromFile('SessionParas/Bath_2011_0200_OTT_Paras.txt')

        # Calculate the Analysis area and Study area coords

        sb.calcAreaCoords()

        # Load the Background from an Ascii grid file

        sb.loadBackgroundFromFile('BckGrnds/' + sb.projParams.background)

        # Load the origin csv, get max/min coords

        sb.loadOriginFromFile('Origins/' + sb.projParams.origin)

        # Run some basic checks on the coordinate geometry

        sb.geometryChecks()

        # Load all the destination csv files into a combined list / pandas dataframe

        sb.loadDestinationData('Dests/')

        # load TimeSeries from an Excel file

        sb.loadTimeSeriesFromFile('TimeSeries/' + sb.projParams.timeseries)

        # Everything is now loaded, run the model with some sample values

        try:
            ageband = '18_64'
            run_date = datetime.date(2020, 2, 25)  # currently not used
            run_time = datetime.time(9, 40, 0)

            sb.runSBModel(ageband, run_date, run_time,
                          # optional testing parameters, for quick model run, remove or set to 1 for normal operation
                          destination_sample_rate = 20,  # process 1 in N rows of each destination dataset
                          origin_sample_rate = 1         # process 1 in N rows of origin data
                          )

            sb.createGridData()

            # Create a suitable path/file prefix for saving the files, based on model run parameters
            file_prefix = 'Results/Bath_2011_0200_OTT_Paras_' \
                          + ageband + '_' \
                          + str(run_time.hour) + '_' + str(run_time.minute) + '_'

            sb.saveGridData(file_prefix)

            logging.info('\nRun successful.')

        except Exception as err:
            logging.info('Problem running model: ' + str(err))

    except Exception as e:
        logging.error(e)


# Executed when the program runs

if __name__ == '__main__':
    main()
