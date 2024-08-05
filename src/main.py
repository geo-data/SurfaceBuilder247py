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
    """
    Demo of a model run for a single age group
    """

    try:

        # instantiate an object of the SB247 class
        #   project directory parameter will prefix all file location references
        sb = SB247('./Data')

        # instance variables of the class are accessible as required
        print('Project directory: ' + sb.projDir)

        # Project parameters can be loaded from a Dictionary or (VB compatible) file

        # Load Project parameters from a dictionary

        proj_dict = {
            'analysisarray': [373000, 160000, 40, 40, 200],  # BL_E, BL_N, nrows, ncols, cellsize
            'buffer': 8000,
            'background': 'rastp6_monfri_10_16_2011.txt',
            'timeseries': 'TimeSeries.xls',
            'origin': 'Origin_Eng_OTT_2011.csv',
            'destarray': ['Dest_Eng_Accom_OTT_2011.csv',
                          'Dest_Eng_Agri+Fish_OTT_2011.csv',
                          'Dest_Eng_Education_UniPGR_OTT_2011.csv',
                          'Dest_Eng_Healthcare_OTT_2011.csv',
                          'Dest_Eng_Mine+Transp_OTT_2011.csv',
                          'Dest_Eng_Public+Office_OTT_2011.csv',
                          'Dest_Eng_Retail+Arts_OTT_2011.csv',
                          'Dest_Eng_Service_OTT_2011.csv']            
        }

        sb.loadProjectParamsFromDict(proj_dict)

        # Load Project parameters from a file

        # sb.loadProjectParamsFromFile('SessionParas/Bath_2011_0200_OTT_Paras.txt')

        # Calculate the Analysis area and Study area coords

        sb.calcAreaCoords()

        # Load the Background from an Ascii grid file
        #   Optional threshold parameter selection of lowest value to use for inTravel dispersion
        #   0 includes all data, 0.0001 will speed things up for large areas, but can cause problems with dispersing population
        #    in areas with few transport links

        sb.loadBackgroundFromFile('BckGrnds/' + sb.projParams.background,
                                  threshold = 0)

        # Load the origin csv, get max/min coords

        sb.loadOriginFromFile('Origins/' + sb.projParams.origin)

        # Run some basic checks on the coordinate geometry

        sb.geometryChecks()

        # Load all the destination csv files into a combined list / pandas dataframe

        sb.loadDestinationData('Dests/')

        # load TimeSeries from an Excel file

        sb.loadTimeSeriesFromFile('TimeSeries/' + sb.projParams.timeseries)

        # Everything is now loaded, run the model with some sample values

    except Exception as err:
        logging.error('\nError loading data: ' + str(err))
        exit(1)

    try:
        ageband = 'OV65'  # OV65 or 18_64
        run_date = datetime.date(2020, 2, 25)  # currently not used
        run_time = datetime.time(14, 0, 0)     # 2pm

        sb.runSBModel(ageband, run_date, run_time,
                      # optional testing parameters, for quick model run, remove or set to 1 for normal operation
                      destination_sample_rate = 1,   # process 1 in N rows of each destination dataset
                      origin_sample_rate = 1         # process 1 in N rows of origin data
                      )

        sb.createGridData(create_non_LD = True,
                          cressman_power = 1)

        # Create a suitable path/file prefix for saving the files, based on model run parameters
        file_prefix = 'Results/Bath_2011_0200_OTT_Paras_{}_{:02}_{:02}_'.format(ageband,
                                                                                run_time.hour,
                                                                                run_time.minute)

        sb.saveOutputData(file_prefix, save_CSV_files = False, save_destination_file_grids = False)

        logging.info('\nRun successful.')
        exit(0)

    except Exception as err:
        logging.info('Problem running model: ' + str(err))



def main_allages():
    """
    Demo of a model run with looping through the age groups
    All age groups present in the origin file wil be processed.
    """

    try:

        # instantiate an object of the SB247 class
        #   project directory parameter will prefix all file location references
        sb = SB247('./Data')

        # instance variables of the class are accessible as required
        print('Project directory: ' + sb.projDir)

        # Project parameters can be loaded from a Dictionary or (VB compatible) file

        # Load Project parameters from a dictionary

        proj_dict = {
            'analysisarray': [373000, 160000, 40, 40, 200],  # BL_E, BL_N, nrows, ncols, cellsize
            'buffer': 8000,
            'background': 'rastp6_monfri_10_16_2011.txt',
            'timeseries': 'TimeSeries.xls',
            'origin': 'Origin_Eng_OTT_2011.csv',
            'destarray': ['Dest_Eng_Accom_OTT_2011.csv',
                          'Dest_Eng_Agri+Fish_OTT_2011.csv',
                          'Dest_Eng_Education_UniPGR_OTT_2011.csv',
                          'Dest_Eng_Healthcare_OTT_2011.csv',
                          'Dest_Eng_Mine+Transp_OTT_2011.csv',
                          'Dest_Eng_Public+Office_OTT_2011.csv',
                          'Dest_Eng_Retail+Arts_OTT_2011.csv',
                          'Dest_Eng_Service_OTT_2011.csv']
        }

        sb.loadProjectParamsFromDict(proj_dict)

        # Load Project parameters from a file

        # sb.loadProjectParamsFromFile('SessionParas/Bath_2011_0200_OTT_Paras.txt')

        # Calculate the Analysis area and Study area coords

        sb.calcAreaCoords()

        # Load the Background from an Ascii grid file
        #   Optional threshold parameter selection of lowest value to use for inTravel dispersion
        #   0 includes all data, 0.0001 will speed things up for large areas, but can cause problems with dispersing population
        #    in areas with few transport links

        sb.loadBackgroundFromFile('BckGrnds/' + sb.projParams.background,
                                  threshold = 0)

        # Load the origin csv, get max/min coords

        sb.loadOriginFromFile('Origins/' + sb.projParams.origin)

        # Run some basic checks on the coordinate geometry

        sb.geometryChecks()

        # Load all the destination csv files into a combined list / pandas dataframe

        sb.loadDestinationData('Dests/')

        # load TimeSeries from an Excel file

        sb.loadTimeSeriesFromFile('TimeSeries/' + sb.projParams.timeseries)

        # Everything is now loaded, run the model with some sample values

    except Exception as err:
        logging.error('\nError loading data: ' + str(err))
        exit(1)

    try:
        agebands = sb.projParams.origin_data['subgroup_names'] # Read from the origin file header
        # agebands = ['0_4','5_9','10_15','16_17','18_64','OV65']  # List of age groups to loop through

        run_date = datetime.date(2020, 2, 25)  # currently not used
        run_time = datetime.time(14, 0, 0)     # 2pm

        for ageband in agebands:
            sb.runSBModel(ageband, run_date, run_time,
                          # optional testing parameters, for quick model run, remove or set to 1 for normal operation
                          destination_sample_rate = 1,   # process 1 in N rows of each destination dataset
                          origin_sample_rate = 1         # process 1 in N rows of origin data
                          )

            sb.createGridData(create_non_LD = False,
                              cressman_power = 1)

            # Create a suitable path/file prefix for saving the files, based on model run parameters
            file_prefix = 'Results/Bath_2011_0200_OTT_Paras_{}_{:02}_{:02}_'.format(ageband,
                                                                                    run_time.hour,
                                                                                    run_time.minute)

            sb.saveOutputData(file_prefix, save_CSV_files = False, save_destination_file_grids = False)

            # Aggregate the outputs from this run of the model into a set of running totals
            sb.aggregateOutputData()

            logging.info('\n{} Run successful.'.format(ageband))

        file_prefix = 'Results/Bath_2011_0200_OTT_Paras_{}_{:02}_{:02}_'.format('Aggregated',
                                                                                run_time.hour,
                                                                                run_time.minute)

        sb.saveAggregatedOutputData(file_prefix)

        logging.info('\nRun successful for Allages.')

        exit(0)

    except Exception as err:
        logging.info('Problem running model: ' + str(err))


# Executed when the program runs

if __name__ == '__main__':
    # main() # Single age band

    main_allages() # Loop through multiple age bands
