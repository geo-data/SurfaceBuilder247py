#!/usr/bin/env python3
#
# Python version of Surface Builder 24/7
#
# Jan 2022
# GeoData Institute
# University of Southampton
# on behalf of ONS

# unit_tests.py

#   See https://docs.python.org/3/library/unittest.html

import unittest
import datetime
import hashlib
import os

from sb247 import SB247

pre = 'Test: '
post = ' failed'


class MyTestCases(unittest.TestCase):

    def test_1_proj_params_file_function(self):
        # Test values from a project parameters file

        print("test_1_proj_params_file_function: read in values")

        # sb.loadProjectParamsFromFile('SessionParas/Bath_2011_0200_OTT_Paras.txt')
        proj_dict = {
            'analysisarray': [373000,160000,40,40,200],
            'buffer': 8000,
            'background': 'rastp1_monfri_00_06_2011.txt',
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

        msg = pre + 'Analysis Array contains expected values' + post
        assert sb.projParams.analysisarray == [373000, 160000, 40, 40, 200], msg

        msg = pre + 'Buffer = 8000' + post
        assert sb.projParams.buffer == 8000, msg

        msg = pre + 'Background filename is rastp1_monfri_00_06_2011.txt' + post
        assert sb.projParams.background == 'rastp1_monfri_00_06_2011.txt', msg

        msg = pre + 'Timeseries filename is TimeSeries.xls' + post
        assert sb.projParams.timeseries == 'TimeSeries.xls', msg

        msg = pre + 'Origin filename is Origin_Eng_OTT_2011.csv' + post
        assert sb.projParams.origin == 'Origin_Eng_OTT_2011.csv', msg

        msg = pre + 'Destination Array contains expected values' + post
        assert sb.projParams.destarray == ['Dest_Eng_Accom_OTT_2011.csv', 'Dest_Eng_Agri+Fish_OTT_2011.csv',
                                           'Dest_Eng_Education_UniPGR_OTT_2011.csv', 'Dest_Eng_Healthcare_OTT_2011.csv',
                                           'Dest_Eng_Mine+Transp_OTT_2011.csv', 'Dest_Eng_Public+Office_OTT_2011.csv',
                                           'Dest_Eng_Retail+Arts_OTT_2011.csv', 'Dest_Eng_Service_OTT_2011.csv'], msg

    def test_2_calc_area_coordinates(self):
        # Check Analysis area and Study area coords

        print("test_2_calc_area_coordinates: check calculated values")

        sb.calcAreaCoords()

        # We can just check calculated values, which will be wrong if anything else is!

        msg = pre + 'Analysis Area TR Easting = 381000' + post
        assert sb.projParams.aarea_tr_east == 381000, msg

        msg = pre + 'Analysis Area TR Northing = 160000' + post
        assert sb.projParams.aarea_tr_north == 168000, msg

        msg = pre + 'Study Area TR Easting = 389000' + post
        assert sb.projParams.sarea_tr_east == 389000, msg

        msg = pre + 'Study Area TR Northing = 176000' + post
        assert sb.projParams.sarea_tr_north == 176000, msg

    def test_3_background_file(self):
        # Check values from the background file header

        print("test_3_background_file: check background file header values")

        sb.loadBackgroundFromFile('BckGrnds/' + sb.projParams.background)

        # We can just check calculated values, which will be wrong if anything else is!

        msg = pre + 'Background file TR Easting = 389200' + post
        assert sb.projParams.background_tr_east == 389200, msg

        msg = pre + 'Background file TR Northing = 176200' + post
        assert sb.projParams.background_tr_north == 176200, msg

    def test_4_origin_file(self):
        # Check values from the origin file

        print("test_4_origin_file: check origin file values")

        sb.loadOriginFromFile('Origins/' + sb.projParams.origin)

        # Check count and calculated values, which will be wrong if anything else is!

        msg = pre + 'Origin file Easting/Northing count within study area  = 998' + post
        assert len(sb.projParams.origin_data['northings']) == 988 \
               and len(sb.projParams.origin_data['northings']) == 988, msg

        msg = pre + 'Origin last row ID = E00173096' + post
        assert sb.projParams.origin_data['ID'][987] == 'E00173096', msg

        msg = pre + 'Origin file Min Easting = 365005' + post
        assert sb.projParams.origin_eastings_min == 365006, msg

        msg = pre + 'Origin file Max Northing = 175950' + post
        assert sb.projParams.origin_northings_max == 175950, msg

        msg = pre + 'Origin file last row Pop = 342' + post
        assert sb.projParams.origin_data['pop_data'][987] == 342, msg

        msg = pre + 'Origin file final pop subgroup name = OV65' + post
        assert sb.projParams.origin_data['subgroup_names'][5] == 'OV65', msg

        msg = pre + 'Origin file final pop subgroup final value = 1.169591' + post
        assert sb.projParams.origin_data['subgroups_data']['OV65'][987] == 1.169591, msg

    def test_5_destination_data(self):
        # Check values from the destination data files

        print("test_5_destination_data: check destination files values")

        sb.loadDestinationData('Dests/')

        # Check count and calculated values, which will be wrong if anything else is!

        dest_check = sb.projParams.destination_data[7]

        msg = pre + 'Destination final file time profile = TS01.WKG.2SERV' + post
        assert dest_check['time_profiles'][0] == 'TS01.WKG.2SERV', msg

        msg = pre + 'Destination final file OA/Easting/Northing/WAD count = 651' + post
        assert len(dest_check['ID']) == 306 \
               and len(dest_check['eastings']) == 306 \
               and len(dest_check['northings']) == 306 \
               and len(dest_check['WAD']) == 306, msg

        msg = pre + 'Destination final file last row OA = E33050776' + post
        assert dest_check['ID'][305] == 'E33050776', msg

        msg = pre + 'Destination final file last row Easting = 361601' + post
        assert dest_check['eastings'][305] == 382831, msg

        msg = pre + 'Destination final file last row Northing = 174207' + post
        assert dest_check['northings'][305] == 154797, msg

        msg = pre + 'Destination final file last row Pop = 16' + post
        assert dest_check['pop_data'][305] == 16, msg

        msg = pre + 'Destination final file final pop subgroup name = OV65' + post
        assert dest_check['subgroup_names'][5] == 'OV65', msg

        msg = pre + 'Destination final file final pop subgroups final value = 5.070422535' + post
        assert dest_check['subgroups_data']['OV65'][305] == 5.070422535, msg

        msg = pre + 'Destination final file WAD last row, second tuple radius squared = 2000^2 ' + post
        assert dest_check['WAD'][305][1][0] == 2000 ** 2, msg

        msg = pre + 'Destination file WAD last row, third tuple percent = 12 ' + post
        assert dest_check['WAD'][305][2][1] == 12, msg

    def test_6_timeseries_file(self):
        # Check count and value from the timeseries file

        print("test_6_timeseries_file: check timeseries file values")

        sb.loadTimeSeriesFromFile('TimeSeries/' + sb.projParams.timeseries)

        msg = pre + 'Timeseries file Key count = 18' + post
        assert len(sb.projParams.timeseries_data) == 18, msg

        msg = pre + 'Timeseries file first time value = 00:00:00 ' + post
        first_time = next(iter(sb.projParams.timeseries_data['TS01.WKG.6Offb'.upper()]['OnSite']))
        assert str(first_time) == '00:00:00', msg

        msg = pre + 'Timeseries first final weight value = 2.68065 ' + post
        assert sb.projParams.timeseries_data['TS01.WKG.6Offb'.upper()]['OnSite'][first_time] == 2.68065, msg

    def test_7_run_model(self):
        # run the model, check the ouptput grids

        print("test_7_run_model: check model output gridded values")

        ageband = 'OV65'
        run_date = datetime.date(2020, 2, 25)  # currently not used
        run_time = datetime.time(9, 30, 0)

        sb.runSBModel(ageband, run_date, run_time,
                      destination_sample_rate = 20,
                      origin_sample_rate = 5)

        sb.createGridData()

        msg = pre + 'Sum of all remaining origins locally dispersed = 53639' + post
        assert int(sum(sum(sb.modelRun.grid_origins_remain_LD))) == 53639, msg

        msg = pre + 'Sum of all on site destinations locally dispersed = 153' + post
        assert int(sum(sum(sb.modelRun.grid_dest_onSite_LD))) == 153, msg

        msg = pre + 'Sum of all in travel destinations = 13' + post
        assert int(sum(sum(sb.modelRun.grid_dest_inTravel))) == 13, msg

        msg = pre + 'Sum of all immobile origins locally dispersed = 689' + post
        assert int(sum(sum(sb.modelRun.grid_origins_immob_LD))) == 689, msg

    def test_8_save_data(self):
        # Save the data (temporarily), check the file checksums!

        print("test_8_save_data: Save and check output data")

        sb.saveOutputData('Results/Unit_tests_')

        ck1 = self.file_checksum(sb.projDir + 'Results/Unit_tests_dest_inTravel.asc')[:16]
        ck2 = self.file_checksum(sb.projDir + 'Results/Unit_tests_dest_onSite_LD.asc')[:16]
        ck3 = self.file_checksum(sb.projDir + 'Results/Unit_tests_origins_immob_LD.asc')[:16]
        ck4 = self.file_checksum(sb.projDir + 'Results/Unit_tests_origins_remain_LD.asc')[:16]
        ck5 = self.file_checksum(sb.projDir + 'Results/Unit_tests_results_LD.csv')[:16]
        #print (ck1 + ck2 + ck3 + ck4 + ck5)

        ckall = 'c68e4c91b2c93d5e8939efef7972d0db74ba30ae216ba442262e4de5dd47ccbf1352dba2860ede0e'
        msg = pre + 'Check all saved file checksums are as expected' + post
        assert ck1 + ck2 + ck3 + ck4 + ck5 == ckall, msg

        # remove the files
        os.remove(sb.projDir + 'Results/Unit_tests_dest_inTravel.asc')
        os.remove(sb.projDir + 'Results/Unit_tests_dest_onSite_LD.asc')
        os.remove(sb.projDir + 'Results/Unit_tests_origins_immob_LD.asc')
        os.remove(sb.projDir + 'Results/Unit_tests_origins_remain_LD.asc')
        os.remove(sb.projDir + 'Results/Unit_tests_results_LD.csv')

    def file_checksum(self, filename):
        try:
            with open(filename, 'r') as file_opened:
                content = file_opened.read().encode('ascii', 'ignore')
                hl = hashlib.sha256()
                hash = hl.update(content)
                digest = hl.hexdigest()
                return digest

        except IOError as e:
            assert 'File read error', e


# Executed when the program runs

if __name__ == "__main__":
    print('Starting SB247 Unit tests...')

    # create a single Surface Builder object to avoid repeatedly starting from scratch

    sb = SB247('./Data')

    # run all tests

    unittest.main()

    # set up a test suite to run subsets of tests
    """
    def suite():
        suite = unittest.TestSuite()
        suite.addTest(MyTestCases('test_1_proj_params_file_function'))
        suite.addTest(MyTestCases('test_2_calc_area_coordinates'))
        suite.addTest(MyTestCases('test_3_background_file'))
        suite.addTest(MyTestCases('test_4_origin_file'))
        suite.addTest(MyTestCases('test_5_destination_data'))
        suite.addTest(MyTestCases('test_6_timeseries_file'))
        suite.addTest(MyTestCases(''))
        suite.addTest(MyTestCases(''))
        return suite

    
    # run just 1 test - disable the test all above and enable these 2 lines
    runner = unittest.TextTestRunner()    
    runner.run(suite())
    """
