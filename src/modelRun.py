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
import time
import math
import numpy as np
import copy
import itertools

from locationIndex import LocationIndex
from gridCreate import GridCreate

DEST_DEBUG_LIMIT = -1    # limit the number of rows we process (for each dest collection) for testing
ORIG_DEBUG_LIMIT = -1    # set to a number or -1 for all of them


# A class for carrying out SB247 model runs

class ModelRun:

    def __init__(self, ageBand, runDate, runTime,
                 dest_sample_rate,
                 orig_sample_rate):

        self.ageband = ageBand
        self.date = runDate
        self.time = runTime

        logging.info('  Age band (.modelRun.ageband): ' + str(self.ageband))
        logging.info('  Date     (.modelRun.date):    ' + str(self.date))
        logging.info('  Time     (.modelRun.time):    ' + str(self.time))

        # for testing model runs more quickly by sampling destinations / origins
        # not sensical for real data modelling
        self.dest_sample_rate = dest_sample_rate
        self.orig_sample_rate = orig_sample_rate

        if dest_sample_rate > 1:
            logging.info('  Sampling 1 in ' + str(dest_sample_rate) + ' destinations')

        if orig_sample_rate > 1:
            logging.info('  Sampling 1 in ' + str(orig_sample_rate) + ' origins')

    def runModel(self, sb):

        # loop through each destination collection
        loop_count = 0
        initialTime = time.time()

        # make a copy of the relevant origin populations
        self.originPopData = sb.projParams.origin_data['subgroups_pop'][self.ageband].copy()

        # transfer immobile population
        self.originPopDataImmob = []
        for origin in range(len(self.originPopData)):
            immob = self.originPopData[origin] * sb.projParams.origin_data['subgroups_mob'][self.ageband][origin]
            self.originPopDataImmob.append(immob)
            self.originPopData[origin] -= immob

        logging.info(  '\n  Immobile population removed: {}'.format(round(sum(self.originPopDataImmob),3)) )

        originInitialPop = sum(self.originPopData)  # record initial total pop in this ageband
        destIncrease = 0  # (for checking) a record of how many dest pop transfers were made

        fail_count = 0  # for reporting how many destinations couldn't be provided with all requested population
        fail_total = 0  # and how much total population wasn't provided

        # create dictionary of arrays to store all of the destination data values needed (and their grid indexes)
        # create a list of dictionaries for each destination file containing all the information for each destination
        self.destination_data_array = {}


        for destdata in sb.projParams.destination_data:

            destination_data = {}
            destination_data['inTravel'] = []
            destination_data['onSite'] = []
            destination_data['eastings'] = []  # list of Eastings
            destination_data['northings'] = []  # list of Eastings
            destination_data['XY'] = []
            destination_data['WAD'] = []
            destination_data['wadStrings'] = []
            destination_data['LD'] = []
            destination_data['hash'] = []

            logging.info('\n  New destination collection (' + destdata['Filename'] + ') ...')

            current_time_profile = ''  # only recalculate the inTravel and onSite percentages when things change

            # loop through each destination row in the collection

            for dest in range(len(destdata['pop_data'])):

                #if dest == DEST_DEBUG_LIMIT:  # fewer rows for debugging
                #    break

                if dest % self.dest_sample_rate != 0:  # sample the dest
                    continue

                # grab the population for the specified age category
                dest_pop = destdata['subgroups_pop'][self.ageband][dest]

                if dest_pop == 0: # No demand for this age group.
                    continue

                # grab the time profile and calculate the percentages (if things change)
                time_profile = destdata['time_profiles'][dest]

                if time_profile != current_time_profile:

                    (inTravel_pc, onSite_pc) = self.timeProfileLookup(sb, time_profile)

                    if inTravel_pc == 0 and onSite_pc == 0:
                        # No demand for this time period - move on to next
                        continue

                    logging.info('\n  New time profile: ' + time_profile + '\n')
                    logging.info('    dest percentages for inTravel: '
                                 + str(inTravel_pc) + '  onSite: ' + str(onSite_pc))

                    # if we need to distribute unused origin pop to more distant WADs
                    dest_inTravel_ratio = inTravel_pc / (inTravel_pc + onSite_pc)
                    dest_onSite_ratio = onSite_pc / (inTravel_pc + onSite_pc)

                    current_time_profile = time_profile  # remember for next time

                # calculate the required OnSite Pop / InTravel Pop
                dest_inTravel_pop = dest_pop * inTravel_pc / 100
                dest_onSite_pop = dest_pop * onSite_pc / 100
                dest_req_pop = dest_inTravel_pop + dest_onSite_pop
                destIncrease += dest_req_pop

                dest_E = destdata['eastings'][dest]
                dest_N = destdata['northings'][dest]

                # save values and grid index for our output grids
                destination_data['inTravel'].append(dest_inTravel_pop)
                destination_data['onSite'].append(dest_onSite_pop)
                destination_data['eastings'].append(dest_E)
                destination_data['northings'].append(dest_N)
                destination_data['XY'].append(destdata['XY'][dest])
                destination_data['WAD'].append(destdata['WAD'][dest])
                destination_data['wadStrings'].append(destdata['wadStrings'][dest])
                destination_data['LD'].append(destdata['LD'][dest])
                destination_data["hash"].append(destdata['hash'][dest])                

                # logging.info('\n    Dest ' + str(dest)
                #              + '. E: ' + str(dest_E) + ' N: ' + str(dest_N)
                #              + ' Pop: ' + str(round(dest_pop,3))
                #              + '  inTravel: ' + str(round(dest_inTravel_pop,3))
                #              + '  onSite: ' + str(round(dest_onSite_pop,3))
                #              + '  total req: ' + str(round(dest_req_pop,3)))

                dest_remove_check = 0

                # Are there any Major Flows for this dest record?
                #   if so, loop through each of them
                #      grab ID pattern and percentage
                #      calculate destination inTravel/onSite to remove based on mf percent
                #      find the list of origins matching and the total population they hold
                #      add origins to list of MF origins affected by this dest
                #          for each origin in the list
                #             calc ratio of origin pop to pop total
                #             subtract proportion from origins

                mf_list = destdata['major_flows'][self.ageband][dest]
                mf_origins_dest = []  # for checking against later
                mf_total = 0
                mf_total_inTravel = 0
                mf_total_onSite = 0
                if mf_list:
                    for mf in mf_list:
                        mf_ID = mf[0]
                        mf_pc = mf[1]
                        mf_dest_remove_inTravel = dest_inTravel_pop * mf_pc / 100  # do we need separate values?
                        mf_dest_remove_onSite = dest_onSite_pop * mf_pc / 100
                        mf_dest_remove_total = dest_req_pop * mf_pc / 100

                        mf_total += mf_dest_remove_total  # update totals
                        mf_total_inTravel += mf_dest_remove_inTravel
                        mf_total_onSite += mf_dest_remove_onSite

                        (mf_origins, mf_origin_pop) = self.mf_origin_list(sb, mf_ID)
                        mf_origins_dest.extend(mf_origins)  # add to master list
                        # now loop through each origin, removing the proportion from each
                        for origin in mf_origins:
                            origin_remove = self.originPopData[origin] / mf_origin_pop * mf_dest_remove_total
                            self.originPopData[origin] -= origin_remove

                    # reduce the remaining required destination amounts
                    dest_inTravel_pop -= mf_total_inTravel
                    dest_onSite_pop -= mf_total_onSite
                    dest_remove_check += mf_total

                    logging.info('      Major Flows: {} - Population removed {:.3f} from {} Origins'.format(len(mf_list), mf_total, len(mf_origins_dest)))

                # loop through each WAD pair (assume nearest is always first)

                # Extract the radius and perc, and then append 0 and an empty list - benchmarking suggests that this is 10x faster than deepcopy
                dest_WAD = []
                for w in destdata['WAD'][dest]:
                    dest_WAD.append(w + [0,[]])
                
                # find list of origins which are within the radius

                largest_radius = math.sqrt(dest_WAD[len(dest_WAD) - 1][0])
                potential_origins_array = sb.projParams.originLocationIndex.possible_locations(dest_E, dest_N, largest_radius)
                #logging.info('      Max radius: {:.3f} containing {} potential Origins'.format(largest_radius,
                #                                                                               len(potential_origins_array)))

                for potential_origins in potential_origins_array:
                    for origin in potential_origins:

                        #if origin == ORIG_DEBUG_LIMIT:  # fewer rows for debugging
                        #    break

                        if origin % self.orig_sample_rate != 0:  # sample the origins
                            continue

                        if origin in mf_origins_dest:  # we have Major Flowed this origin already
                            continue

                        loop_count += 1

                        orig_E = sb.projParams.origin_data['eastings'][origin]
                        orig_N = sb.projParams.origin_data['northings'][origin]

                        orig_pop = self.originPopData[origin]

                        # if the pop is zero, though unlikely, we might as well stop (continue) here?
                        #  test more later when using more destinations, appears to slow things down
                        #if orig_pop == 0:
                        #    continue

                        # pythagoras gives us the distance between origin and destination
                        dist_sq = (dest_E - orig_E) ** 2 + (dest_N - orig_N) ** 2

                        #logging.debug('      Orig ' + str(origin)
                        #             + '. E: ' + str(orig_E) + ' N: ' + str(orig_N)
                        #             + ' Pop: ' + str(round(orig_pop, 2))
                        #             + ' distance: ' + str(round(dist, 2)))

                        # which wad is this distance relevant to
                        for wad in dest_WAD:
                            if dist_sq <= wad[0]:  # within range
                                wad[2] += orig_pop  # store the total origin population
                                wad[3].append(origin)  # add the origin index to our list
                                break

                # The dest wad is now fully populated with origin indexes and origin pop total

                residue_inTravel = 0   # keep track of unmet wad pop requirement, to pass up to next radius
                residue_onSite = 0

                for wad in dest_WAD:

                    # origin_wad_pop = wad[2]
                    available_pop = wad[2]

                    # figure out how many people need pulling for this WAD
                    #rad = wad[0]
                    pc = wad[1]
                    wad_inTravel = dest_inTravel_pop * pc / 100
                    wad_onSite = dest_onSite_pop * pc / 100
                    #wad_total = wad_inTravel + wad_onSite  # to remove from these origins

                    # add to our current balance of required pop for dest
                    residue_inTravel += wad_inTravel
                    residue_onSite += wad_onSite
                    residue_total = residue_inTravel + residue_onSite

                    """
                    logging.debug('      ' + str(pc) + '%  within ' + str(rad) + 'm -> '
                                 + ' inTravel ' + str(round(wad_inTravel, 3))
                                 + ', onSite ' + str(round(wad_onSite, 3))
                                 + ' total ' + str(round(wad_total,3))
                                 + ' from origins (' + str(round(origin_wad_pop,3)) + ' available)')
                    logging.debug('              residue: '
                                 + ' inTravel ' + str(round(residue_inTravel, 3))
                                 + ', onSite ' + str(round(residue_onSite, 3))
                                 + ' total ' + str(round(residue_total,3)))
                    """

                    orig_remove_check = 0

                    if available_pop > 0:  # some origin population is available to take

                        available_origins = wad[3]

                        if available_pop > residue_total:
                            # enough to satisfy requirement fully, satisfy all residue
                            wad_remove_total = residue_total
                            # breakdowns not used, probably not needed
                            # wad_remove_inTravel = residue_inTravel
                            # wad_remove_onSite = residue_onSite
                            residue_inTravel = 0
                            residue_onSite = 0
                        else:
                            # not enough origin pop in this WAD, use it all up and update residues
                            wad_remove_total = available_pop
                            wad_remove_inTravel = wad_remove_total * dest_inTravel_ratio
                            wad_remove_onSite = wad_remove_total * dest_onSite_ratio
                            residue_inTravel -= wad_remove_inTravel
                            residue_onSite -= wad_remove_onSite
                            #logging.debug('not enough origin pop in this WAD!')

                        # let's do some removing of pops from origins
                        # code is much optimised by using a simple multiplier
                        mult = (available_pop - wad_remove_total) / available_pop
                        for origin in available_origins:
                            # go through each origin index, remove in proportion with origin pop
                            #pop = self.originPopData[origin]
                            #orig_remove_total = pop / available_pop * wad_remove_total
                            # breakdowns not used, probably not needed
                            #orig_remove_inTravel = pop / available_pop * wad_remove_inTravel
                            #orig_remove_onSite = pop / available_pop * wad_remove_onSite
                            #self.originPopData[origin] -= orig_remove_total
                            self.originPopData[origin] *= mult
                            #logging.debug('        removed '
                            #             + str(round(orig_remove_total,3))
                            #             + ' from ' + str(origin) + ' (' + str(round(pop,3)) + ')')
                            #orig_remove_check += orig_remove_total

                        available_pop -= wad_remove_total  # update total pop available
                        orig_remove_check = wad_remove_total

                        #logging.debug('        Orig remove check: ' + str(round(orig_remove_check,3)))
                        dest_remove_check += orig_remove_check

                # destination is complete, check we removed the full amount
                # if round(dest_req_pop,3) == round(dest_remove_check,3):
                if dest_req_pop - dest_remove_check < 1:
                    pass
                    # logging.info('      Dest remove check SUCCESS: {}'.format(str(round(dest_remove_check,3))))
                else:
                    logging.info('      Dest {} remove check FAIL: {} removed, {} outstanding'.format(
                        destdata['ID'][dest],
                        round(dest_remove_check, 3),
                        round(dest_req_pop - dest_remove_check, 3)))
                    fail_count += 1
                    fail_total += dest_req_pop - dest_remove_check

            #save the destination_data dict to the by-dest file dict
            self.destination_data_array[destdata['Filename']] = destination_data
        originFinalPop = sum(self.originPopData)  # record initial total pop in this ageband

        logging.info('\n  Run Complete - Loop count: ' + str(loop_count)
                     + ' in ' + str(round(time.time() - initialTime,1)) + ' seconds')

        logging.info('\n  Origin pop initial / final / diff: {} / {} / {}\n  Dest pop requested: {}'.format(
            round(originInitialPop, 3), round(originFinalPop, 3),
            round(originInitialPop - originFinalPop, 3), round(destIncrease, 3)))

        if fail_count > 0:
            logging.info('\n  {} unsatisfied destinations, {} requested'.format(fail_count, round(fail_total, 3)))

    def concatenateDestinationDataArrays(self, destdata_dict):

        concat_dict = {}

        #get all the array namess that should be in a destdata dict from the first one
        array_names = destdata_dict[next(iter(destdata_dict))].keys()

        #and combine all of the versions of each of them them into one list and put it
        # in the dict
        for array_name in array_names:

            list_of_lists = [destdata_dict[filename][array_name] for filename in destdata_dict]
            concat_dict[array_name] = list(itertools.chain(*list_of_lists))

        return concat_dict

    def timeProfileLookup(self, sb, time_profile):
        # lookup a time profile for the required time, return a tuple of inTravel and onSite percents

        if self.time in sb.projParams.timeseries_data[time_profile]['InTravel']:
            inTravel_pc = sb.projParams.timeseries_data[time_profile]['InTravel'][self.time]
        else:
            # try going backwards until we find a match, e.g 9.45 falls into 9.40-9.50 slot
            inTravel_pc = 0  # fallback value
            try_mins = 1
            while try_mins <= 60:  # up to an hour BEFORE
                try_time = self.addMins(self.time, -(try_mins))
                if try_time in sb.projParams.timeseries_data[time_profile]['InTravel']:
                    inTravel_pc = sb.projParams.timeseries_data[time_profile]['InTravel'][try_time]
                    break
                try_mins += 1

        if self.time in sb.projParams.timeseries_data[time_profile]['OnSite']:
            onSite_pc = sb.projParams.timeseries_data[time_profile]['OnSite'][self.time]
        else:
            onSite_pc = 0
            try_mins = 1
            while try_mins <= 60:
                try_time = self.addMins(self.time, -(try_mins))
                if try_time in sb.projParams.timeseries_data[time_profile]['OnSite']:
                    onSite_pc = sb.projParams.timeseries_data[time_profile]['OnSite'][try_time]
                    break
                try_mins += 1

        return(inTravel_pc, onSite_pc)


    def addMins(self, tm, mins):
        fulldate = datetime.datetime(100, 1, 1, tm.hour, tm.minute, 0)
        fulldate = fulldate + datetime.timedelta(minutes=mins)
        return fulldate.time()

    def mf_origin_list(self, sb, mf_ID):

        # find a list of origins matching this ID and the total origin population for them

        mf_list = []
        mf_pop = 0
        mf_len = len(mf_ID)

        for origin in range(len(sb.projParams.origin_data['eastings'])):
            if sb.projParams.origin_data['ID'][origin][:mf_len] == mf_ID:
                mf_list.append(origin)
                mf_pop += self.originPopData[origin]

        return(mf_list, mf_pop)


    def createGridData(self, sb, create_non_LD, cressman_power):

        # Create grids for each required source of data
        # Always create Local Dispersion grids
        #   optionally (create_non_LD = true) create non dispersed grids for testing/comparison
        # Optional cressman_power parameter (default 1) raises the weighting to the power of N

        if cressman_power != 1:
            logging.info('   Cressman weightings scaled to the power of {}\n'.format(cressman_power))

        rows = sb.projParams.background_rows
        cols = sb.projParams.background_cols
        gridCreator = GridCreate()

        logging.info('   Origins immobile     (.modelRun.grid_origins_immob) ...')
        if create_non_LD:
            self.grid_origins_immob = gridCreator.createGrid(rows, cols, sb.projParams.origin_data['XY'],self.originPopDataImmob)
        else:
            self.grid_origins_immob = None

        self.grid_origins_immob_LD = gridCreator.createGrid_LD(sb, sb.projParams.origin_data,  # bounds / locations
                                                               self.originPopDataImmob,        # the data to be spread
                                                               sb.projParams.originLocationIndex,  # location index
                                                               cressman_power)

        logging.info('\n   Origins remaining     (.modelRun.grid_origins_remain) ...')
        if create_non_LD:
            self.grid_origins_remain = gridCreator.createGrid(rows, cols, sb.projParams.origin_data['XY'],self.originPopData)
        else:
            self.grid_origins_remain = None

        self.grid_origins_remain_LD = gridCreator.createGrid_LD(sb, sb.projParams.origin_data,
                                                                self.originPopData, sb.projParams.originLocationIndex,
                                                                cressman_power)


        self.destination_file_grids = {}
        self.destination_file_grids_LD = {}
        for destination_file in self.destination_data_array:

            destination_data_partial = self.destination_data_array[destination_file]

            logging.info(f'\n   Populating Destination Location Index for {destination_file}...')
            destinationLocationIndex_partial = LocationIndex(sb.projParams, destination_data_partial) 
            
            logging.info(f'\n   Destinations onSite {destination_file}   (.modelRun.grid_dest_onSite)...')
            if create_non_LD:
                self.destination_file_grids[destination_file] = gridCreator.createGrid(rows, cols,
                                                            destination_data_partial['XY'], destination_data_partial['onSite'])
            else:
                self.grid_dest_onSite = None

            self.destination_file_grids_LD[destination_file] = gridCreator.createGrid_LD(sb, destination_data_partial,
                                                                    destination_data_partial['onSite'],
                                                                    destinationLocationIndex_partial,
                                                                    cressman_power)

        # create an index for quick access to Destination locations

        self.destination_data = self.concatenateDestinationDataArrays(self.destination_data_array)

        logging.info('\n   Populating Destination Location Index...')
        destinationLocationIndex = LocationIndex(sb.projParams, self.destination_data)

        logging.info('\n   Destinations onSite   (.modelRun.grid_dest_onSite)...')
        if create_non_LD:
            self.grid_dest_onSite = gridCreator.createGrid(rows, cols,
                                                    self.destination_data['XY'], self.destination_data['onSite'])
        else:
            self.grid_dest_onSite = None

        self.grid_dest_onSite_LD = gridCreator.createGrid_LD(sb, self.destination_data,
                                                             self.destination_data['onSite'],
                                                             destinationLocationIndex,
                                                             cressman_power)

        logging.info('\n   Destinations inTravel (.modelRun.grid_dest_inTravel)...')
        self.grid_dest_inTravel = gridCreator.createGrid_inTravel(sb, self.destination_data)


    def saveGridData(self, sb, file_prefix, save_destination_file_grids):

        # Save grid data to files, including option for both Locally Dispersed
        #   and not dispersed (probably never wanted, just for comparison)

        # create a grid file header
        header = 'ncols        {}\nnrows        {}\nxllcorner    {}\nyllcorner    {}\ncellsize     {}\nNODATA_value  -9999'.format(
            sb.projParams.aarea_cols, sb.projParams.aarea_rows,
            sb.projParams.aarea_bl_east, sb.projParams.aarea_bl_north,
            sb.projParams.aarea_csize)

        # clip grids to the Analysis area
        start_row = int((sb.projParams.background_rows - 1 - sb.projParams.aarea_rows) / 2)
        start_col = int((sb.projParams.background_cols - 1 - sb.projParams.aarea_cols) / 2)
        end_row = start_col + sb.projParams.aarea_rows
        end_col = start_col + sb.projParams.aarea_cols

        # use flipud (flip up down) as ASCII Grids are written top to bottom

        if self.grid_origins_immob is not None:
            self.grid_origins_immob = self.grid_origins_immob[start_row:end_row, start_col:end_col]
            filename = sb.projDir + file_prefix + 'origins_immob.asc'
            np.savetxt(filename, np.flipud(self.grid_origins_immob), fmt='%.4f', comments='', header=header)
            logging.info('   Written: ' + filename)

        if self.grid_origins_immob_LD is not None:
            self.grid_origins_immob_LD = self.grid_origins_immob_LD[start_row:end_row, start_col:end_col]
            filename = sb.projDir + file_prefix + 'origins_immob_LD.asc'
            np.savetxt(filename, np.flipud(self.grid_origins_immob_LD), fmt='%.4f', comments='', header=header)
            logging.info('   Written: ' + filename)

        if self.grid_origins_remain is not None:
            self.grid_origins_remain = self.grid_origins_remain[start_row:end_row, start_col:end_col]
            filename = sb.projDir + file_prefix + 'origins_remain.asc'
            np.savetxt(filename, np.flipud(self.grid_origins_remain), fmt='%.4f', comments='', header=header)
            logging.info('   Written: ' + filename)

        if self.grid_origins_remain_LD is not None:
            self.grid_origins_remain_LD = self.grid_origins_remain_LD[start_row:end_row, start_col:end_col]
            filename = sb.projDir + file_prefix + 'origins_remain_LD.asc'
            np.savetxt(filename, np.flipud(self.grid_origins_remain_LD), fmt='%.4f', comments='', header=header)
            logging.info('   Written: ' + filename)

        if self.grid_dest_inTravel is not None:
            self.grid_dest_inTravel = self.grid_dest_inTravel[start_row:end_row, start_col:end_col]
            filename = sb.projDir + file_prefix + 'dest_inTravel.asc'
            np.savetxt(filename, np.flipud(self.grid_dest_inTravel), fmt='%.4f', comments='', header=header)
            logging.info('   Written: ' + filename)

        if save_destination_file_grids:
            for dest_filename, destination_grid in self.destination_file_grids.items():

                destination_grid_filename = dest_filename.split(".")[0] # remove the csv part of the filename
                if destination_grid is not None:
                    destination_grid = destination_grid[start_row:end_row, start_col:end_col]
                    filename = sb.projDir + file_prefix + f'{destination_grid_filename}.asc'
                    np.savetxt(filename, np.flipud(destination_grid), fmt='%.4f', comments='', header=header)
                    logging.info('   Written: ' + filename)

            for dest_filename, destination_grid_LD in self.destination_file_grids_LD.items():

                destination_grid_filename = dest_filename.split(".")[0] # remove the csv part of the filename
                if destination_grid_LD is not None:
                    destination_grid_LD = destination_grid_LD[start_row:end_row, start_col:end_col]
                    filename = sb.projDir + file_prefix + f'{destination_grid_filename}_LD.asc'
                    np.savetxt(filename, np.flipud(destination_grid_LD), fmt='%.4f', comments='', header=header)
                    logging.info('   Written: ' + filename)


        if self.grid_dest_onSite is not None:
            self.grid_dest_onSite = self.grid_dest_onSite[start_row:end_row, start_col:end_col]
            filename = sb.projDir + file_prefix + 'dest_onSite.asc'
            np.savetxt(filename, np.flipud(self.grid_dest_onSite), fmt='%.4f', comments='', header=header)
            logging.info('   Written: ' + filename)

        if self.grid_dest_onSite_LD is not None:
            self.grid_dest_onSite_LD = self.grid_dest_onSite_LD[start_row:end_row, start_col:end_col]
            filename = sb.projDir + file_prefix + 'dest_onSite_LD.asc'
            np.savetxt(filename, np.flipud(self.grid_dest_onSite_LD), fmt='%.4f', comments='', header=header)
            logging.info('   Written: ' + filename)

        # save origin_all files (sum the component grids here)
        if self.grid_origins_remain is not None and self.grid_origins_immob is not None:
            self.grid_origins_all = self.grid_origins_remain + self.grid_origins_immob
            filename = sb.projDir + file_prefix + 'origins_all.asc'
            np.savetxt(filename, np.flipud(self.grid_origins_all), fmt='%.4f', comments='', header=header)
            logging.info('   Written: ' + filename)
        else:
            self.grid_origins_all = None

        if self.grid_origins_remain_LD is not None and self.grid_origins_immob_LD is not None:
            self.grid_origins_all_LD = self.grid_origins_remain_LD + self.grid_origins_immob_LD
            filename = sb.projDir + file_prefix + 'origins_all_LD.asc'
            np.savetxt(filename, np.flipud(self.grid_origins_all_LD), fmt='%.4f', comments='', header=header)
            logging.info('   Written: ' + filename)
        else:
            self.grid_origins_all_LD = None

        # save totals files (sum the component grids here)
        if self.grid_origins_remain is not None and self.grid_origins_immob is not None \
                and self.grid_dest_onSite is not None and self.grid_dest_inTravel is not None:
            self.grid_total = self.grid_origins_all + self.grid_dest_onSite + self.grid_dest_inTravel
            filename = sb.projDir + file_prefix + 'total.asc'
            np.savetxt(filename, np.flipud(self.grid_total), fmt='%.4f', comments='', header=header)
            logging.info('   Written: ' + filename)
        else:
            self.grid_total = None

        if self.grid_origins_remain_LD is not None and self.grid_origins_immob_LD is not None \
                and self.grid_dest_onSite_LD is not None and self.grid_dest_inTravel is not None:
            self.grid_total_LD = self.grid_origins_all_LD + self.grid_dest_onSite_LD + self.grid_dest_inTravel
            filename = sb.projDir + file_prefix + 'total_LD.asc'
            np.savetxt(filename, np.flipud(self.grid_total_LD), fmt='%.4f', comments='', header=header)
            logging.info('   Written: ' + filename)
        else:
            self.grid_total_LD = None

    def aggregateOutputData(self, sb):

        # Aggregate the data from a single model run into a set of running totals
        # these can be output at the end as a full set ot totals for all runs

        if self.grid_origins_immob is not None:
            try:
                sb.aggregate_grid_origins_immob += self.grid_origins_immob
            except:
                sb.aggregate_grid_origins_immob = np.copy(self.grid_origins_immob)

        if self.grid_origins_immob_LD is not None:
            try:
                sb.aggregate_grid_origins_immob_LD += self.grid_origins_immob_LD
            except:
                sb.aggregate_grid_origins_immob_LD = np.copy(self.grid_origins_immob_LD)

        if self.grid_origins_remain is not None:
            try:
                sb.aggregate_grid_origins_remain += self.grid_origins_remain
            except:
                sb.aggregate_grid_origins_remain = np.copy(self.grid_origins_remain)

        if self.grid_origins_remain_LD is not None:
            try:
                sb.aggregate_grid_origins_remain_LD += self.grid_origins_remain_LD
            except:
                sb.aggregate_grid_origins_remain_LD = np.copy(self.grid_origins_remain_LD)

        if self.grid_dest_inTravel is not None:
            try:
                sb.aggregate_grid_dest_inTravel += self.grid_dest_inTravel
            except:
                sb.aggregate_grid_dest_inTravel = np.copy(self.grid_dest_inTravel)

        if self.grid_dest_onSite is not None:
            try:
                sb.aggregate_grid_dest_onSite += self.grid_dest_onSite
            except:
                sb.aggregate_grid_dest_onSite = np.copy(self.grid_dest_onSite)

        if self.grid_dest_onSite_LD is not None:
            try:
                sb.aggregate_grid_dest_onSite_LD += self.grid_dest_onSite_LD
            except:
                sb.aggregate_grid_dest_onSite_LD = np.copy(self.grid_dest_onSite_LD)

        if self.grid_origins_all is not None:
            try:
                sb.aggregate_grid_origins_all += self.grid_origins_all
            except:
                sb.aggregate_grid_origins_all = np.copy(self.grid_origins_all)

        if self.grid_origins_all_LD is not None:
            try:
                sb.aggregate_grid_origins_all_LD += self.grid_origins_all_LD
            except:
                sb.aggregate_grid_origins_all_LD = np.copy(self.grid_origins_all_LD)

        if self.grid_total is not None:
            try:
                sb.aggregate_grid_total += self.grid_total
            except:
                sb.aggregate_grid_total = np.copy(self.grid_total)

        if self.grid_total_LD is not None:
            try:
                sb.aggregate_grid_total_LD += self.grid_total_LD
            except:
                sb.aggregate_grid_total_LD = np.copy(self.grid_total_LD)

    def saveCSVData(self, sb, file_prefix):

        # The non-Locally Dispersed data (probably never wanted)

        if self.grid_origins_remain is not None and self.grid_origins_immob is not None \
                and self.grid_dest_inTravel is not None and self.grid_dest_onSite is not None:

            filename = sb.projDir + file_prefix + 'results.csv'

            header = 'E, N, OriginRemain, OriginImmob, Intravel, OnSite, Total\n'
            halfcell = int(sb.projParams.background_csize / 2)

            try:
                with open(filename, 'w') as file_opened:

                    file_opened.write(header)

                    for x, y in [(x, y) for x in range(sb.projParams.aarea_cols) for y in range(sb.projParams.aarea_rows)]:
                        total = self.grid_origins_remain[y,x] \
                            + self.grid_origins_immob[y,x] \
                            + self.grid_dest_inTravel[y,x] \
                            + self.grid_dest_onSite[y,x]

                        if total > 0:
                            E = sb.projParams.aarea_bl_east + (sb.projParams.aarea_csize * x) + halfcell
                            N = sb.projParams.aarea_bl_north + (sb.projParams.aarea_csize * y) + halfcell

                            data = '{}, {}, {}, {}, {}, {}, {}\n'.format(E, N,
                                                                       round(self.grid_origins_remain[y,x],7),
                                                                       round(self.grid_origins_immob[y,x],7),
                                                                       round(self.grid_dest_inTravel[y,x],7),
                                                                       round(self.grid_dest_onSite[y,x],7),
                                                                       round(total,7))
                            file_opened.write(data)

                    logging.info('   Written: ' + filename)

            except IOError as e:
                logging.error(e)

        # the Locally Dispersed data

        if self.grid_origins_remain_LD is not None and self.grid_origins_immob_LD is not None \
                and self.grid_dest_inTravel is not None and self.grid_dest_onSite_LD is not None:

            filename = sb.projDir + file_prefix + 'results_LD.csv'

            header = 'E, N, OriginRemain, OriginImmob, Intravel, OnSite, Total\n'
            halfcell = int(sb.projParams.background_csize / 2)

            try:
                with open(filename, 'w') as file_opened:

                    file_opened.write(header)

                    for x, y in [(x, y) for x in range(sb.projParams.aarea_cols) for y in
                                 range(sb.projParams.aarea_rows)]:
                        total = self.grid_origins_remain_LD[y, x] \
                                + self.grid_origins_immob_LD[y, x] \
                                + self.grid_dest_inTravel[y, x] \
                                + self.grid_dest_onSite_LD[y, x]

                        if total > 0:
                            E = sb.projParams.aarea_bl_east + (sb.projParams.aarea_csize * x) + halfcell
                            N = sb.projParams.aarea_bl_north + (sb.projParams.aarea_csize * y) + halfcell

                            data = '{}, {}, {}, {}, {}, {}, {}\n'.format(E, N,
                                                                         round(self.grid_origins_remain_LD[y, x],7),
                                                                         round(self.grid_origins_immob_LD[y, x],7),
                                                                         round(self.grid_dest_inTravel[y, x],7),
                                                                         round(self.grid_dest_onSite_LD[y, x],7),
                                                                         round(total,7))
                            file_opened.write(data)

                    logging.info('   Written: ' + filename)

            except IOError as e:
                logging.error(e)
                
    def saveAggregatedOutputData(self, sb, file_prefix):

        # Save aggregated grid data to files, including option for both Locally Dispersed
        #   and not dispersed (probably never wanted, just for comparison)

        # create a grid file header
        header = 'ncols        {}\nnrows        {}\nxllcorner    {}\nyllcorner    {}\ncellsize     {}\nNODATA_value  -9999'.format(
            sb.projParams.aarea_cols, sb.projParams.aarea_rows,
            sb.projParams.aarea_bl_east, sb.projParams.aarea_bl_north,
            sb.projParams.aarea_csize)

        # aggregated grids are already clipped, so no cutting of the arrays needed

        # use flipud (flip up down) as ASCII Grids are written top to bottom

        try:
            sb.aggregate_grid_origins_immob
            filename = sb.projDir + file_prefix + 'origins_immob.asc'
            np.savetxt(filename, np.flipud(sb.aggregate_grid_origins_immob), fmt='%.4f', comments='', header=header)
            logging.info('   Written: ' + filename)
        except:
            pass

        try:
            sb.aggregate_grid_origins_immob_LD
            filename = sb.projDir + file_prefix + 'origins_immob_LD.asc'
            np.savetxt(filename, np.flipud(sb.aggregate_grid_origins_immob_LD), fmt='%.4f', comments='', header=header)
            logging.info('   Written: ' + filename)
        except:
            pass

        try:
            sb.aggregate_grid_origins_remain
            filename = sb.projDir + file_prefix + 'origins_remain.asc'
            np.savetxt(filename, np.flipud(sb.aggregate_grid_origins_remain), fmt='%.4f', comments='', header=header)
            logging.info('   Written: ' + filename)
        except:
            pass

        try:
            sb.aggregate_grid_origins_remain_LD
            filename = sb.projDir + file_prefix + 'origins_remain_LD.asc'
            np.savetxt(filename, np.flipud(sb.aggregate_grid_origins_remain_LD), fmt='%.4f', comments='', header=header)
            logging.info('   Written: ' + filename)
        except:
            pass

        try:
            sb.aggregate_grid_dest_inTravel
            filename = sb.projDir + file_prefix + 'dest_inTravel.asc'
            np.savetxt(filename, np.flipud(sb.aggregate_grid_dest_inTravel), fmt='%.4f', comments='', header=header)
            logging.info('   Written: ' + filename)
        except:
            pass

        try:
            sb.aggregate_grid_dest_onSite
            filename = sb.projDir + file_prefix + 'dest_onSite.asc'
            np.savetxt(filename, np.flipud(sb.aggregate_grid_dest_onSite), fmt='%.4f', comments='', header=header)
            logging.info('   Written: ' + filename)
        except:
            pass

        try:
            sb.aggregate_grid_dest_onSite_LD
            filename = sb.projDir + file_prefix + 'dest_onSite_LD.asc'
            np.savetxt(filename, np.flipud(sb.aggregate_grid_dest_onSite_LD), fmt='%.4f', comments='', header=header)
            logging.info('   Written: ' + filename)
        except:
            pass

        try:
            sb.aggregate_grid_origins_all
            filename = sb.projDir + file_prefix + 'origins_all.asc'
            np.savetxt(filename, np.flipud(sb.aggregate_grid_origins_all), fmt='%.4f', comments='', header=header)
            logging.info('   Written: ' + filename)
        except:
            pass

        try:
            sb.aggregate_grid_origins_all_LD
            filename = sb.projDir + file_prefix + 'origins_all_LD.asc'
            np.savetxt(filename, np.flipud(sb.aggregate_grid_origins_all_LD), fmt='%.4f', comments='', header=header)
            logging.info('   Written: ' + filename)
        except:
            pass

        try:
            sb.aggregate_grid_total
            filename = sb.projDir + file_prefix + 'total.asc'
            np.savetxt(filename, np.flipud(sb.aggregate_grid_total), fmt='%.4f', comments='', header=header)
            logging.info('   Written: ' + filename)
        except:
            pass

        try:
            sb.aggregate_grid_total_LD
            filename = sb.projDir + file_prefix + 'total_LD.asc'
            np.savetxt(filename, np.flipud(sb.aggregate_grid_total_LD), fmt='%.4f', comments='', header=header)
            logging.info('   Written: ' + filename)
        except:
            pass
