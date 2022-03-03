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
        self.originPopData = sb.projParams.origin_subgroups_pop[self.ageband].copy()
        #originPopData = np.array(sb.projParams.origin_subgroups_pop[self.ageband])

        originInitialPop = sum(self.originPopData)  # record initial total pop in this ageband
        destIncrease = 0  # (for checking) a record of how many dest pop transfers were made

        # create arrays to store all of the destination data values needed (and their grid indexes)
        self.dest_inTravel = []
        self.dest_onSite = []
        self.dest_XY = []

        for destdata in sb.projParams.destination_data:

            # grab the time profile for this collection

            time_profile = destdata['time_profile']
            logging.info('\n  destData: ' + time_profile + '\n')

            #     find the on site / in travel % for the specified time

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

            logging.info('    dest percentages for inTravel: '
                         + str(inTravel_pc) + '  onSite: ' + str(onSite_pc))

            # if we need to distribute unused origin pop to more distant WADs
            dest_inTravel_ratio = inTravel_pc / (inTravel_pc + onSite_pc)
            dest_onSite_ratio = onSite_pc / (inTravel_pc + onSite_pc)

            # loop through each destination row in the collection

            for dest in range(0, len(destdata['pop_data'])):

                if dest == DEST_DEBUG_LIMIT:  # fewer rows for debugging
                    break

                if dest % self.dest_sample_rate != 0:  # sample the dest
                    continue

                # grab the population for the specified age category

                dest_pop = destdata['subgroups_pop'][self.ageband][dest]

                # calculate the required OnSite Pop / InTravel Pop
                dest_inTravel_pop = dest_pop * inTravel_pc / 100
                dest_onSite_pop = dest_pop * onSite_pc / 100
                dest_req_pop = dest_inTravel_pop + dest_onSite_pop
                destIncrease += dest_req_pop

                # save values and grid index for our output grids
                self.dest_inTravel.append(dest_inTravel_pop)
                self.dest_onSite.append(dest_onSite_pop)
                self.dest_XY.append(destdata['XY'][dest])

                dest_E = destdata['eastings'][dest]
                dest_N = destdata['northings'][dest]

                logging.info('\n    Dest ' + str(dest)
                             + '. E: ' + str(dest_E) + ' N: ' + str(dest_N)
                             + ' Pop: ' + str(round(dest_pop,3))
                             + '  inTravel: ' + str(round(dest_inTravel_pop,3))
                             + '  onSite: ' + str(round(dest_onSite_pop,3))
                             + '  total req: ' + str(round(dest_req_pop,3)))

                dest_remove_check = 0

                # loop through each WAD pair (assume nearest is always first)

                dest_wad = destdata['WAD'][dest]

                for wad in dest_wad:
                    # add holders for extra data
                    wad.append(0)  # population count for the origin list
                    wad.append([])  # empty list for origins

                # find list of origins which are within the radius
                # (could we store the indexes of these for next model run?)

                for origin in range(0, len(sb.projParams.origin_eastings)):

                    if origin == ORIG_DEBUG_LIMIT:  # fewer rows for debugging
                        break

                    if origin % self.orig_sample_rate != 0:  # sample the dest
                        continue

                    loop_count += 1

                    orig_E = sb.projParams.origin_eastings[origin]
                    orig_N = sb.projParams.origin_northings[origin]

                    orig_pop = self.originPopData[origin]

                    # pythagoras gives us the distance between origin and destination
                    dist = math.sqrt((dest_E - orig_E) ** 2 + (dest_N - orig_N) ** 2)

                    logging.debug('      Orig ' + str(origin)
                                 + '. E: ' + str(orig_E) + ' N: ' + str(orig_N)
                                 + ' Pop: ' + str(round(orig_pop, 2))
                                 + ' distance: ' + str(round(dist, 2)))

                    # which wad is this distance relevant to
                    for wad in dest_wad:
                        if dist <= wad[0] or wad[0] == 0:  # within range or final zero catch all
                            wad[2] += orig_pop  # store the total origin population
                            wad[3].append(origin)  # add the origin index to our list
                            break

                # The dest wad is now fully populated with origin indexes and origin pop total

                residue_inTravel = 0   # keep track of unmet wad pop requirement, to pass up to next radius
                residue_onSite = 0

                available_pop = 0    # keep track of wad pop available (prev and current WADs)
                available_origins = []  # and list of those origins who will provide it

                for wad in dest_wad:

                    origin_wad_pop = wad[2]
                    available_pop += origin_wad_pop

                    # figure out how many people need pulling for this WAD
                    rad = wad[0]
                    pc = wad[1]
                    wad_inTravel = dest_inTravel_pop * pc / 100
                    wad_onSite = dest_onSite_pop * pc / 100
                    wad_total = wad_inTravel + wad_onSite  # to remove from these origins

                    # add to our current balance of required pop for dest
                    residue_inTravel += wad_inTravel
                    residue_onSite += wad_onSite
                    residue_total = residue_inTravel + residue_onSite

                    logging.debug('      ' + str(pc) + '%  within ' + str(rad) + 'm -> '
                                 + ' inTravel ' + str(round(wad_inTravel, 3))
                                 + ', onSite ' + str(round(wad_onSite, 3))
                                 + ' total ' + str(round(wad_total,3))
                                 + ' from origins (' + str(round(origin_wad_pop,3)) + ' available)')
                    logging.debug('              residue: '
                                 + ' inTravel ' + str(round(residue_inTravel, 3))
                                 + ', onSite ' + str(round(residue_onSite, 3))
                                 + ' total ' + str(round(residue_total,3)))
                    orig_remove_check = 0

                    if available_pop > 0:  # some origin population is available to take

                        available_origins.extend(wad[3]) # add these origins to the available list

                        if available_pop > residue_total:
                            # enough to satisfy requirement fully, satisfy all residue
                            wad_remove_total = residue_total
                            wad_remove_inTravel = residue_inTravel
                            wad_remove_onSite = residue_onSite
                            residue_inTravel = 0
                            residue_onSite = 0
                        else:
                            # not enough origin pop in this WAD, use it all up and update residues
                            wad_remove_total = available_pop
                            wad_remove_inTravel = wad_remove_total * dest_inTravel_ratio
                            wad_remove_onSite = wad_remove_total * dest_onSite_ratio
                            residue_inTravel -= wad_remove_inTravel
                            residue_onSite -= wad_remove_onSite
                            logging.debug('not enough origin pop in this WAD!')

                        # let's do some removing of pops from origins
                        for origin in available_origins:
                            # go through each origin index, remove in proportion with origin pop
                            pop = self.originPopData[origin]
                            orig_remove_total = pop / available_pop * wad_remove_total
                            orig_remove_inTravel = pop / available_pop * wad_remove_inTravel
                            orig_remove_onSite = pop / available_pop * wad_remove_onSite
                            self.originPopData[origin] -= orig_remove_total
                            logging.debug('        removed '
                                         + str(round(orig_remove_total,3))
                                         + ' from ' + str(origin) + ' (' + str(round(pop,3)) + ')')
                            orig_remove_check += orig_remove_total

                        available_pop -= wad_remove_total  # update total pop available

                        if available_pop == 0.0:  # may need some rounding, or v small val comparison?
                            available_origins = []  # all used up, continue with empty array of origins

                        logging.debug('        Orig remove check: ' + str(round(orig_remove_check,3)))
                        dest_remove_check += orig_remove_check

                # destination is complete, check we removed the full amount
                if round(dest_req_pop,3) == round(dest_remove_check,3):
                    logging.info('      Dest remove check SUCCESS: ' + str(round(dest_remove_check,3)))
                else:
                    logging.info('      Dest remove check FAIL: ' + str(round(dest_remove_check, 3)))

        originFinalPop = sum(self.originPopData)  # record initial total pop in this ageband

        logging.info('\n  Run Complete - Loop count: ' + str(loop_count)
                     + ' in ' + str(round(time.time() - initialTime,1)) + ' seconds')

        logging.info('\n  Origin pop initial / final / diff: '
                     + str(round(originInitialPop,3)) + ' / ' + str(round(originFinalPop,3))
                     + ' / ' + str(round(originInitialPop - originFinalPop,3))
                     + '\n  Dest pop requested: ' + str(round(destIncrease,3)))

        # raise ValueError('sorry, not good')

    def addMins(self, tm, mins):
        fulldate = datetime.datetime(100, 1, 1, tm.hour, tm.minute, 0)
        fulldate = fulldate + datetime.timedelta(minutes=mins)
        return fulldate.time()

    def createGridData(self, sb):

        # create grids for each required source of data

        rows = sb.projParams.background_rows
        cols = sb.projParams.background_cols

        logging.info('   Origins remaining     (.modelRun.grid_origins) ...')
        self.grid_origins = self.createGrid(rows, cols, sb.projParams.origin_XY,self.originPopData)

        logging.info('   Destinations inTravel (.modelRun.grid_dest_inTravel)...')
        self.grid_dest_inTravel = self.createGrid(rows, cols, self.dest_XY, self.dest_inTravel)

        logging.info('   Destinations onSite   (.modelRun.grid_dest_onSite)...')
        self.grid_dest_onSite = self.createGrid(rows, cols, self.dest_XY, self.dest_onSite)

    def saveGridData(self, sb, file_prefix):

        # save grid data to files

        # create a grid file header
        header = ''
        for (param, value) in sb.projParams.background_header.items():
            header = header + param + ' ' + str(value) + '\n'

        # use flipud (flip up down) as ASCII Grids are written top to bottom

        filename = sb.projDir + file_prefix + 'origins.asc'
        np.savetxt(filename, np.flipud(self.grid_origins), fmt='%.4f', comments='', header=header)
        logging.info('   Written: ' + filename)

        filename = sb.projDir + file_prefix + 'dest_inTravel.asc'
        np.savetxt(filename, np.flipud(self.grid_dest_inTravel), fmt='%.4f', comments='', header=header)
        logging.info('   Written: ' + filename)

        filename = sb.projDir + file_prefix + 'dest_onSite.asc'
        np.savetxt(filename, np.flipud(self.grid_dest_onSite), fmt='%.4f', comments='', header=header)
        logging.info('   Written: ' + filename)

        logging.info('\n  Model Data Saved.')

    def createGrid(self, rows, cols, XY_array, vals_array):

        grid = np.zeros((rows,cols))

        minX = 10  # just for checking, we remove/comment out later
        maxX = 10
        minY = 10
        maxY = 10

        for row in range(0, len(XY_array)):
            (X, Y) = XY_array[row]

            if X < minX:
                minX = X
            if X > maxX:
                maxX = X
            if Y < minY:
                minY = Y
            if Y > maxY:
                maxY = Y

            if X >= rows or Y >= cols:
                logging.info('     Ignoring out of bounds value at row '
                             + str(row) + ' ('+ str(X)+','+str(Y)+')')
            else:
                val = vals_array[row]
                grid[Y,X] = val

        return grid