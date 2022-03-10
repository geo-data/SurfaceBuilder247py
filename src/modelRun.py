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
        self.originPopData = sb.projParams.origin_data['subgroups_pop'][self.ageband].copy()

        # transfer immobile population
        self.originPopDataImmob = []
        for origin in range(0, len(self.originPopData)):
            immob = self.originPopData[origin] * sb.projParams.origin_data['subgroups_mob'][self.ageband][origin]
            self.originPopDataImmob.append(immob)
            self.originPopData[origin] -= immob

        logging.info(  '\n  Immobile population removed: {}'.format(round(sum(self.originPopDataImmob),3)) )

        originInitialPop = sum(self.originPopData)  # record initial total pop in this ageband
        destIncrease = 0  # (for checking) a record of how many dest pop transfers were made

        # create arrays to store all of the destination data values needed (and their grid indexes)
        self.dest_inTravel = []
        self.dest_onSite = []
        self.dest_EN = []  # list of Easting, Northing tuples
        self.dest_XY = []
        self.dest_WAD = []

        # only need to make the full list once
        all_origins = range(0, len(sb.projParams.origin_data['eastings']))

        for destdata in sb.projParams.destination_data:

            logging.info('\n  New destination collection (' + destdata['Filename'] + ') ...')

            current_time_profile = ''  # only recalculate the inTravel and onSite percentages when things change

            # loop through each destination row in the collection

            for dest in range(0, len(destdata['pop_data'])):

                if dest == DEST_DEBUG_LIMIT:  # fewer rows for debugging
                    break

                if dest % self.dest_sample_rate != 0:  # sample the dest
                    continue

                # grab the time profile and calculate the percentages (if things change)
                time_profile = destdata['time_profiles'][dest]

                if time_profile != current_time_profile:
                    logging.info('\n  New time profile: ' + time_profile + '\n')

                    (inTravel_pc, onSite_pc) = self.timeProfileLookup(sb, time_profile)

                    logging.info('    dest percentages for inTravel: '
                                 + str(inTravel_pc) + '  onSite: ' + str(onSite_pc))

                    # if we need to distribute unused origin pop to more distant WADs
                    dest_inTravel_ratio = inTravel_pc / (inTravel_pc + onSite_pc)
                    dest_onSite_ratio = onSite_pc / (inTravel_pc + onSite_pc)

                    current_time_profile = time_profile  # remember for next time

                # grab the population for the specified age category

                dest_pop = destdata['subgroups_pop'][self.ageband][dest]

                # calculate the required OnSite Pop / InTravel Pop
                dest_inTravel_pop = dest_pop * inTravel_pc / 100
                dest_onSite_pop = dest_pop * onSite_pc / 100
                dest_req_pop = dest_inTravel_pop + dest_onSite_pop
                destIncrease += dest_req_pop

                dest_E = destdata['eastings'][dest]
                dest_N = destdata['northings'][dest]

                # save values and grid index for our output grids
                self.dest_inTravel.append(dest_inTravel_pop)
                self.dest_onSite.append(dest_onSite_pop)
                self.dest_EN.append((dest_E, dest_N))
                self.dest_XY.append(destdata['XY'][dest])
                self.dest_WAD.append(destdata['WAD'][dest])

                logging.info('\n    Dest ' + str(dest)
                             + '. E: ' + str(dest_E) + ' N: ' + str(dest_N)
                             + ' Pop: ' + str(round(dest_pop,3))
                             + '  inTravel: ' + str(round(dest_inTravel_pop,3))
                             + '  onSite: ' + str(round(dest_onSite_pop,3))
                             + '  total req: ' + str(round(dest_req_pop,3)))

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

                dest_wad = destdata['WAD'][dest]

                for wad in dest_wad:
                    # reset the holders for extra data
                    wad[2] = 0  # population count for the origin list
                    wad[3] = []  # empty list for origins

                # find list of origins which are within the radius

                # index method:
                largest_radius = dest_wad[len(dest_wad)-1][0]
                if largest_radius == 0:
                    largest_radius = dest_wad[len(dest_wad) - 2][0]
                potential_origins = sb.projParams.originLocationIndex.possible_locations(dest_E, dest_N, largest_radius)

                #for origin in all_origins:
                for origin in potential_origins:

                    if origin == ORIG_DEBUG_LIMIT:  # fewer rows for debugging
                        break

                    if origin % self.orig_sample_rate != 0:  # sample the origins
                        continue

                    if origin in mf_origins_dest:  # we have Major Flowed this origin already
                        continue

                    loop_count += 1

                    orig_E = sb.projParams.origin_data['eastings'][origin]
                    orig_N = sb.projParams.origin_data['northings'][origin]

                    orig_pop = self.originPopData[origin]

                    # TODO - if the pop is zero, though unlikely,we might as well stop (continue) here?

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
                            # breakdowns not used, probably not needed
                            #orig_remove_inTravel = pop / available_pop * wad_remove_inTravel
                            #orig_remove_onSite = pop / available_pop * wad_remove_onSite
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

        for origin in range(0, len(sb.projParams.origin_data['eastings'])):
            if sb.projParams.origin_data['ID'][origin][:mf_len] == mf_ID:
                mf_list.append(origin)
                mf_pop += self.originPopData[origin]

        return(mf_list, mf_pop)


    def createGridData(self, sb):

        # create grids for each required source of data

        rows = sb.projParams.background_rows
        cols = sb.projParams.background_cols

        logging.info('   Origins remaining     (.modelRun.grid_origins) ...')
        self.grid_origins = self.createGrid(rows, cols, sb.projParams.origin_data['XY'],self.originPopData)

        logging.info('\n   Destinations inTravel (.modelRun.grid_dest_inTravel)...')
        self.grid_dest_inTravel = self.createGrid_inTravel(sb)

        logging.info('\n   Destinations onSite   (.modelRun.grid_dest_onSite)...')
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

            if X >= cols or Y >= rows:
                logging.info('     Ignoring out of bounds value at row '
                             + str(row) + ' ('+ str(X)+','+str(Y)+')')
            else:
                val = vals_array[row]
                grid[Y,X] = val

        return grid

    def createGrid_inTravel(self, sb):

        # for each dest array in travel value
        #    for each wad (rad:pc)
        #       allocate background cells (within radius) and total bg amount into each wad
        #
        #    for each wad (pc)
        #       -> amount within this radius
        #            for each contained background cell
        #                add to same place in new grid dest pop * wad pc * background weighting /  total weighting

        rows = sb.projParams.background_rows
        cols = sb.projParams.background_cols

        grid = np.zeros((rows,cols))

        loop_count = 0
        initialTime = time.time()

        for row in range(0, len(self.dest_inTravel)):
            (E, N) = self.dest_EN[row]
            dest_WAD = self.dest_WAD[row]
            dest_pop = self.dest_inTravel[row]

            for wad in dest_WAD:
                if wad[1] > 0:  # any data (pc > 0) to be held in here at all
                    # reset holders for extra data
                    wad[2] = 0     # background grid total
                    wad[3] = []   # empty list for background grid cell indexes

            bgindex = 0

            for bgcell in sb.projParams.background_values:

                    bg_E = bgcell[2]
                    bg_N = bgcell[3]
                    bg_val = bgcell [4]

                    # pythagoras gives us the distance between origin and destination
                    dist = math.sqrt((bg_E - E) ** 2 + (bg_N - N) ** 2)

                    # which wad is this distance relevant to
                    for wad in dest_WAD:
                        if dist <= wad[0] or wad[0] == 0:  # within range or final zero catch all
                            if wad[1] > 0:  # any data (pc > 0) to be held in here at all
                                wad[2] += bg_val  # store the total origin population
                                wad[3].append(bgindex)  # add the origin index to our list
                                break
                            # otherwise it will get added to the next wad outward
                        loop_count += 1

                    bgindex += 1

            # the dest wad is now fully populated with grid cell indexes and total amounts
            # loop through it again, spreading the dest inTravel pop into the relevant grid cells

            bg_tot = 0  # start with no background data
            bg_cells = []  # and an empty list of background cell indexes

            for wad in dest_WAD:
                pc = wad[1]

                if pc > 0:  # any data in here?
                    bg_tot += wad[2]  # accumulate the background values
                    bg_cells.extend(wad[3])  # add the background cells to the list to spread population

                    for bg in bg_cells:  # loop through the background cell indexes
                        bgcell = sb.projParams.background_values[bg]
                        X = bgcell[0]
                        Y = bgcell[1]
                        val = bgcell[4]
                        # add the share of the pop (dest * wad pc) to the relevant grid cell (grid amount / grid tot)
                        # Row - Y, Col - X

                        if X >= cols or Y >= rows:
                            logging.info('     Ignoring out of bounds value at row '
                                         + str(loop_count) + ' (' + str(X) + ',' + str(Y) + ')')
                        else:
                            grid[Y, X] = val

                        grid[Y,X] += dest_pop * pc / 100 * val / bg_tot
                        loop_count += 1

        logging.info('\n     created - Loop count: ' + str(loop_count)
                     + ' in ' + str(round(time.time() - initialTime,1)) + ' seconds')

        return grid