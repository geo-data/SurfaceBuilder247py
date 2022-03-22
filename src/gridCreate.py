#
# Python version of Surface Builder 24/7
#
# Jan 2022
# GeoData Institute
# University of Southampton
# on behalf of ONS

# Import core modules

import logging
import time
import math
import copy

# Import additional modules (which will need to be installed)

import numpy as np

from locationIndex import LocationIndex


# A class for creating Grids from point based population data

class GridCreate:

    def createGrid(self, rows, cols, XY_array, vals_array):

        # create a grid just using XY locations and values, no redistribution of any kind

        grid = np.zeros((rows,cols))

        for row in range(len(XY_array)):
            (X, Y) = XY_array[row]

            if X >= cols or Y >= rows:
                logging.info('     Ignoring out of bounds value at row '
                             + str(row) + ' ('+ str(X)+','+str(Y)+')')
            else:
                val = vals_array[row]
                grid[Y,X] += val

        return grid

    def createGrid_inTravel(self, sb, destination_data):

        # for each dest array in travel value
        #    for each wad (rad:pc)
        #       allocate background cells (within radius) and total bg amount into each wad
        #
        #    for each wad (pc)
        #       -> amount within this radius
        #            for each contained background cell
        #                add to same place in new grid dest pop * wad pc * background weighting / total weighting

        rows = sb.projParams.background_rows
        cols = sb.projParams.background_cols

        grid = np.zeros((rows,cols))

        loop_count = 0
        lost_pop = 0
        initialTime = time.time()

        logging.info('\n     Populating Background Location Index...')
        backgroundLocationIndex = LocationIndex(sb.projParams, sb.projParams.background_data)

        for row in range(len(destination_data['inTravel'])):
            E = destination_data['eastings'][row]
            N = destination_data['northings'][row]
            dest_WAD = copy.deepcopy(destination_data['WAD'][row])
            # we use deepcopy here to start with a full copy of the wad, with zero pop count and empty location list
            dest_pop = destination_data['inTravel'][row]
            dest_pop_leftover = 0  # for storing pop amounts that have nowhere to go within a WAD

            if row % 1000 == 0:
                #    print('.', end='', flush=True)
                logging.info('     {} Destinations'.format(row))

            largest_radius = math.sqrt(dest_WAD[len(dest_WAD) - 1][0])
            potential_background_vals_array = backgroundLocationIndex.possible_locations(E, N, largest_radius)

            for potential_background_vals in potential_background_vals_array:
                for bgcell in potential_background_vals:

                    bg_E = sb.projParams.background_data['eastings'][bgcell]
                    bg_N = sb.projParams.background_data['northings'][bgcell]
                    bg_val = sb.projParams.background_data['bgval'][bgcell]

                    # pythagoras gives us the distance between origin and destination
                    dist_sq = (bg_E - E) ** 2 + (bg_N - N) ** 2

                    # which wad is this distance relevant to
                    for wad in dest_WAD:
                        if dist_sq <= wad[0]:  # within range
                            if wad[1] > 0:  # any data (pc > 0) to be held in here at all
                                wad[2] += bg_val  # store the total origin population
                                wad[3].append(bgcell)  # add the background index to our list
                                break
                            # otherwise it will get added to the next wad outward
                        loop_count += 1

            # the dest wad is now fully populated with grid cell indexes and total amounts
            # loop through it again, spreading the dest inTravel pop into the relevant grid cells

            for wad in dest_WAD:
                pc = wad[1]

                if pc > 0:  # any data in here?
                    bg_tot = wad[2]  # available background values
                    bg_cells = wad[3]  # available background cells

                    if bg_tot > 0:
                        pop = dest_pop_leftover + (dest_pop * pc / 100)  # pop to disperse amongst these cells
                        dest_pop_leftover = 0
                        mult = pop / bg_tot

                        for bg in bg_cells:  # loop through the background cell indexes
                            (X, Y) = sb.projParams.background_data['XY'][bg]
                            bgval = sb.projParams.background_data['bgval'][bg]
                            loop_count += 1
                            # add the share of the pop (dest * wad pc) to the relevant grid cell (grid amount / grid tot)
                            # Row - Y, Col - X
                            #if X >= cols or Y >= rows:
                            #    logging.info('     Ignoring out of bounds value at row '
                            #                 + str(loop_count) + ' (' + str(X) + ',' + str(Y) + ')')
                            #else:
                            grid[Y, X] += bgval * mult

                    else:
                        # there are no non-zero bg cells in this WAD, but we still have pop to disperse
                        dest_pop_leftover += dest_pop * pc / 100

            if dest_pop_leftover > 0:
                # we've been through all of our WADs and there is still undistributed population (hopefully unlikely)
                logging.info('     Destination {} had {:.3f} unallocated population'.format(row, dest_pop_leftover))
                lost_pop += dest_pop_leftover

        logging.info('\n     created - Loop count: {} in {} seconds'.format(loop_count,
                                                                            round(time.time() - initialTime, 1)))
        logging.info('     Values total: {:.3f}  Dispersed Grid total: {:.3f}  Unallocated: {:.3f}'.format(sum(destination_data['inTravel']),
                                                                                      sum(sum(grid)), lost_pop))

        return grid

    def createGrid_LD(self, sb, locations, vals_array, locationIndex, cressman_power):

        # create a grid using Local Dispersion

        halfcell = int(sb.projParams.background_csize / 2)
        csize = sb.projParams.background_csize
        rows = sb.projParams.background_rows
        cols = sb.projParams.background_cols
        grid = np.zeros((rows, cols))

        # For each location:

        for loc in range(len(vals_array)):

            popval = vals_array[loc]
            if popval == 0:
                continue  # nothing here to disperse!

            #  1. Create a list of other (populated) data points within the LDParam radius
            #   [ (E, N, dist) ... ]

            loc_E = locations['eastings'][loc]
            loc_N = locations['northings'][loc]

            LD = locations['LD'][loc]
            possible_neighbours_array = locationIndex.possible_locations(loc_E, loc_N, LD)
            neighbours = []
            for possible_neighbours in possible_neighbours_array:
                for poss in possible_neighbours:
                    if poss != loc:  # not itself!
                        poss_E = locations['eastings'][poss]
                        poss_N = locations['northings'][poss]
                        dist = math.sqrt((loc_E - poss_E) ** 2 + (loc_N - poss_N) ** 2)
                        if dist > 0 and dist <= LD:
                            neighbours.append((poss_E, poss_N, dist))  # add a tuple to the list

            #  2. Calculate AVI (Inter Centroid Distance)
            #       if none -> LDradius, if 1 -> distance, otherwise average distance

            if len(neighbours) == 0:
                avi = LD
            elif len(neighbours) == 1:
                avi = neighbours[0][2]  # use the single neighbour's distance
            else:
                avi = sum([elem[2] for elem in neighbours]) / len(neighbours)

            #  3. Create a list of cells whose centroids fall within the AVI
            #    insert weight into each item
            #   [ (X, Y, Ecen, Ncen, W) .... ]

            # 3a. define a bounding box to search for cell centroids potentially within range
            # add or subtract 1m to ensure we include (potential) exact centroid matches
            X1 = int(((loc_E - avi + halfcell - 1) - sb.projParams.background_bl_east) / csize)
            Y1 = int(((loc_N - avi + halfcell - 1) - sb.projParams.background_bl_north) / csize)
            if X1 < 0:
                X1 = 0
            if Y1 < 0:
                Y1 = 0
            X2 = int(((loc_E + avi - halfcell + 1) - sb.projParams.background_bl_east) / csize)
            Y2 = int(((loc_N + avi - halfcell + 1) - sb.projParams.background_bl_north) / csize)
            if X2 >= cols:
                X2 = cols - 1
            if Y2 >= rows:
                Y2 = rows - 1

            # 3b Loop through the bounding box, check exact distance to each cell centroid (Xc, Yc)
            cell_list = []
            avi_sq = avi ** 2
            total_weight = 0
            for x, y in [(x, y) for x in range(X1, X2 + 1) for y in range(Y1, Y2 + 1)]:
                Xc = sb.projParams.background_bl_east + (x * csize) + halfcell
                Yc = sb.projParams.background_bl_north + (y * csize) + halfcell
                dist = math.sqrt((loc_E - Xc) ** 2 + (loc_N - Yc) ** 2)
                if dist <= avi:
                    dist_sq = dist ** 2
                    weight = ((avi_sq - dist_sq) / (avi_sq + dist_sq))
                    if weight > 0:  # only add to our list of cells for dispersion if weight is positive
                        weight = weight ** cressman_power
                        total_weight += weight
                        # a tuple of the cell location and weight
                        cell_list.append((x, y, weight))

            #  4. Populate the results grid
            #    each affected cell from 3 receives a portion of the location's population
            #      X,Y += pop * weight(x,y)
            if total_weight > 0:
                mult = popval / total_weight
                for (x, y, weight) in cell_list:
                    grid[y, x] += weight * mult
            else:
                # no centroids are within the AVI, so put all of the pop in its own cell
                x = int((loc_E - sb.projParams.background_bl_east) / csize)
                y = int((loc_N - sb.projParams.background_bl_north) / csize)
                grid[y, x] += popval

        # all locations have been processed, grid is fully populated

        logging.info('     Values total: {:.3f}  Dispersed Grid total: {:.3f}'.format(sum(vals_array), sum(sum(grid))))

        return grid

""" 
# Unused attempt at a more efficient inTravel grid
# turned out 3 times slower than the current method

    def createGrid_inTravel2(self, sb, destination_data):

        # for each destination
        #  loop through BB of largest wad
        #    grab X,Y of BG grid, bg_val
        #      if bg_val > 0
        #        dist_sq = X,Y to dest(x,y)
        #    loop through wads
        #      add (x,y) to wad list, increment wad total
        #
        #  loop though wads
        #    dest bg total for wad
        #      for each x,y in list for wad
        #        inTrav(x,y) += inTrav * bg[x,y] / bg_total

        halfcell = int(sb.projParams.background_csize / 2)
        csize = sb.projParams.background_csize
        rows = sb.projParams.background_rows
        cols = sb.projParams.background_cols

        grid = np.zeros((rows,cols))

        loop_count = 0
        lost_pop = 0
        initialTime = time.time()

        for row in range(len(destination_data['inTravel'])):
            loc_E = destination_data['eastings'][row]
            loc_N = destination_data['northings'][row]
            dest_WAD = copy.deepcopy(destination_data['WAD'][row])
            # we use deepcopy here to start with a full copy of the wad, with zero pop count and empty location list
            dest_pop = destination_data['inTravel'][row]
            dest_pop_leftover = 0  # for storing pop amounts that have nowhere to go within a WAD

            if row % 1000 == 0:
                #    print('.', end='', flush=True)
                logging.info('     {} Destinations'.format(row))

            # largest_radius_sq = dest_WAD[len(dest_WAD) - 1][0]
            largest_radius = math.sqrt(dest_WAD[len(dest_WAD) - 1][0])

            # find the bounding box for this wad
            X1 = int(((loc_E - largest_radius) - sb.projParams.background_bl_east) / csize)
            Y1 = int(((loc_N - largest_radius) - sb.projParams.background_bl_north) / csize)
            if X1 < 0:
                X1 = 0
            if Y1 < 0:
                Y1 = 0
            X2 = int(((loc_E + largest_radius) - sb.projParams.background_bl_east) / csize)
            Y2 = int(((loc_N + largest_radius) - sb.projParams.background_bl_north) / csize)
            if X2 >= cols:
                X2 = cols - 1
            if Y2 >= rows:
                Y2 = rows - 1

            # Loop through the bounding box, check exact distance to each cell centroid (Xc, Yc)
            for x, y in [(x, y) for x in range(X1, X2) for y in range(Y1, Y2)]:

                bg_val = sb.projParams.background_array[y, x]

                if bg_val > 0.0001:  # threshold
                    Xc = sb.projParams.background_bl_east + (x * csize) + halfcell
                    Yc = sb.projParams.background_bl_north + (y * csize) + halfcell
                    dist_sq = (loc_E - Xc) ** 2 + (loc_N - Yc) ** 2

                    # which wad is this distance relevant to
                    for wad in dest_WAD:
                        if dist_sq <= wad[0]:  # within range
                            if wad[1] > 0:  # any data (pc > 0) to be held in here at all
                                wad[2] += bg_val  # store the total origin population
                                wad[3].append((x, y, bg_val))  # add the background index to our list
                                break
                            # otherwise it will get added to the next wad outward
                        loop_count += 1

            # wad is fully populated

            for wad in dest_WAD:
                pc = wad[1]

                if pc > 0:  # any data in here?
                    bg_tot = wad[2]  # available background values
                    bg_cells = wad[3]  # available background cells

                    if bg_tot > 0:
                        pop = dest_pop_leftover + (dest_pop * pc / 100)  # pop to disperse amongst these cells
                        dest_pop_leftover = 0
                        mult = pop / bg_tot

                        for (X, Y, bgval) in bg_cells:  # loop through the background cell indexes

                            loop_count += 1
                            grid[Y, X] += bgval * mult

                    else:
                        # there are no non-zero bg cells in this WAD, but we still have pop to disperse
                        dest_pop_leftover += dest_pop * pc / 100

            if dest_pop_leftover > 0:
                # we've been through all of our WADs and there is still undistributed population (hopefully unlikely)
                logging.info('     Destination {} had {:.3f} unallocated population'.format(row, dest_pop_leftover))
                lost_pop += dest_pop_leftover

        logging.info('\n     created - Loop count: {} in {} seconds'.format(loop_count,
                                                                            round(time.time() - initialTime, 1)))
        logging.info('     Values total: {:.3f}  Dispersed Grid total: {:.3f}  Unallocated: {:.3f}'.format(
            sum(destination_data['inTravel']),
            sum(sum(grid)), lost_pop))

        return grid
"""