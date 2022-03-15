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

# Import additional modules (which will need to be installed)

import numpy as np

from locationIndex import LocationIndex


# A class for creating Grids from point based population data

class GridCreate:

    def createGrid(self, rows, cols, XY_array, vals_array):

        # create a grid just using XY locations and values, no redistribution of any kind

        grid = np.zeros((rows,cols))

        minX = 10  # just for checking, we remove/comment out later
        maxX = 10
        minY = 10
        maxY = 10

        for row in range(len(XY_array)):
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
        #                add to same place in new grid dest pop * wad pc * background weighting /  total weighting

        rows = sb.projParams.background_rows
        cols = sb.projParams.background_cols

        grid = np.zeros((rows,cols))

        loop_count = 0
        lost_pop = 0
        initialTime = time.time() 

        # determine background data selection method, depending on if there is a large study area
        # if wider/higher than twice the common maximum WAD radius it will be highly beneficial
        if (sb.projParams.sarea_tr_east - sb.projParams.sarea_bl_east) >= 120000 \
                or (sb.projParams.sarea_tr_north - sb.projParams.sarea_bl_north) >= 120000:
            logging.info('\n     Large Study Area: Using Background Location Index')
            use_location_index = True
            # create an index for quick access to non-zero Background locations
            backgroundLocationIndex = LocationIndex(sb.projParams, sb.projParams.background_data)
        else:
            use_location_index = False
            all_background_vals = range(len(sb.projParams.background_data['eastings']))

        for row in range(len(destination_data['inTravel'])):
            E = destination_data['eastings'][row]
            N = destination_data['northings'][row]
            dest_WAD = destination_data['WAD'][row]
            dest_pop = destination_data['inTravel'][row]
            dest_pop_leftover = 0  # for storing pop amounts that have nowhere to go within a WAD

            for wad in dest_WAD:
                if wad[1] > 0:  # any data (pc > 0) to be held in here at all
                    # reset holders for extra data
                    wad[2] = 0     # background grid total
                    wad[3] = []   # empty list for background grid cell indexes

            # index method:
            if use_location_index:
                largest_radius = dest_WAD[len(dest_WAD) - 1][0]
                if largest_radius == 0:
                    largest_radius = dest_WAD[len(dest_WAD) - 2][0]
                potential_background_vals = backgroundLocationIndex.possible_locations(E, N, largest_radius)
            else:
                potential_background_vals = all_background_vals

            for bgcell in potential_background_vals:

                    bg_E = sb.projParams.background_data['eastings'][bgcell]
                    bg_N = sb.projParams.background_data['northings'][bgcell]
                    bg_val = sb.projParams.background_data['bgval'][bgcell]

                    # pythagoras gives us the distance between origin and destination
                    dist = math.sqrt((bg_E - E) ** 2 + (bg_N - N) ** 2)

                    # which wad is this distance relevant to
                    for wad in dest_WAD:
                        if dist <= wad[0] or wad[0] == 0:  # within range or final zero catch all
                            if wad[1] > 0:  # any data (pc > 0) to be held in here at all
                                wad[2] += bg_val  # store the total origin population
                                wad[3].append(bgcell)  # add the background index to our list
                                break
                            # otherwise it will get added to the next wad outward
                        loop_count += 1

            # the dest wad is now fully populated with grid cell indexes and total amounts
            # loop through it again, spreading the dest inTravel pop into the relevant grid cells

            bg_tot = 0  # start with no background data
            bg_cells = []  # and an empty list of background cell indexes

            for wad in dest_WAD:
                pc = wad[1]

                if pc > 0:  # any data in here?
                    bg_tot += wad[2]  # accumulate the background values
                    bg_cells.extend(wad[3])  # add the background cells to the list to spread population

                    if bg_tot > 0:
                        pop = dest_pop_leftover + (dest_pop * pc / 100)  # pop to disperse amongst these cells
                        dest_pop_leftover = 0

                        for bg in bg_cells:  # loop through the background cell indexes
                            (X, Y) = sb.projParams.background_data['XY'][bg]
                            bgval = sb.projParams.background_data['bgval'][bg]
                            loop_count += 1
                            # add the share of the pop (dest * wad pc) to the relevant grid cell (grid amount / grid tot)
                            # Row - Y, Col - X
                            if X >= cols or Y >= rows:
                                logging.info('     Ignoring out of bounds value at row '
                                             + str(loop_count) + ' (' + str(X) + ',' + str(Y) + ')')
                            else:
                                grid[Y, X] += pop * bgval / bg_tot

                    else:
                        # there are no non-zero bg cells in this WAD, but we still have pop to disperse
                        dest_pop_leftover += dest_pop * pc / 100

            if dest_pop_leftover > 0:
                # we've been through all of our WADs and there is still undistributed population (hopefully unlikely)
                logging.info('     Destination {} had {} unallocated population',format(row, dest_pop_leftover))
                lost_pop += dest_pop_leftover

        logging.info('\n     created - Loop count: {0} in {1} seconds'.format(loop_count,
                                                                              round(time.time() - initialTime, 1)))
        logging.info('     Values total: {:.3f}  Dispersed Grid total: {:.3f}  Unallocated: {:.3f}'.format(sum(destination_data['inTravel']),
                                                                                      sum(sum(grid)), lost_pop))

        return grid

    def createGrid_LD(self, sb, locations, vals_array, locationIndex):

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
            possible_neighbours = locationIndex.possible_locations(loc_E, loc_N, LD)
            neighbours = []
            for poss in possible_neighbours:
                if poss != loc:  # not itself!
                    poss_E = locations['eastings'][poss]
                    poss_N = locations['northings'][poss]
                    dist = math.sqrt((loc_E - poss_E) ** 2 + (loc_N - poss_N) ** 2)
                    if dist <= LD:
                        neighbours.append((poss_E, poss_N, dist))  # add a tuple to the lsit

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
            X1 = int(((loc_E - avi + halfcell) - sb.projParams.background_bl_east) / csize)
            Y1 = int(((loc_N - avi + halfcell) - sb.projParams.background_bl_north) / csize)
            if X1 < 0:
                X1 = 0
            if Y1 < 0:
                Y1 = 0
            X2 = int(((loc_E + avi - halfcell) - sb.projParams.background_bl_east) / csize)
            Y2 = int(((loc_N + avi - halfcell) - sb.projParams.background_bl_north) / csize)
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
                        total_weight += weight
                        # a tuple of the cell location and weight
                        cell_list.append((x, y, weight))

            #  4. Populate the results grid
            #    each affected cell from 3 receives a portion of the location's population
            #      X,Y += pop * weight(x,y)
            if total_weight > 0:
                for (x, y, weight) in cell_list:
                    grid[y, x] += popval * weight / total_weight
            else:
                # no centroids are within the AVI, so put all of the pop in its own cell
                x = int((loc_E - sb.projParams.background_bl_east) / csize)
                y = int((loc_N - sb.projParams.background_bl_north) / csize)
                grid[y, x] += popval

        # all locations have been processed, grid is fully populated

        logging.info('     Values total: {:.3f}  Dispersed Grid total: {:.3f}'.format(sum(vals_array), sum(sum(grid))))

        return grid