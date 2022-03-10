#
# Python version of Surface Builder 24/7
#
# Jan 2022
# GeoData Institute
# University of Southampton
# on behalf of ONS

# Import core modules

import logging

# A class for indexing and caching candidate locations within an area
#   These will later be checked using pythagoras, we are just cutting down the number to examine

LOC_COLS=20
LOC_ROWS=20

class LocationIndex:

    def __init__(self, projParams):  # could be replaced with BL east, north and origin_data[] to be more generic
        # instantiate a 2D grid across the background to hold lists of origins within each cell
        # and an empty dictionary for caching previous grid lookups
        self.loc_cache = {}
        self.loc_index = [[[] for col in range(LOC_COLS)] for row in range(LOC_ROWS)]

        # bottom left will match study area/background grid
        # calculate the size of each cell, rounded up if not a whole number
        self.size_x = round((projParams.background_tr_east - projParams.background_bl_east) / LOC_COLS + 0.5)
        self.size_y = round((projParams.background_tr_north - projParams.background_bl_north) / LOC_ROWS + 0.5)
        self.bl_east = projParams.background_bl_east
        self.bl_north = projParams.background_bl_north

        # insert all locations into their correct grid cell locations
        for location in range(0, len(projParams.origin_data['eastings'])):
            E = projParams.origin_data['eastings'][location]
            N = projParams.origin_data['northings'][location]
            # cell index:
            X = int((E - self.bl_east) / self.size_x)
            Y = int((N - self.bl_north) / self.size_y)

            # add the location index to the containing cell
            self.loc_index[Y][X].append(location)



    def possible_locations(self, Easting, Northing, Distance):
        # find locations within the bounding box around the location +- distance

        # BL cell index:
        X1 = int(((Easting - Distance) - self.bl_east) / self.size_x)
        Y1 = int(((Northing - Distance) - self.bl_north) / self.size_y)
        if X1 < 0:
            X1 = 0
        if Y1 < 0:
            Y1 = 0

        # TR cell index:
        X2 = int(((Easting + Distance) - self.bl_east) / self.size_x)
        Y2 = int(((Northing + Distance) - self.bl_north) / self.size_y)
        if X2 >= LOC_COLS:
            X2 = LOC_COLS-1
        if Y2 >= LOC_ROWS:
            Y2 = LOC_ROWS-1

        # encode BB as unique string for caching in the dictionary, eg. 2,3,15,19
        bbstr = '{},{},{},{}'.format(X1,X2,Y1,Y2)
        if bbstr in self.loc_cache:
            return self.loc_cache[bbstr]

        # add all stored grid locations from within the bounding box
        loc_list = []
        for x, y in [(x, y) for x in range(X1,X2) for y in range(Y1,Y2)]:
            loc_list.extend(self.loc_index[y][x])

        # add to the cache and return it
        self.loc_cache[bbstr] = loc_list
        return loc_list