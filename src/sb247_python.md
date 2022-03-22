# Python version of Surface Builder 24/7

### Developed by University of Southampton GeoData Institute
### on behalf of ONS

v1.0 Jan-Mar 2022

A Python (3.x and above) redevelopment of the original VB code.

For information about the project visit http://pop247.geodata.soton.ac.uk/

## Documentation of the Python structure and code


### `sb247_modelRun.txt` - Pseudocode descriptions

This file lays out detailed, but readable step-by-step explanations of
the algorithms implemented in the source code,


### `main.py` - A sample entry program

The sample main.py Python script contains all of the steps needed to
load data, run the model and save outputs.

This script should be examined to see how to use the tool and adapted
for the purposes of any specific requirement.


### `sb247.py` - The main SB247 Class

The SB247 Class provides the container for all functionality and data.

The main.py or other entry script instantiates an object of the SB247
Class, implemented within this script and used as the container for
all other components and Subclasses.


### `projectParams.py` - Stores all Project parameters

Objects instantiated from the SB247 Class have a `.projParams`
Subclass, instantiated for storing all project data loaded, and
available to run the model.

The methods to load Origin, Destination, Timeseries and Background
data and the data structures to store them are implemented within the
Subclass defined in this script.


### `modelRun.py` - Run the model algorithms

Once data has been loaded, the SB247 object is ready to run the
SurfaceBuilder model, with a `.modelRun` Subclass instantiated to hold
run parameters and model outputs.

The methods to run the model, and data structures to hold output data,
are implemented within this script.


### `gridCreate.py` - Create gridded data

Once output data has been created by executing a model run, production
of gridded outputs is implemented by the GridCreate Class, with
methods to generate 2D (Numpy) arrays of:

* raw outputs (remaining and immobile Origins, on-site Destinations)

* in-travel Destinations, which are combined with the weighted
background grid

* Locally Dispersed outputs

These are saved to disk as ESRI format ASCII grid (.asc) files and in
csv format.  With the latter, only populated cells are exported, each
row representing a single cell centre, with grid reference and column
values for each output data type.


### `locationIndex.py` - Rapidly identify data within a bounding area 

The grid creation process frequently needs to find all locations of
various types of data which are within a radius (e.g. largest WAD or
Local Dispersion radius).

The LocationIndex Class, whose structures and methods are implemented
in this script optimises this process.  An instance is instantiated
for a particular type of data, with a 2D grid created, containing the
indexes of locations within the area represented.

When the possible_locations() method is called, a bounding box is
calculated across the grid and all possible locations are identified
from the cells populated within.

A Python dictionary is then used to cache and then easily look up the
output from previous requests across the same cell references.


### `unit_tests.py` - Perform a series of validations

The unit tests script runs the entire process of loading data, running
the model for a small test area, and producing all outputs.  Each
stage and output is compared with expected values, ensuring the tool
is operating correctly in a particular environment and has with valid
input data available.