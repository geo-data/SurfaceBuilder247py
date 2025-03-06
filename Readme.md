# Python version of Surface Builder 24/7

### Developed by University of Southampton GeoData Institute
### on behalf of ONS

v1.0 Jan-Mar 2022 - Initial version  
v1.1 Jul 2024 - Usability and file outputs improvements
v1.2 Mar 2025 - Bugfixes and optimisations Jan-Feb 2025

A Python (3.x and above) redevelopment of the original VB code.

For information about the project visit http://pop247.geodata.soton.ac.uk/

A user guide to Surface Builder 24/7 py is available from http://pop247.geodata.soton.ac.uk/software/SurfaceBuilder247Py_User_Guide.pdf

`src/` - sources files

> `sb247_modelRun.txt` - Pseudocode descriptions of the algorithms implemented
>
> `sb247_python.md` - Documentation of the Python structure and code
>
> `main.py` - A sample main Python entry program
> 
> `sb247.py` - The main SB247 Class, the container for all functionality and data
> 
> `projectParams.py` - Stores all Project parameters
> 
> `modelRun.py` - Run the model including the main algorithms
> 
> `locationIndex.py` - Structures and methods to rapidly identify data within a bounding area 
> 
> `gridCreate.py` - Create gridded data from model run outputs
> 
> `unit_tests.py` - Perform a series of validations (requires sample Data)

e.g. run:

```python main.py```

Unit tests:

```python unit_tests.py```

### Using Docker

`docker/` - Dockerfile to build a suitable Python container

e.g. run:

```docker run -it --rm --name sb2472py_docker -v "$PWD":/opt -w /opt sb247-docker:latest ./main.py```

Unit tests (requires sample Data):

```docker run -it --rm --name sb2472py_docker -v "$PWD":/opt -w /opt sb247-docker:latest python ./unit_tests.py```

### Using published Python package

Package home at PyPI: https://pypi.org/project/SurfaceBuilder247/

```pip install SurfaceBuilder247```

In `Main.py`:

```from SurfaceBuilder247 import SB247```
