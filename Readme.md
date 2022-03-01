# Python version of Surface Builder 24/7

### Developed by University of Southampton GeoData Institute
### on behalf of ONS

v1.0 Jan-Mar 2022

A Python (3.x and above) redevelopment of the original VB code.

`src/` - sources files

> `main.py` - A sample main Python entry program
> 
> `sb247.py` - The main SB247 Class, the container for all functionality and data
> 
> `projectParams.py` - Stores all Project parameters
> 
> `modelRun.py` - Run the model including the main algorithms
> 
> `unit_tests.py`- Perform a series of validations (requires sample Data)

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
