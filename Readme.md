# Python version of Surface Builder 24/7

### Developed by University of Southampton GeoData Institute
### on behalf of ONS

v1.0 Jan-Mar 2022

A Python (3.x and above) redevelopment of the original VB code.

`src/` - sources files

> `main.py` - sample main Python entry program
> 
> `sb247.py` - ...
> 
> `projectParams.py` - ...
> 
> `modelRun.py` - ...
> 
> `unit_tests.py`- ... (requires sample Data)

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

