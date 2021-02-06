# prometheus 

We should include some description here.

**NOTE**: This repository make part of the ringer analysis kit (rak).

## Organization:

We should include some description here.

## How to run on singularity:

Download the prometheus local image to your home:
```bash
singularity pull docker://jodafons/ringer:base
```

After donwload it, just execute the `run` command:
```bash
singularity run --nv ringer_base.sif
```
**NOTE**: The argument `--nv` is responsible to attach the gpu into the container if you have it.

Inside of the singularity enviroment, just setup all prometheus dependencies using:
```bash
source /setup_all_here.sh ringer-atlas
```

## Install by hand:

After that, you can download the prometheus repository into your home:

```bash
git clone https://github.com/jodafons/prometheus.git
```
and install it:
```bash
cd prometheus && mkdir build && cd build && cmake .. && make && cd ..
source setup.sh
```

**NOTE**: Need to have ROOT and Gaugi installed.

## Framework status:

|  Branch    | Build Status |
| ---------- | ------------ |
|   Master   | [![Build Status](https://travis-ci.com/jodafons/prometheus.svg?branch=master)](https://travisci.org/jodafons/lorenzetti) |


## Notes:


