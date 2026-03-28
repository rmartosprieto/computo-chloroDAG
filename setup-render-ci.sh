#!/bin/bash

# Update package list
sudo apt-get update

# Numerical libraries required by numpy/scipy and causal methods
sudo apt-get install -y \
    libblas-dev \
    liblapack-dev \
    gfortran