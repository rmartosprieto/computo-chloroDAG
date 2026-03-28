#!/bin/bash

# Update package list
sudo apt-get update

# Basic system libraries (safe baseline)
sudo apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev