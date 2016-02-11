#!/usr/bin/env bash

. venv/bin/activate

# wx
sudo apt-get install \
  dpkg-dev build-essential \
  swig python3-dev \
  libwebkit-dev \
  libjpeg-dev \
  libtiff-dev \
  checkinstall \
  freeglut3 \
  freeglut3-dev \
  libgtk2.0-dev \
  libsdl1.2-dev \
  libgstreamer-plugins-base0.10-dev
pip install --upgrade --trusted-host wxpython.org --pre -f http://wxpython.org/Phoenix/snapshot-builds/ wxPython_Phoenix

# numpy, et al
pip install -r requirements.txt

deactivate
