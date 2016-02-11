#!/usr/bin/env bash

PIP_VERSION='8.0.2'
VENV_NAME='venv'

sudo apt-get install -y python3.4-venv
pyvenv-3.4 --without-pip $VENV_NAME
. ${VENV_NAME}/bin/activate
mkdir .setup_tmp
cd .setup_tmp
curl https://bootstrap.pypa.io/ez_setup.py -o - | python
wget "https://pypi.python.org/packages/source/p/pip/pip-${PIP_VERSION}.tar.gz"
tar xzvf pip-${PIP_VERSION}.tar.gz
cd pip-${PIP_VERSION}
python setup.py install
cd ..
cd ..
deactivate
