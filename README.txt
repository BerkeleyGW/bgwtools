
    BerkeleyGW Tools
    Copyright (C) 2011 Felipe Homrich da Jornada


This python package provides a set of abstraction classes and tools to 
 convert, analyze and hack BerkeleyGW-related data.


INSTALL
-------

Add the following line to your ~/.bashrc or ~/.zshrc:
export PYTHONPATH=${PYTHONPATH}:$HOME/bgw-tools

where $HOME/bgw-tools is the location where this directory.


TEST
----

To test the installation, follow these steps:

1) Open a new terminal (or type ". ~/.bashrc" or ". ~/.zshrc")
2) Go to a directory with a BGW-compatible WFN
3) Lunch python
4) Type
from bgwtools.IO.wfn import wfnIO
wfn = wfnIO('WFN')
print wfn

5) Have fun :-)

