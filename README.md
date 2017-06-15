BerkeleyGW Tools
Copyright (C) 2011 Felipe Homrich da Jornada


This python package provides a set of abstraction classes and tools to 
convert, analyze and hack BerkeleyGW-related data.


INSTALL
-------

1) Extract this repository to a directory DIR. If you wish to clone the
   repository (recommended), do something like:

   $ cd ~/
   $ git clone https://user@bitbucket.org/jornada/berkeleygw-tools.git bgw-tools

   where user is your bitbucket username. In this example, DIR=$HOME/bgw-tools

2) Add the following line to your ~/.bashrc or ~/.zshrc:
   export PYTHONPATH=${PYTHONPATH}:DIR

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
5) If you don't see some neat stuff, see step (2) of the INSTALL part.
6) Happy BerkeleyGW hacking :-)

