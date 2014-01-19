host-para-geo-sim
=================

This is a simple tool for generating host and parasite co-phylogenies
on a spatial grid.

More general documentation will be placed at http://phylo.bio.ku.edu/host-para-geo-sim

Usage
=====

    $ python sim-host-parasite.py

Will quit with an error and display the help message with more
instructions.

    $ python sim-host-parasite.py $RANDOM $RANDOM example.cfg

will launch a simulation with 2 random number seeds for the host and 
parasite pseudo random number generators respectively.

The run_sims.sh in example-scatterplot is a simple bash script
that generates 10 simulations from 4 different scenarios and
summarizes the simulation in a way that can be plotted using
the R script plot-grand-summary-4cond.R in that directory.
That invocation requires DendroPy as described on the main page mentioned above.


Authorship
==========

Software written by Mark T. Holder (2014) mtholder@gmail.com
as part of a collaboration with Dr. Janine Caira, Dr. Kirsten Jensen,
and Dr. Gavin J. P. Naylor.

Released under the GNU General Public License. V3. See LICENSE.TXT

