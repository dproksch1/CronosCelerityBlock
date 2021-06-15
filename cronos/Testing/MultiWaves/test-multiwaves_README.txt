==================================
Manual for MultiWaves test problem
==================================
Created:      25sep2017 jk
Last changed: 12feb2018 jk
Reference: Section 8.5 of the Cronos reference paper (Kissmann+2018)

Procedure to run test case:
 1. create or link to data file
 2. set $poub via 'export poub=data' (or equivalent)
 3. default k range is 0.5 to 10 in steps of 0.50.
    To change, edit line #2 of file 'bat_spawn-catfiles.sh'
    (see also comments therein)
 4. execute file 'bat_spawn-catfiles.sh' to generate (in the default
    case) 20 separate *.cat files within 'data'
 5. execute file 'bat_run-all-sims.sh <np>' with <np> equal to the number of
    cores that you wish to use simulateneously, and confirm safety question
    to do just that. On my machine, this may take about 50 sec per simulation.
 6. execute file 'bat_condense-data.sh'. Data tripels (k, v_phase, damping)
    are now listed in file "multiwaves-results.dat" for further processing.
    Note: This step uses the Python program "analysis/wavefit.py", which may
    also be employed stand-alone to extract and/or plot (using the -showplot
    option) the By_max(t) profile from a single simulation.

Caution: Temporary files of name pattern _*.{cat|dat|txt} will be created,
	 and then deleted at the end, possibly overwriting pre-existing
	 files of the same name!

For debugging, use 'python analysis/wavefit.py <pname> [-showplot]'
