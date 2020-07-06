# Python Mean-Squared Displacement Code

Usage: Requires python 3.7.4 or greater version (for pickle protocol 4)

usage: msd.py [-h] [-fname FNAME] [-comp COMP] [-N N] [-corr\_len CORR\_LEN]
              [-or_split OR_SPLIT] [-timestep TIMESTEP] [-skip SKIP]
              [-nframes NFRAMES] [-Lfile LFILE] [-op OP] [-unwrap UNWRAP]
              [-start_conf START_CONF] [-write_pckl WRITE_PCKL]

optional arguments:
  -h, --help            show this help message and exit
  -fname FNAME          The file name for the trajectory
  -comp COMP            This is the name of each component (repeatable)
  -N N                  The number of each component (repeatable)
  -corr\_len CORR\_LEN    Correlation length
  -or\_split OR\_SPLIT    Separation of origins
  -timestep TIMESTEP    Timestep (in ps)
  -skip SKIP            Lines to skip at beginning of each frame
  -nframes NFRAMES      Number of total frames
  -Lfile LFILE          File with box length
  -op OP                Options
  -unwrap UNWRAP        [1] if need to unwrap, [0] if already unwrapped
  -start\_conf START\_CONF
                        Starting configuration
  -write\_pckl WRITE\_PCKL
                        Write a pckl file with configs, [1] yes, [0] no



