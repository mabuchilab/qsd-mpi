# qsd-mpi #

A C++ library using quantum trajectories to solve quantum master equations.

© 2015-2016 Michael Goerz

`qsd-mpi` continues the development of [`qsd 1.3.5`][QSD],
© 1996-2004 Todd Brun, Rüdiger Schack

QSD is described in [Comput.Phys.Commun. 102 (1997) 210-228][paper],
also contained in the `doc` directory.

The further development of `qsd-mpi` has two long-term goals:

*   Introduce MPI-based parallelization to the code base
*   Adapt QSD to serve as a numerical backend for the [QNET][] library


[QSD]: https://www.ma.rhul.ac.uk/quantum-trajectories
[paper]: http://arxiv.org/abs/quant-ph/9608004
[QNET]: https://github.com/mabuchilab/QNET#qnet

## Installation and Usage ##

Run `make help` inside the project folder to receive help on compilation. You
may adapt the `Makefile` to your specific platform.

The `qsd-mpi` contains a number of example programs, which can be compiled with

    make <progname>

The `testprog` program specifically is intended to test the core functionality
of `qsd-mpi`. See [`TEST`](TEST) for details.

For use in your own program, `qsd-mpi` should first be compiled into a static
library, against which your program can then be linked. For example:

    ~/qsd-mpi> make libqsd.a
    ~/qsd-mpi> make install PREFIX=~/local/
    ~/qsd-mpi> cd ~
    ~> cp ./qsd-mpi/onespin.cc myprog.cc
    ~> g++ -O2 -I$HOME/local/include/qsd -o myprog myprog.cc -L$HOME/local/lib -lqsd

## Development ##

The development of `qsd-lib` is organized at
<https://github.com/mabuchilab/qsd-mpi>. Please submit bug reports or pull
requests there.

The project uses [semantic versioning](http://semver.org) and the [git-flow][]
branching model.

[git-flow]: http://nvie.com/posts/a-successful-git-branching-model/

If you encounter problems, or if you have questions, comments or suggestions
about QSD, please contact:

*   Todd Brun (<tbrun@usc.edu>)
*   Rüdiger Schack (<r.schack@rhul.ac.uk>)

For the further development of `qsd-mpi`, contact:

*   Michael Goerz (<goerz@stanford.edu>)
