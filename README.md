libpointmatcher is a modular point cloud meshing on top of [libpointmatcher].

libpointmatcher depends on:

 * [Eigen], a modern C++ matrix and linear-algebra library,
 * [libnabo], a fast K Nearest Neighbour library for low-dimensional spaces.
 * [libpointmatcher], a modular ICP library,
 
libpointmatcher is being developed by Andreas Breitenmoser and Stéphane Magnenat as part of our work at [ASL-ETH](http://www.asl.ethz.ch).


Compilation
-----------

libpointmesher uses [CMake] as build system.
Just create a directory, go inside it and type:

	cmake SRC_DIR
    
where `SRC_DIR` is the top-level directory of libpointmesher's sources.
If the dependencies are not installed system wide, you might have to tell [CMake] where to find them.
In case of doubt, read the [CMake documentation].

You first need to fetch and compile [libnabo] and [libpointmatcher].
To do so, you need [cmake], [git] and [Eigen].
On [Ubuntu], you can install these with `apt-get`:

	sudo apt-get install git-core cmake cmake-qt-gui libeigen2-dev

Then, you need to clone and build [libnabo]:

	git clone http://github.com/ethz-asl/libnabo
	mkdir build
	cd build
	cmake ..

Then, you need to clone and build [libpointmatcher]:

	git clone http://github.com/ethz-asl/libpointmatcher
	mkdir build
	cd build
	cmake ..

Launch `cmake-gui .` and specify the location of [libnabo]'s headers and static library.

Then, in the directory in which you are building libpointmesher, launch `cmake-gui .` and specify the location of [libnabo]'s and [libpointmatcher]'s headers and static libraries.


Test
----


Use [Paraview] to view the results.
On [Ubuntu], you can install [Paraview] with `apt-get`:
	sudo-apt get install paraview


Bug reporting
-------------

Please use [github's issue tracker](http://github.com/ethz-asl/libpointmesher/issues) to report bugs.


License
-------

libpointmesher is released under a permissive BSD license.

[Ubuntu]: http://www.ubuntu.com
[CMake]: http://www.cmake.org
[CMake documentation]: http://www.cmake.org/cmake/help/cmake2.6docs.html
[git]: http://git-scm.com
[Eigen]: http://eigen.tuxfamily.org
[libnabo]: http://github.com/ethz-asl/libnabo
[libpointmatcher]: http://github.com/ethz-asl/libpointmatcher
[Paraview]: http://www.paraview.org/
