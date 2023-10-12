.. _Installation and usage:

Installation and usage
======================

This tutorial is going to focus **only** and how to install micm and/or include it 
into your project in multiple different ways. Any API specific details will be elided and will be covered 
in another tutorial.


Github codespace
----------------

If you want to play around with micm without figuring out how to include it in your local setup, this is the **easiest** 
option. Github codespaces offers a cloud-hosted version of Visual Studio Code configured to work with specific projects.
We've set one up for micm. It'll allow you to instantly run the tests and make changes. Please note that there is a cap on
the number of hours your personal github account has each month. Follow these instructions to see your 
`github codespace usage <https://docs.github.com/en/billing/managing-billing-for-github-codespaces/viewing-your-github-codespaces-usage>`_.
At the time of this writing, there was a maximum of 120 core-houres allowed for github codespaces for free accounts.

To set this up, on the github page for `micm <https://github.com/NCAR/micm>`_, poke on the green code button and choose the
codespaces tab and select "create codespace on main". This will open up a new tab and start building a cloud environment 
running an instance of VSCode with the micm repository displayed and all of the tests prebuilt.

.. image:: images/codespaces.png

The first time that you open up a codespace, it will spend some time building the image and then compiling the test files.
Onces that's done, you can move into the build direcotry and run the tests.


.. code-block:: bash

  cd build
  make test

Installing
----------

From an archive
^^^^^^^^^^^^^^^

All versions of micm are associated with a github `release <https://github.com/NCAR/micm/releases>`_. 
Each release includes a tarall and zip that you can use to grab the code.

Find a release of micm that you want to build and download that archive. You can either do this with the browser by
poking on the desired file or with the commands below.

If you intend to use cmake to install micm, you can choose the install location when you configure
cmake: ``cmake -D CMAKE_INSTALL_PREFIX=/Users/me/Documents ..``.

Zip
~~~
.. code-block:: bash

  wget https://github.com/NCAR/micm/archive/refs/tags/v3.2.0.zip
  unzip v3.2.0.zip
  cd micm-3.2.0 
  mkdir build && cd build
  cmake ..
  make -j 8
  make test
  sudo make install

Tarball
~~~~~~~
.. code-block:: bash

  wget https://github.com/NCAR/micm/archive/refs/tags/v3.2.0.tar.gz
  tar -xf v3.2.0.tar.gz
  cd micm-3.2.0 
  mkdir build && cd build
  cmake ..
  make -j 8
  make test
  sudo make install

Cloning from github
^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

  git clone https://github.com/NCAR/micm.git
  cd micm
  mkdir build && cd build
  cmake ..
  make -j 8
  make test
  sudo make install

Usage after installation
^^^^^^^^^^^^^^^^^^^^^^^^

micm installs itself in a location typical on your system, like ``/usr/local``. It does so under it's own
directory with the version appended, e.g. ``/usr/local/micm-3.2.0``. It installs header files and files suitable
for use with cmake's `find_package <https://cmake.org/cmake/help/latest/command/find_package.html>`_.

::

  $ tree /usr/local/micm-3.2.0 -L 2
  /usr/local/micm-3.2.0
  ├── cmake
  │   ├── micmConfig.cmake
  │   ├── micmConfigVersion.cmake
  │   └── micm_Exports.cmake
  └── include
      └── micm

Specify include path
~~~~~~~~~~~~~~~~~~~~

To compile micm code, it's as simple as adding the include path to your compile command ``-I/usr/local/micm-3.2.0/include``
or ``export CPPFLAGS="-I/usr/local/micm-3.2.0/include"``. If you changed the install location when configuring cmake, you'll
need to set that path instead.

cmake with find_package
~~~~~~~~~~~~~~~~~~~~~~~



Cmake
-----

micm is developed with cmake support. This makes the inclusion of micm into projects that use cmake especially easy.

Fetch content
^^^^^^^^^^^^^

If you can use cmake 3.11+, the easiest way to include micm is with the 
`FetchContent  module <https://cmake.org/cmake/help/latest/module/FetchContent.html>`_.
If you must use a lower version, you'll either need to install the files on your system, or properly set
the include flags of your cmake project to point to the micm header files if you don't need GPU support.
You may also want to look into `ExternalProject <https://cmake.org/cmake/help/latest/module/ExternalProject.html>`_.

To use micm with fetch content, you'll need to include the module and point it to the micm repository and a commit
or tag that you want to use. Then you make the content available and link your cmake target to micm.

.. code-block:: cmake

  include(FetchContent)

  FetchContent_Declare(micm
    GIT_REPOSITORY https://github.com/NCAR/micm.git
    GIT_TAG        0996e5848b097e77ccbb2819f22c49844154f3e3
  )

  FetchContent_MakeAvailable(micm)

  add_executable(my_target my_target.cpp)

  target_link_libraries(my_target 
    PUBLIC 
      musica::micm
  )


Debugging
---------

VS Code
^^^^^^^

Xcode
^^^^^

Visual Studio
^^^^^^^^^^^^^