.. _Contributing:

Contributing
============

For all proposed changes (bug fixes, new feature, documentation upates, etc.) please file
an `issue <https://github.com/NCAR/micm/issues/new/choose>`_ detailing what your ask is or what you intend to do.

The NCAR software developers will work with you on that issue to help recommend implementations, work through questions,
or provide other background knowledge.

Testing
-------

Any code that you contribute must include associated unit tests and an integration test of some kind. Our tests run across
various platforms and with different configurations automatically. By including tests that exercise your code, you help to
minimize the amount of time debugging platform portability issues. Further, new features must include at least an example
in our documentation.

MICM collects code coverage statistics which is displayed on our homepage in a little badge

.. image:: https://codecov.io/gh/NCAR/micm/branch/main/graph/badge.svg?token=ATGO4DKTMY
    :target: https://codecov.io/gh/NCAR/micm
    :alt: codecov badge

This coverage is also reported from the code cov bot in most PRs, like `this one <https://github.com/NCAR/micm/pull/105>`_.
The coverage amount from any PR must not drop because of code that you add. This means that added code is required to be tested.
We don't require your code to have 100% test coverage since some edge cases (like poor configurations) or platforms (like GPUs)
may not be testable in the github runners, but at minimum you must not remove test coverage where it exists.

Style guide
-----------
We (mostly) follow the `Google C++ style guide <https://google.github.io/styleguide/cppguide.html>`_. Please attempt to do 
the same to minimize comments on PRs. However, we are not dogmatic and reasonable exceptions are allowed, especially where they
help to simplify code or API intefaces.

After each PR is pulled, we have an automated github action that runs ``clang-format`` over our codebase, so don't worry too much
about formatting your code. Some of that format may be removed anyway. No format configuration is liked by everyone, but the use
of ``clang-format`` enforces some consistency which means moving from file to file shouldn't have a large congitive load on the developer.

Setting up developer environments
---------------------------------

Most developers work on micm in VS Code. For instructions on how to setup micm with Xcode or Visual Studio,
please see the :ref:`Debugging` section of the :ref:`Installation and usage` guide.

Building the documentation
--------------------------

All of our docs are stored in the ``docs`` directory. There are several python dependencies that need to install. To make 
this easy, we provide a ``requirements.txt`` file that you can use to install. 

Virtualenv
^^^^^^^^^^

From the root directory of micm:

.. code-block:: bash

  pip install virtualenv 
  virtualenv micm_env
  source micm_env/bin/activate
  cd docs 
  pip install -r requirements.txt
  cd .. && mkdir build && cd build
  cmake -DMICM_BUILD_DOCS=ON ..
  make docs
  open docs/sphinx/index.html 

Venv
^^^^

.. code-block:: bash

  python -m venv micm_env
  source micm_env/bin/activate
  cd docs 
  pip install -r requirements.txt
  cd .. && mkdir build && cd build
  cmake -DMICM_BUILD_DOCS=ON ..
  make docs
  open docs/sphinx/index.html 

Conda
^^^^^

.. code-block:: bash

  conda create --name micm_env python=3.8 -y
  conda activate micm_env
  cd docs 
  pip install -r requirements.txt
  cd .. && mkdir build && cd build
  cmake -DMICM_BUILD_DOCS=ON ..
  make docs
  open docs/sphinx/index.html 