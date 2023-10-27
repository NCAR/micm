.. _JIT:

A JIT-compiled solver
=====================

Solving a chemical mechanism is a computationally intensive task. There are many reasons for this such as
multiple iterations by a solver to achieve an integration, stiff systems requiring internal substepping to 
acheive a numerically stable solution, and cache misses, among others. 

This tutorial focuses on alleviating the cache misses. A popular method for handling cache misses is to pre-compute the indices. 
This method, which may be referred to as ahead-of-time (AOT) compilation, is used in applications such as KPP :cite:`Damian2002`.
Pre-computed methods require code preprocessors and preclude runtime configurable software, which is a goal of micm.

MICM uses just-in-time (JIT) compiled functions built with `LLVM JIT <https://llvm.org/docs/tutorial/BuildingAJIT1.html>`_ libraries
to supply runtime-configurable chemistry to avoid cache misses in important chemistry functions.

Up until now, a :cpp:class:`micm::RosenbrockSolver` has been used. This is a special class in micm which builds all
of the componenets needed to solve chemistry in memory, including the forcing function and the jacobian. MICM
also provides a :cpp:class:`micm::JitRosenbrockSolver` which builds and compiles several important chemistry functions at runtime.

What are we JIT-compilng?
-------------------------

A list of compiled functions is below. How they are compiled and the ellided operations are beyond the scope of this tutorial, but you are free to inspect the source code yourself.

- :cpp:func:`micm::RosenbrockSolver::AlphaMinusJacobian` is JIT compiled and called with :cpp:func:`micm::JitRosenbrockSolver::AlphaMinusJacobian`
- :cpp:func:`micm::LuDecomposition::Decompose` is JIT compiled and called with :cpp:func:`micm::JitLuDecomposition::Decompose`
- :cpp:func:`micm::LinearSolver::Factor` and :cpp:func:`micm::LinearSolver::Solve` are JIT compiled and called with :cpp:func:`micm::JitLinearSolver::Factor` and :cpp:func:`micm::JitLinearSolver::Solve`
- :cpp:func:`micm::ProcessSet::AddForcingTerms` and :cpp:func:`micm::ProcessSet::AddJacobianTerms` are JIT compiled and called with :cpp:func:`micm::JitProcessSet::AddForcingTerms` and :cpp:func:`micm::JitProcessSet::AddJacobianTerms`

So what does this gain me?
--------------------------

Runtime configuraiton of chemical mechanisms that are *fast*. Let's compare the speed of :cpp:class:`micm::RosenbrockSolver` 
and the :cpp:class:`micm::JitRosenbrockSolver` classes applied to the same problem. We'll use this system.

.. math::

  A &\longrightarrow B, &k_{1, \mathrm{user\ defined}} \\
  2B &\longrightarrow B + C, &k_{2, \mathrm{user\ defined}} \\
  B + C &\longrightarrow A + C, \qquad &k_{3, \mathrm{user\ defined}} \\
