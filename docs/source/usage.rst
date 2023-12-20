Usage
=====

.. _installation:

Installation
------------
neutralNEMO is available to install through conda

.. code-block:: console

   (.venv) $ conda install -c conda-forge neutralNEMO

or can be installed through pip

.. code-block:: console

   (.venv) $ pip install neutralNEMO


Example of use 
--------------

.. code-block:: Python

   from neutralNEMO.grid import build_nemo_hgrid, load_hgriddata, load_zgriddata






To retrieve a list of random ingredients,
you can use the ``lumache.get_random_ingredients()`` function:

.. autofunction:: lumache.get_random_ingredients

The ``kind`` parameter should be either ``"meat"``, ``"fish"``,
or ``"veggies"``. Otherwise, :py:func:`lumache.get_random_ingredients`
will raise an exception.

.. autoexception:: lumache.InvalidKindError

For example:

>>> import lumache
>>> lumache.get_random_ingredients()
['shells', 'gorgonzola', 'parsley']

