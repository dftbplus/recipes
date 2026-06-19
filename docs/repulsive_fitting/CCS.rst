*********************************************
Fitting 2-body repulsive potentials using CCS
*********************************************

When using the Curvature Constrained Splines method (`CCS <https://pubs.acs.org/doi/10.1021/acs.jctc.0c01156>`_) 
to fit a repulsive potential we fit a 2-body potential
to, as closely as possible, reproducy the energy difference between a set of DFT reference
energies and the corresponding *electronic energies* from DFTB.

[Input: `recipes/docs/_archives/recipes/repulsives/ccs`]

This  tutorial revolves around a couple of scripts to handle the various tasks
involved in the fitting. 
Here we are assuming that the time consuming part of
selecting and computing the DFT reference data has allready been carried out. 
Furtermore, we also assume that the electronic parametes used to generate the 
co-called Slater-Koster table have been optimized (or been obtained by some other means).

Two data-sets are provided in the archive. One for Si and one ZnO. In both cases,
the fitting invloves the following steps:


.. tip:: 
  :class: info
   
  The folders ``Si`` and ``ZnO`` contains scripts to automatically execute the various steps and can be used to skip over certain part. Thses scripts are also useful as templates when using CCS to fit repulsive potentials for a new set of systems  


1. Calculte the electronic energies for each structure in the DFT training-set.
===============================================================================
We first need to perform single point DFTB+ calculations using the pre-calculated Slater-Koster table.
A convinient way to generate the geometry 
input file is to use the ``master_converter.py`` script provided with CCS 
which takes two arguements, the filename of the file to read and the filname 
of the file to write, e.g:

.. code-block::

  master_converter.py OUTCAR in.gen


.. caution::

  For the procedure used here, you must make sure to use a so-called 
  dummy spline in the Slater-Koster file (it is this block of the file that 
  we will obtain in steps below). 

  Example of a dummy spline:

  .. code-block:: 

    Spline
    1 0.2
    0. 0.1 0. 0. 0. 0. 0. 0.
    0.1 0.2 0. 0. 0. 0. 0. 0.


.. admonition:: Task
  :class: info

  Go to each sub-folder were there is DFT data (the DFT-data can be found in a folder ``RAW_DATA`` located in bot the ``Si`` and ``ZnO`` folders), convert the geometry to a
  DFTB+ readable format and perform a single-point DFTB+ calculation.



2. Collect the DFT en DFTB electroinc energies.
===============================================
Using the script ``ccs_build_db`` we can build two ase data-bases  
containing the DFT and DFTB data. These data-bases are used in the next step. 
We need a list pointing to the DFT data and the corresponding DFTB data. 
The code assume that the we provide locations of the DFT data (e.g. ``OUTCAR.gz`` or any other ASE readable file containing a geometry and corresponding energy) in the 
first column and the DFTB data in the second column (should be a ``results.tag`` file). If the data are stored  
in the same folders the first and second columns are the same.


.. admonition:: Task
  :class: info

  Create a folder with a file pointing to the training-data (``list``) and the
  execute the command:

  .. code-block::

    ccs_build_db DFTB list DFT.db DFTB.db 

  Here it is instructive to make several fits (in separate folders) using 
  diffrent sub-sets of the data. e.g seperate fits for each polymorphs or
  pairs of polymorhs.

  Example of a ``list`` file::

    ../../RAW_DATA/DIMER/3.20/VASP/out.xyz    ../../RAW_DATA/DIMER/3.20/DFTB/results.tag
    ../../RAW_DATA/DIMER/3.21/VASP/out.xyz    ../../RAW_DATA/DIMER/3.21/DFTB/results.tag
    ../../RAW_DATA/DIMER/3.22/VASP/out.xyz    ../../RAW_DATA/DIMER/3.22/DFTB/results.tag
    ../../RAW_DATA/DIMER/3.23/VASP/out.xyz    ../../RAW_DATA/DIMER/3.23/DFTB/results.tag
    ../../RAW_DATA/DIMER/3.24/VASP/out.xyz    ../../RAW_DATA/DIMER/3.24/DFTB/results.tag
    ../../RAW_DATA/DIMER/3.25/VASP/out.xyz    ../../RAW_DATA/DIMER/3.25/DFTB/results.tag
    ../../RAW_DATA/DIMER/3.26/VASP/out.xyz    ../../RAW_DATA/DIMER/3.26/DFTB/results.tag
    ../../RAW_DATA/DIMER/3.27/VASP/out.xyz    ../../RAW_DATA/DIMER/3.27/DFTB/results.tag
    ../../RAW_DATA/DIMER/3.28/VASP/out.xyz    ../../RAW_DATA/DIMER/3.28/DFTB/results.tag
    ../../RAW_DATA/DIMER/3.29/VASP/out.xyz    ../../RAW_DATA/DIMER/3.29/DFTB/results.tag
    ../../RAW_DATA/DIMER/3.30/VASP/out.xyz    ../../RAW_DATA/DIMER/3.30/DFTB/results.tag

.. tip::
  You can get a quick overview of the training-data using the Atomic Simulation 
  Environment (`ASE <https://wiki.fysik.dtu.dk/ase/>`_) by issuing the command: 

  .. code-block::
 
    ase-gui DFT.db 
   
  (see figure below)

  .. figure:: ase.png
      :alt: ase
      :width: 600
      :align: center

      With the Atomic Simulation Environment we can get an quick overview 
      of the training-data. We can, amongst much other, browse the structures 
      and plot the energies.    

.. tip:: 
  :class: info
   
  Examples files for the fitting task can be found in the folders ``EXAMPLE_FIT`` located in the ``Si`` and ``ZnO`` folders. 


3. Produce the specific training-file for CCS.
==============================================
We collect pair-wise distances from the structures stored in the two 
data-bases and create a file called ``structures.json`` that CCS 
use for the fitting.

.. admonition:: Task
  :class: info

  Stay in the same folder and execute:

  .. code-block::

    ccs_fetch DFTB 6.0 all DFT.db DFTB.db

  The arguments corresponds to, in order: 
  
  ``MODE cutoff_radius(Å) No_of_structures DFT_DATABASE DFTB_DATABASE``
  
  For repulsive potential fitting set ``MODE=DFTB``.

.. caution::

  Never use a cut-off radius that is smaller than used in the fitting (see next step).

4. Now we can do fitting! 
=========================
We provide the setting in a file ``input.json`` where we speicify the cut-off radius
the resolution of the spline and the type of constraints (rep = stricktly repulsive, 
sw=attractive at long distance and repulsive at short distance).

example for Si:

input.json::

  {
        "General": {
                "interface": "DFTB"
        },
        "Twobody": {
                "Si-Si": {
                        "Rcut": 3.5,
                        "Resolution": 0.05,
                        "Swtype": "sw"
                }
        }
  }

example for ZnO:

input.json::

    {
        "General": {
                "interface": "DFTB"
        },
        "Twobody": {
                "O-Zn": {
                        "Rcut": 6.0,
                        "Resolution": 0.25,
                        "Swtype": "rep"
                },
                "Zn-Zn": {
                        "Rcut": 6.0,
                        "Resolution": 0.25,
                        "Swtype": "rep"
                },
                "O-O": {
                        "Rcut": 6.0,
                        "Resolution": 0.25,
                        "Swtype": "rep"
                }
        }
    }




.. admonition:: Task
  :class: info

  Write the file ``input.json`` and execute:

  .. code-block::

   ccs_fit 

.. caution::

  Rcut must be smaller than the cut-off radius in the previus step!   

5. Enjoy succes!(?)  
===================
The quallity of the fit is provided in ``error.out`` and the resulting
parameters in ``CCS_params.json``.


6. Validation  
===================
The quallity of the generated CCS potential can be validated to unseen data using 
the ``ccs_validate`` command. 
The validation requires us to point out two ASE data-bases, one with the DFT and 
one with the corresponing DFTB data. 
You can use the same procedure as in Step 2 to create these data-bases.
 
.. admonition:: Task
  :class: info

  After collecting the validation data into two data-bases you can issue the command:

  .. code-block::

   ccs_validate DFTB 

A summary from the validation is stored in a file calles ``CCS_validation.dat``. 
A ASE data-base called ``CCS_validation.db`` is also generated and allows for a quick overview.

7. Convert to DFTB+ Slater-Koster format.
=========================================
DFTB+ have a specific format for the 2-body potential, a cubic 
spline-table appended at the end of the Slater-Koster file. We need
to convert the ``CCS_params.json`` file to this format.

.. admonition:: Task
  :class: info

  Execute: 

  .. code-block::

     ccs_export_sktable CCS_params.json

  The result are printed to files ``X-Y.spl`` where ``X`` and ``Y`` are
  the corresponding elements in the 2-body potential, e.g  
  ``X=Zn, Y=O``.

.. tip::

  You can use the ``splvalue`` executable provided with the dftb+ code to 
  create a a date-file presenting the 2-body spline in a format feasiable 
  for plotting using the following command: 

  .. code-block::

     splvalue Si-Si.spl > Spline.dat
    
  .. figure:: Reps.png
      :alt: ase
      :width: 600
      :align: center

      Comparative plot showing 2-body spline repulsives individually fitted to
      Si polymorphs of varying coordinations (2C-8C).    

8. Use the new parameters.
==========================
Replace the dummy-spline in the ``.skf`` file with the data from the ``.spl`` file generated in step 6 and 
you are good to go.

.. tip::
  You can use the ``ccs_build_db`` command to build a data-base with the new DFTB+ data for a quick overview
  of the performance of the new parameters.

.. caution::

   You should always validate your parameters by recalculating the energies for the structures in your traning-set *and* 
   for strucutres that was not included in the traning-set. The transferabily will have some limits which is especially 
   apparent in the case of Si using the electronic parameters provided with this recipie.




