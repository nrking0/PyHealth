pyhealth.datasets.SEERDataset
=============================

Overview
--------

The `Surveillance, Epidemiology, and End Results (SEER) Program
<https://seer.cancer.gov/>`_ of the US National Cancer Institute records the
incidence and survival of cancer cases. ``SEERDataset`` loads SEER research
data with one record per tumor: a patient (``patient_id``, the SEER
``Patient ID``) has one ``seer`` event per primary tumor, identified by
``case_id`` (patient ID and record number). Each event holds diagnosis-time
attributes and the outcomes used to build survival labels.

The columns follow the variables used by Hegselmann et al., `Reproducible
Survival Prediction with SEER Cancer Data
<https://proceedings.mlr.press/v85/hegselmann18a.html>`_ (MLHC 2018), whose
code is at `stefanhgm/MLHC2018-reproducible-survival-seer
<https://github.com/stefanhgm/MLHC2018-reproducible-survival-seer>`_.

.. code-block:: python

    from pyhealth.datasets import SEERDataset
    from pyhealth.datasets.seer import LUNG_SITE_RECODE

    # From a SEER*Stat case-listing export (builds root/seer.csv once).
    dataset = SEERDataset(
        root="/data/seer",
        export_path="/data/seer/lung_2004_2009.csv",
        year_range=(2004, 2009),
        site_recodes=[LUNG_SITE_RECODE],
    )

    # The synthetic fixture, no SEER access needed.
    dataset = SEERDataset(root="test-resources/seer")
    dataset.stats()

Getting the data
----------------

SEER data cannot be downloaded or redistributed by PyHealth.

1. Each user requests their own SEER Research account at
   https://seer.cancer.gov/data/access.html and agrees to the SEER Research
   Data Use Agreement. Patient-level data must not be shared, linked to other
   data, or committed to a repository. County and state-county codes are only
   in SEER Research Plus, which needs an institutional account.
2. In SEER*Stat (a Windows application), make a case listing of the cohort,
   for example lung and bronchus cancer diagnosed 2004-2009, with the
   variables in the column table below. The authors' SEER*Stat sessions
   describe their cohort selection.
3. Export the listing as CSV or tab-separated text **with codes, not
   labels**, and keep it outside version control.
4. Build ``seer.csv`` with :func:`~pyhealth.datasets.seer.prepare_seer` or by
   passing ``export_path`` to ``SEERDataset``. A case list
   (``Patient ID``, ``Record number``) can be passed as ``cases_path`` to
   keep exactly the authors' cases.

SEER ASCII incidence files with a SAS layout program (the format of the
November 2016 release used by the paper) can be read with ``incidence_path``
and ``sas_path``. SEER no longer distributes that release or ASCII files, so
newer releases are read from SEER*Stat exports. Variable names, codes and
cohort sizes can differ between releases; record the release used and
compare cohort counts with the paper.

Columns of ``seer.csv``
-----------------------

All values are SEER codes stored as integers, except ``patient_id`` and
``case_id``. Blanks are empty. Special codes are kept as they are, for
example ``999`` for unknown age and ``9999`` for unknown survival months, so
that the task decides how to treat them. Columns that a source lacks (for
example ``state_county`` in SEER Research exports) are written empty, with a
warning.

.. list-table::
   :header-rows: 1
   :widths: 26 18 38 18

   * - Column
     - SAS name (Nov 2016)
     - SEER*Stat name
     - Role
   * - ``patient_id``
     - ``PUBCSNUM``
     - Patient ID
     - identifier
   * - ``case_id``
     - —
     - (patient ID and record number)
     - identifier
   * - ``record_number``
     - ``REC_NO``
     - Record number
     - identifier
   * - ``sequence_number``
     - ``SEQ_NUM``
     - Sequence number
     - identifier
   * - ``year_of_diagnosis``
     - ``YEAR_DX``
     - Year of diagnosis
     - cohort
   * - ``site_recode``
     - ``SITERWHO``
     - Site recode ICD-O-3/WHO 2008
     - cohort
   * - ``primary_site``
     - ``PRIMSITE``
     - Primary Site
     - cohort
   * - ``registry``
     - ``REG``
     - SEER registry
     - geography
   * - ``state_county``
     - ``ST_CNTY``
     - State-county
     - geography
   * - ``sex``
     - ``SEX``
     - Sex
     - diagnosis
   * - ``age_at_diagnosis``
     - ``AGE_DX``
     - Age at diagnosis
     - diagnosis
   * - ``race``
     - ``RACE1V``
     - Race recode (W, B, AI, API)
     - diagnosis
   * - ``marital_status``
     - ``MAR_STAT``
     - Marital status at diagnosis
     - diagnosis
   * - ``histology``
     - ``HISTO3V``
     - Histologic Type ICD-O-3
     - diagnosis
   * - ``behavior``
     - ``BEHO3V``
     - Behavior code ICD-O-3
     - diagnosis
   * - ``grade``
     - ``GRADE``
     - Grade
     - diagnosis
   * - ``laterality``
     - ``LATERAL``
     - Laterality
     - diagnosis
   * - ``cs_tumor_size``
     - ``CSTUMSIZ``
     - CS tumor size (2004+)
     - diagnosis
   * - ``cs_extension``
     - ``CSEXTEN``
     - CS extension (2004+)
     - diagnosis
   * - ``cs_lymph_nodes``
     - ``CSLYMPHN``
     - CS lymph nodes (2004+)
     - diagnosis
   * - ``cs_mets_at_dx``
     - ``CSMETSDX``
     - CS mets at dx (2004+)
     - diagnosis
   * - ``regional_nodes_positive``
     - ``EOD10_PN``
     - Regional nodes positive (1988+)
     - diagnosis
   * - ``regional_nodes_examined``
     - ``EOD10_NE``
     - Regional nodes examined (1988+)
     - diagnosis
   * - ``ajcc_stage_6th``
     - ``DAJCCSTG``
     - Derived AJCC Stage Group, 6th ed (2004+)
     - diagnosis
   * - ``first_malignant_primary``
     - ``FIRSTPRM``
     - First malignant primary indicator
     - diagnosis
   * - ``survival_months``
     - ``SRV_TIME_MON``
     - Survival months
     - outcome
   * - ``survival_months_flag``
     - ``SRV_TIME_MON_FLAG``
     - Survival months flag
     - outcome
   * - ``cause_of_death_class``
     - ``VSRTSADX``
     - SEER cause-specific death classification
     - outcome
   * - ``other_cause_of_death_class``
     - ``ODTHCLASS``
     - SEER other cause of death classification
     - outcome
   * - ``vital_status``
     - ``STAT_REC``
     - Vital status recode (study cutoff used)
     - outcome

Only ``cohort``, ``geography`` and ``diagnosis`` columns are candidate
predictors. Identifiers are metadata, and outcome columns build labels and
must never be used as predictors. Treatment variables (surgery, radiation)
and other post-diagnosis variables are not loaded. The names above should
be checked against the release actually used.

Synthetic fixture
-----------------

``test-resources/seer`` holds invented data only: a fixed-width
``incidence/SYNTH.TXT`` with its layout ``read.seer.synthetic.sas``, a case
list ``CASES.csv``, and the resulting ``seer.csv`` for the lung 2004-2009
cohort (five patients, six tumors). It covers a patient with two tumors, a
cancer death at 5 months, survival of exactly 12 months, a death from
another cause at 6 months, censoring at 3 months, unknown survival
(``9999``) and age (``999``), blank values, and, in the raw file only, a case
diagnosed in 2012 and a breast cancer case.

.. autoclass:: pyhealth.datasets.SEERDataset
    :members:
    :undoc-members:
    :show-inheritance:

.. autofunction:: pyhealth.datasets.seer.prepare_seer

.. autofunction:: pyhealth.datasets.seer.parse_sas_layout

.. autoclass:: pyhealth.datasets.seer.SEERField

.. autoclass:: pyhealth.datasets.seer.SASColumn
