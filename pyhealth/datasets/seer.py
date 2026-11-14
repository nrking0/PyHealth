"""SEER cancer registry dataset.

The Surveillance, Epidemiology, and End Results (SEER) Program of the US
National Cancer Institute collects incidence and survival data for cancer
cases. Each record is one tumor (case); a patient with several primary
tumors has several records. This module turns a SEER*Stat case-listing
export, or SEER ASCII incidence files with their SAS layout, into one
canonical table, ``seer.csv``, which :class:`SEERDataset` loads.

The canonical columns follow the variables used by Hegselmann et al. (2018),
"Reproducible Survival Prediction with SEER Cancer Data", MLHC 2018.
"""

from __future__ import annotations

import logging
import re
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .base_dataset import BaseDataset

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class SEERField:
    """One canonical column of ``seer.csv``.

    Attributes:
        name: Canonical column name in ``seer.csv``.
        sas_name: Variable name in the SEER SAS layout (Nov 2016 research
            release), used to read ASCII incidence files.
        labels: SEER*Stat display names the column may have in a case-listing
            export. Matching ignores case and repeated whitespace.
        role: One of ``"identifier"``, ``"cohort"``, ``"geography"``,
            ``"diagnosis"`` or ``"outcome"``. Only ``"cohort"``,
            ``"geography"`` and ``"diagnosis"`` columns are candidate
            predictors; identifiers are metadata and outcomes are only used
            to build labels.
        required: Whether a source must provide the column. Optional columns
            that a source lacks (for example ``state_county``, which only
            SEER Research Plus includes) are written as empty.
    """

    name: str
    sas_name: str
    labels: tuple[str, ...]
    role: str
    required: bool = False


#: Canonical columns of ``seer.csv``, in file order (after ``case_id``).
#: SAS names follow the SEER 1973-2014 (November 2016) research data layout;
#: display names follow SEER*Stat. Both should be checked against the release
#: actually used, and differences documented.
SEER_FIELDS: tuple[SEERField, ...] = (
    # Identifiers (metadata, never predictors).
    SEERField("patient_id", "PUBCSNUM", ("Patient ID",), "identifier", True),
    SEERField("record_number", "REC_NO", ("Record number",), "identifier", True),
    SEERField("sequence_number", "SEQ_NUM", ("Sequence number",), "identifier"),
    # Cohort selection.
    SEERField("year_of_diagnosis", "YEAR_DX", ("Year of diagnosis",), "cohort", True),
    SEERField(
        "site_recode",
        "SITERWHO",
        ("Site recode ICD-O-3/WHO 2008", "Site recode ICD-O-3 2023 Revision"),
        "cohort",
    ),
    SEERField("primary_site", "PRIMSITE", ("Primary Site",), "cohort"),
    # Geography.
    SEERField("registry", "REG", ("SEER registry",), "geography"),
    SEERField("state_county", "ST_CNTY", ("State-county",), "geography"),
    # Diagnosis-time attributes.
    SEERField("sex", "SEX", ("Sex",), "diagnosis"),
    SEERField("age_at_diagnosis", "AGE_DX", ("Age at diagnosis",), "diagnosis"),
    SEERField("race", "RACE1V", ("Race recode (W, B, AI, API)",), "diagnosis"),
    SEERField(
        "marital_status", "MAR_STAT", ("Marital status at diagnosis",), "diagnosis"
    ),
    SEERField("histology", "HISTO3V", ("Histologic Type ICD-O-3",), "diagnosis"),
    SEERField("behavior", "BEHO3V", ("Behavior code ICD-O-3",), "diagnosis"),
    SEERField("grade", "GRADE", ("Grade",), "diagnosis"),
    SEERField("laterality", "LATERAL", ("Laterality",), "diagnosis"),
    SEERField("cs_tumor_size", "CSTUMSIZ", ("CS tumor size (2004+)",), "diagnosis"),
    SEERField("cs_extension", "CSEXTEN", ("CS extension (2004+)",), "diagnosis"),
    SEERField("cs_lymph_nodes", "CSLYMPHN", ("CS lymph nodes (2004+)",), "diagnosis"),
    SEERField("cs_mets_at_dx", "CSMETSDX", ("CS mets at dx (2004+)",), "diagnosis"),
    SEERField(
        "regional_nodes_positive",
        "EOD10_PN",
        ("Regional nodes positive (1988+)",),
        "diagnosis",
    ),
    SEERField(
        "regional_nodes_examined",
        "EOD10_NE",
        ("Regional nodes examined (1988+)",),
        "diagnosis",
    ),
    SEERField(
        "ajcc_stage_6th",
        "DAJCCSTG",
        ("Derived AJCC Stage Group, 6th ed (2004+)",),
        "diagnosis",
    ),
    SEERField(
        "first_malignant_primary",
        "FIRSTPRM",
        ("First malignant primary indicator",),
        "diagnosis",
    ),
    # Outcomes (only used to build labels).
    SEERField("survival_months", "SRV_TIME_MON", ("Survival months",), "outcome", True),
    SEERField(
        "survival_months_flag",
        "SRV_TIME_MON_FLAG",
        ("Survival months flag",),
        "outcome",
    ),
    SEERField(
        "cause_of_death_class",
        "VSRTSADX",
        (
            "SEER cause-specific death classification",
            "SEER cause of death classification",
        ),
        "outcome",
        True,
    ),
    SEERField(
        "other_cause_of_death_class",
        "ODTHCLASS",
        ("SEER other cause of death classification",),
        "outcome",
    ),
    SEERField(
        "vital_status",
        "STAT_REC",
        ("Vital status recode (study cutoff used)",),
        "outcome",
    ),
)

#: Every column of ``seer.csv``, in file order.
SEER_COLUMNS: tuple[str, ...] = ("patient_id", "case_id") + tuple(
    field.name for field in SEER_FIELDS if field.name != "patient_id"
)
#: Geographic columns, removed jointly in the geography ablation.
SEER_GEOGRAPHY_COLUMNS: tuple[str, ...] = tuple(
    field.name for field in SEER_FIELDS if field.role == "geography"
)
#: Outcome columns. They build labels and must never be predictors.
SEER_OUTCOME_COLUMNS: tuple[str, ...] = tuple(
    field.name for field in SEER_FIELDS if field.role == "outcome"
)
#: ``site_recode`` code for lung and bronchus (SEER site recode ICD-O-3/WHO 2008).
LUNG_SITE_RECODE = 22030
#: ``site_recode`` code for breast (SEER site recode ICD-O-3/WHO 2008).
BREAST_SITE_RECODE = 26000

_CSV_NAME = "seer.csv"
_SAS_LINE = re.compile(
    r"^\s*@\s*(\d+)\s+([A-Za-z_][A-Za-z0-9_]*)\s+\$?(?:char)?(\d+)\.\s*(.*)$",
    re.IGNORECASE,
)
_SAS_COMMENT = re.compile(r"/\*\s*(.*?)\s*\*/")


@dataclass(frozen=True)
class SASColumn:
    """One fixed-width column of a SEER ASCII file.

    Attributes:
        name: SAS variable name, e.g. ``"PUBCSNUM"``.
        start: 1-based start position.
        width: Width in characters.
        label: The comment after the variable, e.g. ``"Patient ID"``.
    """

    name: str
    start: int
    width: int
    label: str


def _normalize(name: str) -> str:
    return " ".join(str(name).split()).lower()


def _column_lookup() -> dict[str, str]:
    """Maps normalised canonical, SAS and display names to canonical names."""
    lookup: dict[str, str] = {}
    for field in SEER_FIELDS:
        for alias in (field.name, field.sas_name, *field.labels):
            lookup[_normalize(alias)] = field.name
    return lookup


def parse_sas_layout(sas_path: str | Path) -> list[SASColumn]:
    """Reads column positions from a SEER SAS input statement.

    Each SEER ASCII release ships a SAS program (for example
    ``read.seer.research.nov2016.sas``) whose ``input`` statement lists every
    fixed-width column, one per line::

        @ 1   PUBCSNUM   $char8.  /* Patient ID */

    Lines that do not have this shape are ignored.

    Args:
        sas_path: Path to the SAS program.

    Returns:
        The columns, in file order.

    Raises:
        FileNotFoundError: If ``sas_path`` does not exist.
        ValueError: If no column is found, or two columns share a name.

    Examples:
        >>> from pyhealth.datasets.seer import parse_sas_layout
        >>> parse_sas_layout("read.seer.research.nov2016.sas")[0]  # doctest: +SKIP
        SASColumn(name='PUBCSNUM', start=1, width=8, label='Patient ID')
    """
    path = Path(sas_path)
    if not path.exists():
        raise FileNotFoundError(f"SAS layout {path} not found")
    columns: list[SASColumn] = []
    for line in path.read_text(encoding="latin-1").splitlines():
        match = _SAS_LINE.match(line)
        if match is None:
            continue
        start, name, width, rest = match.groups()
        comment = _SAS_COMMENT.search(rest)
        columns.append(
            SASColumn(
                name=name.upper(),
                start=int(start),
                width=int(width),
                label=comment.group(1) if comment else "",
            )
        )
    if not columns:
        raise ValueError(f"{path}: no '@<start> <name> $char<width>.' lines found")
    names = [column.name for column in columns]
    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        raise ValueError(f"{path}: duplicate SAS variables {duplicates}")
    return columns


def _read_incidence(incidence_path: str | Path, sas_path: str | Path) -> pd.DataFrame:
    """Reads SEER ASCII files into a frame with canonical column names."""
    lookup = _column_lookup()
    wanted = []
    for column in parse_sas_layout(sas_path):
        name = lookup.get(_normalize(column.name)) or lookup.get(
            _normalize(column.label)
        )
        if name is not None:
            wanted.append((name, column))
    path = Path(incidence_path)
    if path.is_dir():
        files = sorted(
            p
            for p in path.iterdir()
            if p.is_file() and re.search(r"\.txt(\.gz)?$", p.name, re.IGNORECASE)
        )
        if not files:
            raise FileNotFoundError(f"No .TXT incidence files in {path}")
    elif path.exists():
        files = [path]
    else:
        raise FileNotFoundError(f"Incidence path {path} not found")
    frames = [
        pd.read_fwf(
            file,
            colspecs=[(c.start - 1, c.start - 1 + c.width) for _, c in wanted],
            names=[name for name, _ in wanted],
            dtype=str,
            header=None,
            compression="infer",
            encoding="latin-1",
        )
        for file in files
    ]
    return pd.concat(frames, ignore_index=True)


def _read_export(export_path: str | Path) -> pd.DataFrame:
    """Reads a SEER*Stat case listing (CSV or tab-separated) into a frame
    with canonical column names. Unknown columns are dropped."""
    path = Path(export_path)
    if not path.exists():
        raise FileNotFoundError(f"SEER*Stat export {path} not found")
    with open(path, encoding="latin-1") as handle:
        header = handle.readline()
    sep = "\t" if header.count("\t") > header.count(",") else ","
    table = pd.read_csv(path, sep=sep, dtype=str, encoding="latin-1")
    lookup = _column_lookup()
    renamed: dict[str, str] = {}
    for column in table.columns:
        name = lookup.get(_normalize(column))
        if name is None:
            continue
        if name in renamed.values():
            raise ValueError(f"{path}: two columns map to {name!r}")
        renamed[column] = name
    return table[list(renamed)].rename(columns=renamed)


def _read_cases(cases_path: str | Path) -> set[tuple[str, int]]:
    """Reads the (patient ID, record number) pairs of a SEER*Stat case list."""
    table = _read_export(cases_path)
    missing = {"patient_id", "record_number"} - set(table.columns)
    if missing:
        raise ValueError(
            f"{cases_path}: case list needs 'Patient ID' and 'Record number' "
            f"columns, missing {sorted(missing)}"
        )
    patients = table["patient_id"].str.strip()
    records = pd.to_numeric(table["record_number"], errors="raise").astype(int)
    return set(zip(patients, records))


def _to_codes(values: pd.Series, name: str, strict: bool) -> pd.Series:
    """Converts a column of SEER codes to numbers, keeping blanks as missing.

    Args:
        values: The column, as strings.
        name: Canonical column name, for error messages.
        strict: Raise if a non-blank value is not a number. Otherwise such
            values become missing.
    """
    text = values.astype("string").str.strip().replace("", pd.NA)
    numbers = pd.to_numeric(text, errors="coerce")
    bad = text.notna() & numbers.isna()
    if bad.any():
        examples = sorted(text[bad].unique())[:3]
        if strict:
            raise ValueError(
                f"{name} has non-numeric values {examples}. Export SEER*Stat "
                "case listings as codes, not labels."
            )
        logger.warning("%s: %d non-numeric values set to missing", name, bad.sum())
    if numbers.dropna().mod(1).eq(0).all():
        return numbers.astype("Int64")
    return numbers


def prepare_seer(
    root: str | Path,
    export_path: str | Path | None = None,
    incidence_path: str | Path | None = None,
    sas_path: str | Path | None = None,
    cases_path: str | Path | None = None,
    year_range: tuple[int, int] | None = None,
    site_recodes: Iterable[int] | None = None,
) -> Path:
    """Builds the canonical ``seer.csv`` from raw SEER files.

    Give exactly one source:

    - ``export_path``: a SEER*Stat case listing saved as CSV (or
      tab-separated text), exported with codes rather than labels. Columns
      are matched by SEER*Stat display name (e.g. ``"Survival months"``), SAS
      name or canonical name; others are dropped.
    - ``incidence_path`` and ``sas_path``: SEER ASCII incidence files (one
      ``.TXT``/``.TXT.gz`` file or a directory of them) and the SAS program
      that describes their layout, as in the November 2016 research release
      used by Hegselmann et al.

    Cases can then be restricted to the authors' cohort with ``cases_path``
    (a SEER*Stat case list with ``Patient ID`` and ``Record number``, as
    exported from their SEER*Stat sessions), and/or with ``year_range`` and
    ``site_recodes``.

    Values are stored as SEER codes. Blanks become empty cells; special
    codes such as ``999`` (unknown age) or ``9999`` (unknown survival
    months) are kept as they are, for the task to handle.

    Args:
        root: Directory to write ``seer.csv`` to.
        export_path: SEER*Stat case-listing export.
        incidence_path: SEER ASCII incidence file or directory.
        sas_path: SAS layout for ``incidence_path``.
        cases_path: Optional case list; only the listed cases are kept.
        year_range: Optional inclusive (first, last) year of diagnosis, e.g.
            ``(2004, 2009)``.
        site_recodes: Optional ``site_recode`` codes to keep, e.g.
            ``[LUNG_SITE_RECODE]``.

    Returns:
        Path to the written ``seer.csv``.

    Raises:
        FileNotFoundError: If an input file is missing.
        ValueError: If the sources are given wrongly, a required column is
            missing, an identifier is blank, a case appears twice, a
            required outcome is not numeric, or a survival time is negative.

    Examples:
        >>> from pyhealth.datasets.seer import LUNG_SITE_RECODE, prepare_seer
        >>> prepare_seer(  # doctest: +SKIP
        ...     "/data/seer",
        ...     export_path="/data/seer/lung_2004_2009.csv",
        ...     year_range=(2004, 2009),
        ...     site_recodes=[LUNG_SITE_RECODE],
        ... )
        PosixPath('/data/seer/seer.csv')
    """
    if (export_path is None) == (incidence_path is None):
        raise ValueError("Give exactly one of export_path or incidence_path")
    if incidence_path is not None:
        if sas_path is None:
            raise ValueError("incidence_path needs sas_path, the SAS layout")
        table = _read_incidence(incidence_path, sas_path)
    else:
        table = _read_export(export_path)

    missing = [f.name for f in SEER_FIELDS if f.required and f.name not in table]
    if missing:
        raise ValueError(f"Source is missing required columns {missing}")
    absent = [f.name for f in SEER_FIELDS if f.name not in table]
    if absent:
        logger.warning("Source has no %s; writing them as empty", absent)

    out = pd.DataFrame(index=table.index)
    out["patient_id"] = table["patient_id"].astype("string").str.strip()
    if out["patient_id"].isna().any() or out["patient_id"].eq("").any():
        raise ValueError("patient_id must not be blank")
    for field in SEER_FIELDS:
        if field.name == "patient_id":
            continue
        if field.name in table:
            out[field.name] = _to_codes(
                table[field.name], field.name, strict=field.required
            )
        else:
            out[field.name] = pd.array([pd.NA] * len(table), dtype="Int64")
    if out["record_number"].isna().any():
        raise ValueError("record_number must not be blank")
    if (out["survival_months"].dropna() < 0).any():
        raise ValueError("survival_months must not be negative")

    if cases_path is not None:
        cases = _read_cases(cases_path)
        keys = list(zip(out["patient_id"], out["record_number"].astype(int)))
        keep = pd.Series([key in cases for key in keys], index=out.index)
        found = {key for key, kept in zip(keys, keep) if kept}
        if len(found) < len(cases):
            logger.warning(
                "%d listed cases are not in the source", len(cases) - len(found)
            )
        out = out[keep]
    if year_range is not None:
        first, last = year_range
        out = out[out["year_of_diagnosis"].between(first, last).fillna(False)]
    if site_recodes is not None:
        codes = list(site_recodes)
        out = out[out["site_recode"].isin(codes).fillna(False)]
    if out.empty:
        logger.warning("No cases left after selection")

    out.insert(
        1,
        "case_id",
        out["patient_id"] + "_" + out["record_number"].astype("string"),
    )
    duplicated = out["case_id"].duplicated(keep=False)
    if duplicated.any():
        raise ValueError(
            "Each (patient_id, record_number) must appear once; duplicated: "
            f"{sorted(out.loc[duplicated, 'case_id'].unique())[:5]}"
        )
    out = out[list(SEER_COLUMNS)].reset_index(drop=True)

    root_path = Path(root)
    root_path.mkdir(parents=True, exist_ok=True)
    csv_path = root_path / _CSV_NAME
    out.to_csv(csv_path, index=False)
    logger.info(
        "Wrote %d SEER cases of %d patients to %s",
        len(out),
        out["patient_id"].nunique(),
        csv_path,
    )
    return csv_path


def _check_columns(csv_path: Path) -> None:
    """Raises ValueError if ``seer.csv`` lacks a canonical column."""
    header = pd.read_csv(csv_path, nrows=0).columns
    missing = [name for name in SEER_COLUMNS if name not in header]
    if missing:
        raise ValueError(
            f"{csv_path} is missing columns {missing}. Build it with "
            "pyhealth.datasets.seer.prepare_seer, which writes every column."
        )


class SEERDataset(BaseDataset):
    """SEER cancer registry dataset, one record per tumor.

    Loads the canonical ``seer.csv`` from ``root``. Each patient
    (``patient_id``, the SEER ``Patient ID``) has one ``seer`` event per
    tumor; ``case_id`` (patient ID and record number) identifies the tumor.
    Each event holds diagnosis-time attributes (registry, age, sex, site,
    histology, stage and extent of disease, ...) and outcomes (survival
    months, SEER cause-specific death classification). All values are SEER
    codes, kept as given, including special codes such as ``9999`` for
    unknown survival months. See :data:`SEER_FIELDS` for every column.

    If ``seer.csv`` is missing, it is built with :func:`prepare_seer` from
    ``export_path`` or from ``incidence_path`` and ``sas_path``.

    Paper:
        Hegselmann S, Greulich L, Varghese J, Dugas M. Reproducible Survival
        Prediction with SEER Cancer Data. MLHC 2018, PMLR 85:49-66.
        https://proceedings.mlr.press/v85/hegselmann18a.html

    Getting the data:
        - Source: SEER Research data, National Cancer Institute,
          https://seer.cancer.gov/data/access.html
        - Access: each user requests their own SEER*Stat account and agrees to
          the SEER Research Data Use Agreement. Patient-level data must not be
          shared or redistributed. County and state-county codes are only in
          SEER Research Plus.
        - Export: in SEER*Stat, make a case listing of the cohort (e.g. lung
          cancer diagnosed 2004-2009) with the variables in
          :data:`SEER_FIELDS`, and export it with codes rather than labels.
          The November 2016 release used by the paper is no longer offered
          in SEER*Stat; newer releases may rename variables.
        - ``root``: the folder that holds (or will hold) ``seer.csv``.

    Args:
        root: Directory holding ``seer.csv``.
        tables: Tables to load. Defaults to ``["seer"]``, the only one.
        dataset_name: Name of the dataset. Defaults to ``"seer"``.
        config_path: Path to the YAML config. Defaults to the bundled
            ``configs/seer.yaml``.
        export_path: SEER*Stat case listing used to build a missing
            ``seer.csv``.
        incidence_path: SEER ASCII incidence file(s) used to build a missing
            ``seer.csv``.
        sas_path: SAS layout for ``incidence_path``.
        cases_path: Optional case list restricting the built ``seer.csv``.
        year_range: Optional inclusive range of diagnosis years for the built
            ``seer.csv``.
        site_recodes: Optional ``site_recode`` codes for the built
            ``seer.csv``.
        **kwargs: Passed to :class:`~pyhealth.datasets.BaseDataset`.

    Raises:
        FileNotFoundError: If ``seer.csv`` is missing and no source is given.
        ValueError: If ``seer.csv`` lacks a canonical column.

    Examples:
        >>> from pyhealth.datasets import SEERDataset
        >>> dataset = SEERDataset(root="test-resources/seer")  # doctest: +SKIP
        >>> dataset.stats()  # doctest: +SKIP
    """

    def __init__(
        self,
        root: str,
        tables: list[str] | None = None,
        dataset_name: str | None = None,
        config_path: str | None = None,
        export_path: str | None = None,
        incidence_path: str | None = None,
        sas_path: str | None = None,
        cases_path: str | None = None,
        year_range: tuple[int, int] | None = None,
        site_recodes: Sequence[int] | None = None,
        **kwargs,
    ) -> None:
        csv_path = Path(root) / _CSV_NAME
        if not csv_path.exists():
            if export_path is None and incidence_path is None:
                raise FileNotFoundError(
                    f"{csv_path} not found. SEER data needs a SEER*Stat account "
                    "(https://seer.cancer.gov/data/access.html). Export a case "
                    "listing and pass export_path=..., or pass incidence_path= "
                    "and sas_path= for ASCII files."
                )
            prepare_seer(
                root,
                export_path=export_path,
                incidence_path=incidence_path,
                sas_path=sas_path,
                cases_path=cases_path,
                year_range=year_range,
                site_recodes=site_recodes,
            )
        _check_columns(csv_path)
        if config_path is None:
            config_path = str(Path(__file__).parent / "configs" / "seer.yaml")
        super().__init__(
            root=root,
            tables=tables or ["seer"],
            dataset_name=dataset_name or "seer",
            config_path=config_path,
            **kwargs,
        )
