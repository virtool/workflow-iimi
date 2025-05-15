import csv
import gzip
import json
import tarfile
from collections import defaultdict
from pathlib import Path

from pydantic import BaseModel
from virtool_workflow.analysis.fastqc import NucleotidePoint





class PredictionCoverage(BaseModel):
    """The RLE coverage data for a sequence."""

    lengths: list[int]
    values: list[int]

class PredictionRaw(BaseModel):
    """A sequence-level prediction result.

    This is a simplified version of the ``Prediction`` class in the
    ``virtool_workflow.analysis.iimi`` module.
    """

    id: str
    coverage: PredictionCoverage
    result: bool
    probability: float
    untrustworthy_ranges: list[tuple[int, int]]

class PredictionSequence(PredictionRaw):
    """A sequence that was hit during map and either accepted or rejected as being present
    by IIMI.
    """
    length: int


class PredictionIsolate(BaseModel):
    """A virus isolate under which member sequences are organized. Isolates are nested
    below viruses.
    """

    id: str
    sequences: list[PredictionSequence]
    source_name: str
    source_type: str


class PredictionOTU(BaseModel):
    """A virus that was hit during mapping and either accepted or rejected by IIMI as being
    present.
    """

    id: str
    abbreviation: str
    isolates: list[PredictionIsolate]
    name: str
    result: bool


def calculate_nucleotide_composition(sequence: str) -> NucleotidePoint:
    """Calculate the nucleotide composition of the given sequence that contains
    the characters A, T, G, and C.

    Return the nucleotide composition as a ``NucleotidePoint`` object.

    :param sequence: the sequence to calculate the nucleotide composition of
    :return: the nucleotide composition
    """
    a = 0
    t = 0
    g = 0
    c = 0

    for char in sequence:
        match char:
            case "A" | "a":
                a += 1
            case "T" | "t":
                t += 1
            case "G" | "g":
                g += 1
            case "C" | "c":
                c += 1

    total = a + t + g + c

    a, t, g, c = map(lambda x: round(x / total, 6), (a, t, g, c))

    return NucleotidePoint(g=g, a=a, t=t, c=c)


def load_coverage(path: Path) -> dict[str, PredictionCoverage]:
    """Load IIMI coverage data from the given path.

    Coverage data is stored as a dictionary of sequence IDs to a
    ``PredictionHitSequenceCoverage``, which contains a list of RLE lengths and values.

    :param path: the path to the coverage CSV file
    :return: a dictionary of coverage data

    """
    coverage = {}

    with open(path) as f:
        reader = csv.reader(f, delimiter=",")
        next(reader)

        for sequence_id, lengths, values in reader:
            coverage[sequence_id] = PredictionCoverage(
                lengths=[int(l) for l in lengths.split(",")],
                values=[int(v) for v in values.split(",")],
            )

    return coverage


def load_untrustworthy_ranges(path: Path) -> dict[str, list[tuple[int, int]]]:
    """Load IIMI untrustworthy ranges from the given path.

    Untrustworthy ranges have been previously identified using a mappability profile
    generated as part of the model training preparation.

    Untrustworthy ranges are stored as a dictionary of sequence IDs to a list of
    tuples of start and end positions.

    :param path: the path to the untrustworthy ranges CSV file
    :return: a dictionary of untrustworthy ranges

    """
    untrustworthy_ranges = {}

    with open(path) as f:
        reader = csv.reader(f, delimiter=",")

        next(reader)

        for sequence_id, ranges in reader:
            sequence_id = str(sequence_id)

            if ranges:
                untrustworthy_ranges[sequence_id] = [
                    tuple(r.split("-")) for r in ranges.split(",")
                ]

    return untrustworthy_ranges


def load_and_format_prediction_results(
    reference_json_path: Path,
    reps_by_sequence_path: Path,
    output_path: Path,
) -> list[PredictionOTU]:
    """Load IIMI results into a list of ``PredictionHit`` objects.

    Results are gathered from the current working directory in the files:
    * ``coverage.csv``
    * ``prediction_sequence.csv``
    * ``prediction_virus.csv``
    * ``untrustworthy.csv``

    The results are annotated with virus annotation data collected from the provided
    ``reference.json.gz`` file.

    :param reference_json_path: the path to the reference.json.gz file
    :param reps_by_sequence_path: the path to the reps_by_sequence.csv file
    :param output_path: the path to the output directory
    :return: a list hits for each virus
    """
    coverage = load_coverage(output_path / "coverage.csv")
    untrustworthy_ranges = load_untrustworthy_ranges(output_path / "untrustworthy.csv")

    sequence_ids_by_reps = defaultdict(set)

    with open(reps_by_sequence_path) as f:
        reader = csv.reader(f, delimiter=",")
        next(reader)

        for row in reader:
            sequence_id, rep_id = row
            sequence_ids_by_reps[rep_id].add(sequence_id)

    predictions_by_seq_id: dict[str, PredictionRaw] = {}

    with open(output_path / "prediction_sequence.csv") as f:
        reader = csv.reader(f, delimiter=",")
        next(reader)

        for _, otu_id, isolate_id, rep_sequence_id, prediction, probability in reader:
            prediction_raw = PredictionRaw(
                id=rep_sequence_id,
                coverage=coverage[rep_sequence_id],
                probability=float(probability),
                result=prediction == "TRUE",
                untrustworthy_ranges=untrustworthy_ranges.pop(rep_sequence_id, []),
            )

            for sequence_id in sequence_ids_by_reps[rep_sequence_id]:
                predictions_by_seq_id[sequence_id] = prediction_raw

    result: list[PredictionOTU] = []

    with gzip.open(reference_json_path) as f:
        reference = json.load(f)

    for otu in reference["otus"]:
        otu_id = otu["_id"]

        prediction_otu = PredictionOTU(
            id=otu_id,
            abbreviation=otu["abbreviation"],
            isolates=[],
            name=otu["name"],
            result=False,
        )

        for isolate in otu["isolates"]:
            isolate_id = isolate["id"]

            prediction_isolate = PredictionIsolate(
                id=isolate_id,
                sequences=[                ],
                source_name=isolate["source_name"],
                source_type=isolate["source_type"],
            )

            had_prediction = False

            for sequence in isolate["sequences"]:
                prediction = predictions_by_seq_id.get(sequence["_id"])

                if prediction is None:
                    prediction = PredictionRaw(
                        id=sequence["_id"],
                        coverage=PredictionCoverage(
                            lengths=[len(sequence["sequence"])],
                            values=[0],
                        ),
                        probability=0.0,
                        result=False,
                        untrustworthy_ranges=[],
                    )
                else:
                    had_prediction = True

                prediction_isolate.sequences.append(
                    PredictionSequence(
                        id=prediction.id,
                        coverage=prediction.coverage,
                        length=len(sequence["sequence"]),
                        probability=prediction.probability,
                        result=prediction.result,
                        untrustworthy_ranges=prediction.untrustworthy_ranges
                    )
                )

            if had_prediction:
                prediction_otu.isolates.append(prediction_isolate)

        if not prediction_otu.isolates:
            continue

        prediction_otu.result = any(
            sequence.result
            for isolate in prediction_otu.isolates
            for sequence in isolate.sequences
        )

        result.append(prediction_otu)

    return result


def write_all_otu_fasta(reference_json_path: Path, output_path: Path):
    """Write a FASTA file containing all OTU sequences from the provided
    reference.json.gz file.

    FASTA records use sequence IDs as headers.

    :param reference_json_path: the path to the reference.json.gz file
    :param output_path: the path to write the FASTA file to

    """
    with gzip.open(reference_json_path, "rt") as f:
        data = json.load(f)

        with open(output_path, "w") as f:
            for otu in data["otus"]:
                for isolate in otu["isolates"]:
                    for sequence in isolate["sequences"]:
                        f.write(f">{sequence['_id']}\n{sequence['sequence']}\n")


def write_iimi_nucleotide_info(reference_json_path: Path, output_path: Path, logger):
    """Write the IIMI nucleotide information for each sequence to a CSV file.

    IIMI needs the following ordered columns to be present in the CSV file:
    * virus name
    * iso_id
    * seg_id
    * A_percent
    * T_percent
    * G_percent
    * GC_percent
    * seg_len

    :param reference_json_path:
    :param output_path:
    :return:
    """
    with gzip.open(reference_json_path, "r") as f_in, open(output_path, "w") as f_out:
        data = json.load(f_in)

        f_out.write(
            "\t".join(
                [
                    "virus name",
                    "iso_id",
                    "seg_id",
                    "A_percent",
                    "T_percent",
                    "G_percent",
                    "GC_percent",
                    "seg_len",
                ],
            ),
        )

        sequence_count = len(
            [s for v in data["otus"] for i in v["isolates"] for s in i["sequences"]],
        )

        logger.info(f"Writing nucleotide info for {sequence_count} sequences")

        for otu in data["otus"]:
            for isolate in otu["isolates"]:
                for sequence in isolate["sequences"]:
                    nucleotide_composition = calculate_nucleotide_composition(
                        sequence["sequence"],
                    )

                    f_out.write(
                        "\n"
                        + "\t".join(
                            str(x)
                            for x in [
                                otu["_id"],
                                isolate["id"],
                                sequence["_id"],
                                nucleotide_composition.a,
                                nucleotide_composition.t,
                                nucleotide_composition.g,
                                nucleotide_composition.c + nucleotide_composition.g,
                                len(sequence["sequence"]),
                            ]
                        ),
                    )


def untar(path: Path, target_path: Path):
    with tarfile.open(path, "r:gz") as tar:
        tar.extractall(target_path)
