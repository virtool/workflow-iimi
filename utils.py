import csv
import gzip
import json
import tarfile
from collections import defaultdict
from pathlib import Path

from pydantic import BaseModel
from virtool_workflow.analysis.fastqc import NucleotidePoint


class PredictionHitSequenceCoverage(BaseModel):
    """The RLE coverage data for a sequence."""

    lengths: list[int]
    values: list[int]


class PredictionHitSequence(BaseModel):
    """
    A sequence that was hit during map and either accepted or rejected as being present
    by IIMI.
    """

    id: str
    coverage: PredictionHitSequenceCoverage
    length: int
    result: bool
    untrustworthy_ranges: list[tuple[int, int]]


class PredictionHitIsolate(BaseModel):
    """
    A virus isolate under which member sequences are organized. Isolates are nested
    below viruses.
    """

    id: str
    sequences: list[PredictionHitSequence]
    source_name: str
    source_type: str


class PredictionHit(BaseModel):
    """
    A virus that was hit during mapping and either accepted or rejected by IIMI as being
    present.
    """

    id: str
    abbreviation: str
    isolates: list[PredictionHitIsolate]
    name: str
    result: bool


def calculate_nucleotide_composition(sequence: str) -> NucleotidePoint:
    """
    Calculate the nucleotide composition of the given sequence that contains
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


def load_coverage(path: Path) -> dict[str, PredictionHitSequenceCoverage]:
    """
    Load IIMI coverage data from the given path.

    Coverage data is stored as a dictionary of sequence IDs to a
    ``PredictionHitSequenceCoverage``, which contains a list of RLE lengths and values.

    :param path: the path to the coverage CSV file
    :return: a dictionary of coverage data

    """
    coverage = {}

    with open(path, "r") as f:
        reader = csv.reader(f, delimiter=",")
        next(reader)

        for sequence_id, lengths, values in reader:
            coverage[sequence_id] = PredictionHitSequenceCoverage(
                lengths=[int(l) for l in lengths.split(",")],
                values=[int(v) for v in values.split(",")],
            )

    return coverage


def load_untrustworthy_ranges(path: Path) -> dict[str, list[tuple[int, int]]]:
    """
    Load IIMI untrustworthy ranges from the given path.

    Untrustworthy ranges have been previously identified using a mappability profile
    generated as part of the model training preparation.

    Untrustworthy ranges are stored as a dictionary of sequence IDs to a list of
    tuples of start and end positions.

    :param path: the path to the untrustworthy ranges CSV file
    :return: a dictionary of untrustworthy ranges

    """
    untrustworthy_ranges = {}

    with open(path, "r") as f:
        reader = csv.reader(f, delimiter=",")

        next(reader)

        for sequence_id, ranges in reader:
            sequence_id = str(sequence_id)

            if ranges:
                untrustworthy_ranges[sequence_id] = [
                    tuple(r.split("-")) for r in ranges.split(",")
                ]

    return untrustworthy_ranges


def load_virus_annotations(
    path: Path,
) -> tuple[dict[str, dict], dict[str, dict], dict[str, dict]]:
    """
    Load the virus annotations from the given path.

    :param path: the path to the virus annotations
    :type path: str

    :return: a dictionary of virus annotations
    :rtype: dict

    """
    otu_annotations = {}
    isolate_annotations = {}
    sequence_annotations = {}

    with gzip.open(path, "rt") as f:
        data = json.load(f)

        for otu in data["otus"]:
            otu_annotations[otu["_id"]] = {
                "abbreviation": otu["abbreviation"],
                "name": otu["name"],
                "schema": otu["schema"],
            }

            for isolate in otu["isolates"]:
                isolate_annotations[isolate["id"]] = {
                    "source_name": isolate["source_name"],
                    "source_type": isolate["source_type"],
                }

                for sequence in isolate["sequences"]:
                    sequence_annotations[sequence["_id"]] = {
                        "accession": sequence["accession"],
                        "definition": sequence["definition"],
                        "length": len(sequence["sequence"]),
                        "segment": sequence.get("segment"),
                    }

        return otu_annotations, isolate_annotations, sequence_annotations


def load_and_format_prediction_results(
    reference_json_path: Path, output_path: Path
) -> list[PredictionHit]:
    """
    Load IIMI results into a list of ``PredictionHit`` objects.

    Results are gathered from the current working directory in the files:
    * ``coverage.csv``
    * ``prediction_sequence.csv``
    * ``prediction_virus.csv``
    * ``untrustworthy.csv``

    The results are annotated with virus annotation data collected from the provided
    ``reference.json.gz`` file.

    :param reference_json_path:
    :param output_path:
    :return: a list hits for each virus
    """
    coverage = load_coverage(output_path / "coverage.csv")
    untrustworthy_ranges = load_untrustworthy_ranges(output_path / "untrustworthy.csv")

    (
        otu_annotations,
        isolate_annotations,
        sequence_annotations,
    ) = load_virus_annotations(reference_json_path)

    sequence_prediction_results = defaultdict(lambda: defaultdict(list))

    with open(output_path / "prediction_sequence.csv", "r") as f:
        reader = csv.reader(f, delimiter=",")

        next(reader)

        for prediction, sequence_id, isolate_id, otu_id, _ in reader:
            annotation = sequence_annotations.pop(sequence_id)

            sequence_prediction_results[otu_id][isolate_id].append(
                PredictionHitSequence(
                    id=sequence_id,
                    accession=annotation["accession"],
                    coverage=coverage.pop(sequence_id),
                    definition=annotation["definition"],
                    length=annotation["length"],
                    result=prediction == "TRUE",
                    segment=annotation["segment"],
                    sequence_id=sequence_id,
                    untrustworthy_ranges=untrustworthy_ranges.pop(sequence_id, []),
                )
            )

    result = []

    with open(output_path / "prediction_virus.csv", "r") as f:
        reader = csv.reader(f, delimiter=",")
        next(reader)

        for _, otu_id, prediction in reader:
            prediction = PredictionHit(
                id=otu_id,
                abbreviation=otu_annotations[otu_id]["abbreviation"],
                isolates=[],
                name=otu_annotations[otu_id]["name"],
                result=prediction == "TRUE",
            )

            prediction.isolates = [
                PredictionHitIsolate(
                    id=isolate_id,
                    sequences=sequence_prediction_results[otu_id][isolate_id],
                    source_name=isolate_annotations[isolate_id]["source_name"],
                    source_type=isolate_annotations[isolate_id]["source_type"],
                )
                for isolate_id in sequence_prediction_results[otu_id]
            ]

            result.append(prediction)

    return result


def write_all_otu_fasta(reference_json_path: Path, output_path: Path):
    """
    Write a FASTA file containing all OTU sequences from the provided
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
    """
    Write the IIMI nucleotide information for each sequence to a CSV file.

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
                ]
            )
        )

        sequence_count = len(
            [s for v in data["otus"] for i in v["isolates"] for s in i["sequences"]]
        )

        logger.info(f"Writing nucleotide info for {sequence_count} sequences")

        for otu in data["otus"]:
            for isolate in otu["isolates"]:
                for sequence in isolate["sequences"]:
                    nucleotide_composition = calculate_nucleotide_composition(
                        sequence["sequence"]
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
                        )
                    )


def untar(path: Path, target_path: Path):
    with tarfile.open(path, "r:gz") as tar:
        tar.extractall(target_path)
