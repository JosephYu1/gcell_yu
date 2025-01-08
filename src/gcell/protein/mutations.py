"""
Module for protein sequence mutations.

Classes
-------
ProteinSequenceManipulator: Class for manipulating protein sequences
SequenceMutation: Class to store information about a sequence mutation
"""

import random
from copy import deepcopy
from dataclasses import dataclass
from enum import Enum


class MutationType(Enum):
    DELETION = "deletion"
    SHUFFLE = "shuffle"
    ALANINE_SCAN = "alanine_scan"
    RANDOM_AA = "random_aa"


@dataclass
class SequenceMutation:
    """Class to store information about a sequence mutation"""

    original_sequence: str
    mutated_sequence: str
    mutation_type: MutationType
    start_idx: int
    end_idx: int
    mutated_region: str
    original_region: str


class ProteinSequenceManipulator:
    """Class for manipulating protein sequences

    Methods
    -------
    validate_sequence: Validate that sequence only contains valid amino acids
    delete_range: Delete a range of amino acids from the sequence
    shuffle_range: Shuffle amino acids within a range
    random_aa_range: Replace a range with random amino acids
    alanine_scan_range: Replace a range with alanines (A)
    """

    # Standard amino acids
    AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")

    def __init__(self, sequence: str):
        """
        Initialize with a protein sequence

        Parameters
        ----------
        sequence (str): Original protein sequence
        """
        self.original_sequence = sequence
        self.validate_sequence(sequence)

    @staticmethod
    def validate_sequence(sequence: str) -> None:
        """
        Validate that sequence only contains valid amino acids

        Parameters
        ----------
        sequence (str): Sequence to validate

        Raises
        ------
            ValueError: If invalid amino acids are found
        """
        invalid_chars = set(sequence.upper()) - set(
            ProteinSequenceManipulator.AMINO_ACIDS
        )
        if invalid_chars:
            raise ValueError(f"Invalid amino acids found in sequence: {invalid_chars}")

    def delete_range(self, start_idx: int, end_idx: int) -> SequenceMutation:
        """
        Delete a range of amino acids from the sequence

        Parameters
        ----------
        start_idx (int): Start index (inclusive)
        end_idx (int): End index (exclusive)

        Returns
        -------
        SequenceMutation: Mutation information
        """
        if not (
            0 <= start_idx < len(self.original_sequence)
            and 0 <= end_idx <= len(self.original_sequence)
            and start_idx < end_idx
        ):
            raise ValueError(
                f"Invalid indices: start_idx={start_idx}, end_idx={end_idx}"
            )

        deleted_region = self.original_sequence[start_idx:end_idx]
        mutated_sequence = (
            self.original_sequence[:start_idx] + self.original_sequence[end_idx:]
        )

        return SequenceMutation(
            original_sequence=self.original_sequence,
            mutated_sequence=mutated_sequence,
            mutation_type=MutationType.DELETION,
            start_idx=start_idx,
            end_idx=end_idx,
            mutated_region="",
            original_region=deleted_region,
        )

    def shuffle_range(self, start_idx: int, end_idx: int) -> SequenceMutation:
        """
        Shuffle amino acids within a range

        Parameters
        ----------
        start_idx (int): Start index (inclusive)
        end_idx (int): End index (exclusive)

        Returns
        -------
        SequenceMutation: Mutation information
        """
        if not (
            0 <= start_idx < len(self.original_sequence)
            and 0 <= end_idx <= len(self.original_sequence)
            and start_idx < end_idx
        ):
            raise ValueError(
                f"Invalid indices: start_idx={start_idx}, end_idx={end_idx}"
            )

        region = list(self.original_sequence[start_idx:end_idx])
        shuffled_region = deepcopy(region)
        while "".join(shuffled_region) == "".join(region):  # Ensure actual shuffle
            random.shuffle(shuffled_region)

        mutated_sequence = (
            self.original_sequence[:start_idx]
            + "".join(shuffled_region)
            + self.original_sequence[end_idx:]
        )

        return SequenceMutation(
            original_sequence=self.original_sequence,
            mutated_sequence=mutated_sequence,
            mutation_type=MutationType.SHUFFLE,
            start_idx=start_idx,
            end_idx=end_idx,
            mutated_region="".join(shuffled_region),
            original_region=self.original_sequence[start_idx:end_idx],
        )

    def random_aa_range(self, start_idx: int, end_idx: int) -> SequenceMutation:
        """
        Replace a range with random amino acids

        Parameters
        ----------
        start_idx (int): Start index (inclusive)
        end_idx (int): End index (exclusive)


        Returns
        -------
        SequenceMutation: Mutation information
        """
        if not (
            0 <= start_idx < len(self.original_sequence)
            and 0 <= end_idx <= len(self.original_sequence)
            and start_idx < end_idx
        ):
            raise ValueError(
                f"Invalid indices: start_idx={start_idx}, end_idx={end_idx}"
            )

        length = end_idx - start_idx
        random_region = "".join(random.choices(self.AMINO_ACIDS, k=length))

        mutated_sequence = (
            self.original_sequence[:start_idx]
            + random_region
            + self.original_sequence[end_idx:]
        )

        return SequenceMutation(
            original_sequence=self.original_sequence,
            mutated_sequence=mutated_sequence,
            mutation_type=MutationType.RANDOM_AA,
            start_idx=start_idx,
            end_idx=end_idx,
            mutated_region=random_region,
            original_region=self.original_sequence[start_idx:end_idx],
        )

    def alanine_scan_range(self, start_idx: int, end_idx: int) -> SequenceMutation:
        """
        Replace a range with alanines (A)

        Parameters
        ----------
        start_idx (int): Start index (inclusive)
        end_idx (int): End index (exclusive)


        Returns
        -------
        SequenceMutation: Mutation information
        """
        if not (
            0 <= start_idx < len(self.original_sequence)
            and 0 <= end_idx <= len(self.original_sequence)
            and start_idx < end_idx
        ):
            raise ValueError(
                f"Invalid indices: start_idx={start_idx}, end_idx={end_idx}"
            )

        length = end_idx - start_idx
        alanine_region = "A" * length

        mutated_sequence = (
            self.original_sequence[:start_idx]
            + alanine_region
            + self.original_sequence[end_idx:]
        )

        return SequenceMutation(
            original_sequence=self.original_sequence,
            mutated_sequence=mutated_sequence,
            mutation_type=MutationType.ALANINE_SCAN,
            start_idx=start_idx,
            end_idx=end_idx,
            mutated_region=alanine_region,
            original_region=self.original_sequence[start_idx:end_idx],
        )


def exhaustive_mutation_scan(
    sequence: str,
    mutation_type: MutationType,
    window_size: int = 30,
    stride: int | None = None,
) -> list[SequenceMutation]:
    """
    Perform exhaustive mutation scan on a protein sequence

    Parameters
    ----------
    sequence (str): Protein sequence to mutate
    mutation_type (MutationType): Type of mutation to perform
    window_size (int): Size of mutation window
    stride (int, optional): Stride for sliding window. If None, uses window_size

    Returns
    -------
    List[SequenceMutation]: List of all mutations
    """
    if stride is None:
        stride = window_size

    manipulator = ProteinSequenceManipulator(sequence)
    mutations = []

    mutation_func = {
        MutationType.DELETION: manipulator.delete_range,
        MutationType.SHUFFLE: manipulator.shuffle_range,
        MutationType.RANDOM_AA: manipulator.random_aa_range,
        MutationType.ALANINE_SCAN: manipulator.alanine_scan_range,
    }[mutation_type]

    for start_idx in range(1, len(sequence), stride):
        end_idx = min(start_idx + window_size, len(sequence))
        mutation = mutation_func(start_idx, end_idx)
        mutations.append(mutation)

    return mutations


# Example usage:
if __name__ == "__main__":
    # Example protein sequence
    example_sequence = "MKLPVRRLVTFVTTQAEGKLV"

    # Create manipulator
    manipulator = ProteinSequenceManipulator(example_sequence)

    # Test individual mutations
    deletion = manipulator.delete_range(3, 6)
    shuffle = manipulator.shuffle_range(3, 6)
    random_aa = manipulator.random_aa_range(3, 6)
    alanine_scan = manipulator.alanine_scan_range(3, 6)

    print("Original sequence:", example_sequence)
    print("After deletion:", deletion.mutated_sequence)
    print("After shuffle:", shuffle.mutated_sequence)
    print("After random AA:", random_aa.mutated_sequence)
    print("After alanine scan:", alanine_scan.mutated_sequence)

    # Test exhaustive scan
    mutations = exhaustive_mutation_scan(
        example_sequence, MutationType.SHUFFLE, window_size=5, stride=3
    )

    print("\nExhaustive scan results:")
    for i, mutation in enumerate(mutations):
        print(f"Mutation {i + 1}:")
        print(f"Region {mutation.start_idx}-{mutation.end_idx}:")
        print(f"Original: {mutation.original_region}")
        print(f"Mutated:  {mutation.mutated_region}")
        print(f"Full sequence: {mutation.mutated_sequence}")
