"""
Classes for peak alignment by dynamic programming.
"""

################################################################################
#                                                                              #
#    PyMassSpec software for processing of mass-spectrometry data              #
#    Copyright (C) 2005-2012 Vladimir Likic                                    #
#    Copyright (C) 2019-2020 Dominic Davis-Foster                              #
#                                                                              #
#    This program is free software; you can redistribute it and/or modify      #
#    it under the terms of the GNU General Public License version 2 as         #
#    published by the Free Software Foundation.                                #
#                                                                              #
#    This program is distributed in the hope that it will be useful,           #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of            #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
#    GNU General Public License for more details.                              #
#                                                                              #
#    You should have received a copy of the GNU General Public License         #
#    along with this program; if not, write to the Free Software               #
#    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                 #
#                                                                              #
################################################################################

# stdlib
import copy
import functools
import math
from typing import Dict, List
import multiprocessing
import time
import pathlib
import logging

# 3rd party
import numpy  # type: ignore

# this package
from pyms.Peak import Peak
import gcms_align.settings

# try:
#     # 3rd party
#     from mpi4py import MPI  # type: ignore
# except ModuleNotFoundError:
#     pass

# this package
from pyms.DPA.Alignment import Alignment
from pyms.DPA.clustering import treecluster
from pyms.Utils.Utils import is_sequence_of

__all__ = [
    "PairwiseAlignment",
    "align",
    "score_matrix",
    "dp",
    "position_similarity",
    "merge_alignments",
    "alignment_similarity",
    "alignment_compare",
    # "score_matrix_mpi",
    "align_with_tree",  # "align_with_tree_mpi",
]


class PairwiseAlignment:
    """
    Models pairwise alignment of alignments.

    :param alignments: A list of alignments.
    :param D: Retention time tolerance parameter for pairwise alignments.
    :param gap: Gap parameter for pairwise alignments.

    :authors: Woon Wai Keen, Vladimir Likic
    """

    def __init__(self, alignments: List[Alignment], D: float, gap: float, pairwise_temp_file: pathlib.Path, config_settings):
        if not is_sequence_of(alignments, Alignment):
            raise TypeError("'alignments' must be a Sequence of Alignment objects")

        if not isinstance(D, float):
            raise TypeError("'D' must be a float")

        if not isinstance(gap, float):
            raise TypeError("'gap' must be a float")

        self.alignments = alignments
        self.D = D
        self.gap = gap
        self._pairwise_temp_file = pairwise_temp_file
        self.config = config_settings

        self._sim_matrix()
        self._dist_matrix()
        self._guide_tree()

    def _sim_matrix(self) -> None:
        """
        Calculates the similarity matrix for the set of alignments.

        :authors: Woon Wai Keen, Vladimir Likic
        """

        n = len(self.alignments)

        total_n = n * (n - 1) // 2

        print(f" Calculating pairwise alignments for {n:d} alignments (D={self.D:.2f}, gap={self.gap:.2f})")

        self.sim_matrix = numpy.zeros((n, n), dtype='f')

        if self.config.align_multiprocess:
            align_source = MultiAlign(self.alignments, self.D, self.gap, self._pairwise_temp_file, self.config)
            align_source.process_tasks(self.sim_matrix)
            
            align_source.process_terminate()
        else:
            old_sim_matrix = numpy.zeros((n, n), dtype='f')
            for i in range(n - 1):
                for j in range(i + 1, n):
                    ma = align(self.alignments[i], self.alignments[j], self.D, self.gap)
                    old_sim_matrix[i, j] = old_sim_matrix[j, i] = ma.similarity
                    total_n = total_n - 1
                    print(f" -> {total_n:d} pairs remaining")
            print(f"{numpy.all(self.sim_matrix == old_sim_matrix)}")

    def _dist_matrix(self) -> None:
        """
        Converts similarity matrix into a distance matrix.

        :authors: Woon Wai Keen, Vladimir Likic
        """

        # change similarity matrix entries (i,j) to max{matrix}-(i,j)
        sim_max = numpy.max(numpy.ravel(self.sim_matrix))
        self.dist_matrix = sim_max - self.sim_matrix

        # set diagonal elements of the similarity matrix to zero
        for i in range(len(self.dist_matrix)):
            self.dist_matrix[i, i] = 0

    def _guide_tree(self) -> None:
        """
        Build a guide tree from the distance matrix.

        :authors: Woon Wai Keen, Vladimir Likic
        """

        n = len(self.dist_matrix)

        print(f" -> Clustering {n * (n - 1):d} pairwise alignments.", end='')
        self.tree = treecluster(data=None, distancematrix=self.dist_matrix, method='a')
        print("Done")


def align(a1: Alignment, a2: Alignment, D: float, gap: float) -> Alignment:
    """
    Aligns two alignments.

    :param a1: The first alignment
    :param a2: The second alignment
    :param D: Retention time tolerance
    :param gap: Gap penalty

    :return: Aligned alignments

    :authors: Woon Wai Keen, Vladimir Likic
    """

    # comm = MPI.COMM_WORLD
    # rank = comm.Get_rank()

    # calculate score matrix for two alignments
    M = score_matrix(a1, a2, D)
    # print("calculated score matrix on rank", rank)

    # run dynamic programming
    result_trace = dp(M, gap)

    # make composite alignment from the results
    ma = merge_alignments(a1, a2, result_trace)

    # calculate the similarity score
    ma.similarity = alignment_similarity(result_trace, M, gap)

    return ma


def score_matrix(a1: Alignment, a2: Alignment, rt_sensitivity: float) -> numpy.ndarray:
    """
    Calculates the score matrix between two alignments.

    :param a1: The first alignment.
    :param a2: The second alignment.
    :param D: Retention time tolerance.

    :return: Aligned alignments.

    :authors: Qiao Wang, Andrew Isaac
    """

    # sim_score = 0
    # Attempt to speed up by only calculating 'in-range' values
    # set tolerance to 1/1000
    _TOL = 0.001
    cutoff = rt_sensitivity * math.sqrt(-2.0 * math.log(_TOL))

    score_matrix = numpy.ones((len(a1.peakalgt), len(a2.peakalgt)), numpy.double)

    for i, algt1pos in enumerate(a1.peakalgt):
        for j, algt2pos in enumerate(a2.peakalgt):
            sim_score = position_similarity(algt1pos, algt2pos, rt_sensitivity, cutoff)
            score_matrix[i][j] = sim_score

    return score_matrix


def compressed_score_matrix(a1: Alignment, a2: Alignment, settings):
    """
    Calculates a partial score matrix between two alignments.

    :param a1: The first alignment.
    :param a2: The second alignment.
    :param D: Retention time tolerance.

    :return: Alignments and positions 
    """

    # set tolerance to 1/1000
    _TOL = 0.001
    rt_sensitivity = settings.rt_sensitivity_s
    start_backtrack = settings.align_score_backtrack
    end_check = settings.align_end_row_check
    max_row_width = end_check*2 + start_backtrack*2
    cutoff = rt_sensitivity * math.sqrt(-2.0 * math.log(_TOL))

    score_matrix = numpy.ones((len(a1.peakalgt), max_row_width+1), numpy.double)
    offset_list = numpy.zeros(len(a1.peakalgt), numpy.int)
    start_index = 0

    for i, algt1pos in enumerate(a1.peakalgt):
        failed_count = 0
        adjust_start = True
        if start_index >= start_backtrack:
            start_index -= start_backtrack

        for j in range(len(a2.peakalgt)):
            algt2pos = a2.peakalgt[start_index+j]
            sim_score = position_similarity(algt1pos, algt2pos, rt_sensitivity, cutoff)
            if sim_score > 0.999:
                if adjust_start:
                    if j-start_index == start_backtrack:
                        start_index = j
                else:
                    failed_count += 1
            else:
                adjust_start = False
                failed_count = 0
            if failed_count > max_row_width or not adjust_start and failed_count > end_check:
                break
            # else: deal with no nonzero entries in the row
            score_matrix[i][j] = sim_score
            offset_list[i] = start_index

    return score_matrix, offset_list


def dp(S, gap_penalty: float) -> Dict:
    """
    Solves optimal path in score matrix based on global sequence alignment.

    :param S: Score matrix
    :param gap_penalty: Gap penalty

    :return: A dictionary of results

    :author: Tim Erwin
    """

    try:
        row_length = len(S[:, 0])
    except IndexError:
        raise IndexError("Zero length alignment found: Samples with no peaks cannot be aligned")

    col_length = len(S[0, :])

    # D contains the score of the optimal alignment
    D = numpy.zeros((row_length + 1, col_length + 1), dtype='d')

    for i in range(1, row_length + 1):
        D[i, 0] = gap_penalty * i
    for j in range(1, col_length + 1):
        D[0, j] = gap_penalty * j

    D[0, 0] = 0.0
    D[1:(row_length + 1), 1:(col_length + 1)] = S

    # Directions for trace
    # 0 - match               (move diagonal)
    # 1 - peaks1 has no match (move up)
    # 2 - peaks2 has no match (move left)
    # 3 - stop
    trace_matrix = numpy.zeros((row_length + 1, col_length + 1))
    trace_matrix[:, 0] = 1
    trace_matrix[0, :] = 2
    trace_matrix[0, 0] = 3

    # TODO: Calculate entries only as required (lazy)
    for i in range(1, row_length + 1):
        for j in range(1, col_length + 1):
            #
            # Needleman-Wunsch Algorithm assuming a score function S(x,x)=0
            #
            #              | D[i-1,j-1] + S(i,j)
            # D[i,j] = min | D(i-1,j] + gap
            #              | D[i,j-1] + gap
            #

            darray = [D[i - 1, j - 1] + S[i - 1, j - 1], D[i - 1, j] + gap_penalty, D[i, j - 1] + gap_penalty]
            D[i, j] = min(darray)
            # Store direction in trace matrix
            trace_matrix[i, j] = darray.index(D[i, j])

    # Trace back from bottom right
    trace = []
    # matches = []
    i = row_length
    j = col_length
    direction = trace_matrix[i, j]
    # p = [row_length - 1]
    # q = [col_length - 1]

    while direction != 3:

        if direction == 0:  # Match
            i = i - 1
            j = j - 1
            # matches.append([i, j])
        elif direction == 1:  # peaks1 has no match
            i = i - 1
        elif direction == 2:  # peaks2 has no match
            j = j - 1
        # p.append(i - 1)
        # q.append(j - 1)
        trace.append(direction)
        direction = trace_matrix[i, j]

        if False:
            direction = get_direction(D, S, i, j, gap_penalty)

    # remove 'stop' entry
    # p.pop()
    # q.pop()
    # reverse the trace back
    # p.reverse()
    # q.reverse()
    trace.reverse()
    # matches.reverse()

    return trace  # 'p': p, 'q': q, , "matches": matches, 'D': D, "phi": trace_matrix}


def d_max(D, S):
    # worst case 1: no matches all gap penalties
    numpy.indices()
    # calculate main diagonal, all matches successful


def get_direction(D, S, i, j, gap_penalty):
    """
    """
    darray = [D[i - 1, j - 1] + S[i - 1, j - 1], D[i - 1, j] + gap_penalty, D[i, j - 1] + gap_penalty]
    D[i, j] = min(darray)
    return darray.index(D[i, j])


def get_D(D, S, i, j, gap_penalty):
    """
    """
    match_direction_penalty = get_D(D, S, i - 1, j - 1, gap_penalty) + S[i - 1, j - 1]
    pass


def position_similarity(pos1, pos2, rt_sens, cutoff) -> float:
    """
    Calculates the similarity between the two alignment positions.

    A score of 0 is best and 1 is worst.

    :param pos1: The position of the first alignment.
    :param pos2: The position of the second alignment.
    :param D: Retention time tolerance.

    :return: The similarity value for the current position.

    :authors: Qiao Wang, Vladimir Likic, Andrew Isaac
    """

    score = 0.0
    count = 0

    for a in pos1:
        if a is not None:
            mass_spect1 = a._mass_spectrum.mass_spec  # TODO: Update method of direct access avoiding copy operation
            art = a.rt
            once = True

            for b in pos2:
                if b is not None:
                    brt = b.rt
                    # in range?

                    if abs(art - brt) > cutoff:
                        score += 1.0  # NB score of 1 is worst
                    else:
                        # Once per b-loop
                        if once:
                            mass_spect1_sum = numpy.sum(mass_spect1**2, axis=0)
                            once = False


                        # TODO: Update mass spec to only have != 0 intensities for sum**2
                        mass_spect2 = b._mass_spectrum.mass_spec
                        mass_spect2_sum = numpy.sum(mass_spect2**2, axis=0)

                        try:
                            top = numpy.dot(mass_spect1, mass_spect2)
                        except ValueError:
                            raise ValueError(
                                """Mass Spectra are of different lengths.
Use `IntensityMatrix.crop_mass()` to set same length for all Mass Spectra"""
                            )
                        all_squared = mass_spect1_sum * mass_spect2_sum
                        if all_squared > 0:
                            cos = top / numpy.sqrt(all_squared)
                            rtime = numpy.exp(-((art - brt) / float(rt_sens))**2 / 2.0)
                            score = score + (1.0 - (cos * rtime))
                        else:
                            score = score + 1.0
                    count += 1

    if count == 0:
        score = 1.0  # NB score of 1 is worst
    else:
        score = score / float(count)

    return score


def merge_alignments(A1: Alignment, A2: Alignment, traces) -> Alignment:
    """
    Merges two alignments with gaps added in from DP traceback.

    :param A1: First alignment.
    :param A2: Second alignment.
    :param traces: DP traceback.

    :return: A single alignment from ``A1`` and ``A2``.

    :authors: Woon Wai Keen, Vladimir Likic, Qiao Wang
    """

    # Create object to hold new merged alignment and fill in its expr_codes
    ma = Alignment(None)
    ma.expr_code = A1.expr_code + A2.expr_code

    # create empty lists of dimension |A1| + |A2|
    dimension = len(A1.peakpos) + len(A2.peakpos)
    merged: List[List[Peak]] = [[] for _ in range(dimension)]

    idx1 = idx2 = 0

    # trace can either be 0, 1, or 2
    # if it is 0, there are no gaps. otherwise, if it is 1 or 2,
    # there is a gap in A2 or A1 respectively.

    for trace in traces:
        if trace in {0, 1}:
            for i, _ in enumerate(A1.peakpos):
                merged[i].append(A1.peakpos[i][idx1])
            idx1 = idx1 + 1

        elif trace == 2:
            for i, _ in enumerate(A1.peakpos):
                merged[i].append(None)  # type: ignore
        # ---

        if trace in {0, 2}:
            for j, peak in enumerate(A2.peakpos):
                merged[1 + i + j].append(peak[idx2])
            idx2 = idx2 + 1

        elif trace == 1:
            for j, _ in enumerate(A2.peakpos):
                merged[1 + i + j].append(None)  # type: ignore

    ma.peakalgt = numpy.transpose(merged)
    # sort according to average peak
    ma.peakalgt = list(ma.peakalgt)
    ma.peakalgt.sort(key=functools.cmp_to_key(alignment_compare))
    ma.peakpos = numpy.transpose(ma.peakalgt)

    return ma


def alignment_similarity(traces, score_matrix, gap: float) -> float:
    """
    Calculates similarity score between two alignments (new method).

    :param traces: Traceback from DP algorithm.
    :param score_matrix: Score matrix of the two alignments.
    :param gap: Gap penalty.

    :return: Similarity score (i.e. more similar => higher score)

    :authors: Woon Wai Keen, Vladimir Likic
    """

    score_matrix = 1.0 - score_matrix
    similarity = 0.
    idx1 = idx2 = 0

    # Trace can either be 0, 1, or 2
    # If it is 0, there is a match and we add to the sum the score between
    # these two aligned peaks.
    #
    # Otherwise, if it is 1 or 2, and there is a gap in A2 or A1
    # respectively. We then subtract the gap penalty from the sum.
    for trace in traces:
        if trace == 0:
            similarity = similarity + score_matrix[idx1][idx2]
            idx1 = idx1 + 1
            idx2 = idx2 + 1
        elif trace == 1:
            similarity = similarity - gap
            idx1 = idx1 + 1
        elif trace == 2:
            similarity = similarity - gap
            idx2 = idx2 + 1

    return similarity


def alignment_compare(x, y):
    """
    A helper function for sorting peak positions in a alignment.

    :param x:
    :param y:
    """

    x = [_.rt for _ in filter(None, x)]
    y = [_.rt for _ in filter(None, y)]

    avg_x = numpy.sum(x) / len(x)
    avg_y = numpy.sum(y) / len(y)

    if avg_x < avg_y:
        return -1
    else:
        return 1


# def score_matrix_mpi(a1: Alignment, a2: Alignment, D: float) -> Alignment:
#     """
#     Calculates the score matrix between two alignments.
#
#     :param a1: The first alignment.
#     :param a2: The second alignment.
#     :param D: Retention time tolerance.
#
#     :return: Aligned alignments
#
#     :authors: Qiao Wang, Andrew Isaac
#     """
#
#     # sim_score = 0
#
#     comm = MPI.COMM_WORLD
#     rank = comm.Get_rank()
#     size = comm.Get_size()
#
#     portion = int(float(len(a1.peakalgt)) / size)
#
#     # print("length of a1.peakalgt =", len(a1.peakalgt))
#     # print("portion size = ", portion)
#
#     # if rank == 0:
#     score_matrix = numpy.zeros((len(a1.peakalgt), len(a2.peakalgt)))
#
#     if rank < size - 1:  # if it's not the last slice
#         score_matrix_part = numpy.zeros((portion, len(a2.peakalgt)))
#         a1_part = a1.peakalgt[rank * portion:(rank + 1) * portion]
#     else:  # if it's the last strip, prob not full portion
#         score_matrix_part = numpy.zeros((
#             len(a1.peakalgt) - (rank * portion),
#             len(a2.peakalgt),
#         ))
#         a1_part = a1.peakalgt[rank * portion:len(a1.peakalgt)]
#
#     for i, algt1pos in enumerate(a1_part):
#         for j, algt2pos in enumerate(a2.peakalgt):
#             sim_score = position_similarity(algt1pos, algt2pos, D)
#             score_matrix_part[i][j] = sim_score
#
#     if rank == 0:
#         score_matrix[0:portion] = score_matrix_part
#         for i in range(1, size):
#             if i == size - 1:
#                 recv_buffer = numpy.zeros((len(a1.peakalgt) - (i * portion), len(a2.peakalgt)))
#                 comm.Recv(recv_buffer, i)
#                 score_matrix[i * portion:len(a1.peakalgt)] = recv_buffer
#             else:
#                 recv_buffer = numpy.zeros((portion, len(a2.peakalgt)))
#                 comm.Recv(recv_buffer, i)
#                 score_matrix[i * portion:(i + 1) * portion] = recv_buffer
#
#     # print("I am rank 0, done!", score_matrix)
#
#     else:
#         # all other process send their result
#         comm.Send(score_matrix_part)
#
#     outputs = []
#     if rank == 0:
#         for rank in range(size):
#             # print rank
#             outputs.append(score_matrix)
#
#     score_matrix = comm.scatter(outputs, root=0)
#
#     # print("before return, rank", rank)
#     return score_matrix


def align_with_tree(T: PairwiseAlignment, min_peaks: int = 1) -> Alignment:
    """
    Aligns a list of alignments using the supplied guide tree.

    :param T: The pairwise alignment object.
    :param min_peaks:

    :return: The final alignment consisting of aligned input alignments.

    :authors: Woon Wai Keen, Vladimir Likic
    """

    print(f" Aligning {len(T.alignments):d} items with guide tree (D={T.D:.2f}, gap={T.gap:.2f})")

    # For everything else, we align according to the guide tree provided by
    # Pycluster. From Pycluster documentation:
    #   Each item and subnode is represented by an integer. For hierarchical
    #   clustering of n items, we number the original items {0, ... , n-1},
    #   nodes are numbered {-1, ... , -(n-1)}. Note that the number of nodes
    #   is one less than the number of items.

    # extend As to length 2n to hold the n items, n-1 nodes, and 1 root
    As: List[Alignment] = copy.deepcopy(T.alignments) + [None for _ in range(len(T.alignments))]  # type: ignore

    # align the alignments into positions -1, ... ,-(n-1)
    total = len(T.tree)
    index = 0

    for node in T.tree[:]:
        index = index - 1
        As[index] = align(As[node.left], As[node.right], T.D, T.gap)
        total = total - 1
        print(f" -> {total:d} item(s) remaining")

    # the final alignment is in the root. Filter min peaks and return
    final_algt: Alignment = As[index]

    # useful for within state alignment only
    if min_peaks > 1:
        final_algt.filter_min_peaks(min_peaks)

    return final_algt


#
# def align_with_tree_mpi(T: Alignment, min_peaks: int = 1) -> Alignment:
#     """
#     Aligns a list of alignments using the supplied guide tree
#
#     :param T: The pairwise alignment object
#     :param min_peaks:
#
#     :return: The final alignment consisting of aligned input alignments
#
#     :authors: Woon Wai Keen, Vladimir Likic
#     """
#
#     try:
#         rank = MPI.COMM_WORLD.Get_rank()
#     except:
#         rank = 0
#
#     if rank == 0:
#         print(f" Aligning {len(T.alignments):d} items with guide tree (D={T.D:.2f}, gap={T.gap:.2f})")
#
#     # For everything else, we align according to the guide tree provided by
#     # Pycluster. From Pycluster documentation:
#     #   Each item and subnode is represented by an integer. For hierarchical
#     #   clustering of n items, we number the original items {0, ... , n-1},
#     #   nodes are numbered {-1, ... , -(n-1)}. Note that the number of nodes
#     #   is one less than the number of items.
#
#     # extend As to length 2n to hold the n items, n-1 nodes, and 1 root
#     As = copy.deepcopy(T.alignments) + [None for _ in range(len(T.alignments))]
#
#     # align the alignments into positions -1, ... ,-(n-1)
#     total = len(T.tree)
#     index = 0
#
#     for node in T.tree[:]:
#         index = index - 1
#         As[index] = align(As[node.left], As[node.right], T.D, T.gap)
#         total = total - 1
#         if rank == 0:
#             print(f" -> {total:d} item(s) remaining")
#
#     # the final alignment is in the root. Filter min peaks and return
#     final_algt = As[index]
#
#     # useful for within state alignment only
#     if min_peaks > 1:
#         final_algt.filter_min_peaks(min_peaks)
#
#     return final_algt


# import getopt, math, random, sys, time, types, wx
#
# from multiprocessing import Process, Queue, cpu_count, current_process, freeze_support


class MultiAlign(object):
    """ """
    def __init__(self, alignments, D, gap, temp_align_file: pathlib.Path, config: gcms_align.settings.Setting):
        numproc = multiprocessing.cpu_count()-3
        self.temp_align_file = temp_align_file
        self.task_queue = multiprocessing.Queue()
        self.done_queue = multiprocessing.Queue()
        self.processes = []
        for _ in range(max(2, numproc)):
            process = multiprocessing.Process(target=MultiAlign.worker, args=(self.task_queue, self.done_queue, D, gap))
            process.start()
            self.processes.append(process)

        self.numprocesses = len(self.processes)

        n = len(alignments)
        self.sim_matrix = numpy.zeros((n, n), dtype='f')

        if temp_align_file.is_file():
            prior_data = self.read_prior_alignment()
        else:
            prior_data = {}

        tasks = []

        if config.align_sparse_mode:
            for i in range(n - 1):
                for j in range(i + 1, n):
                    if j-i <= config.align_diagonal_width or i in config.align_full_compare or j in config.align_full_compare:
                        if (i, j) not in prior_data and (j, i) not in prior_data:
                            tasks.append((alignments[i], alignments[j], (i, j)))
        else:
            for i in range(n - 1):
                for j in range(i + 1, n):
                    if (i, j) not in prior_data and (j, i) not in prior_data:
                        tasks.append((alignments[i], alignments[j], (i, j)))

        self._prior_results = prior_data
        self.tasks = tasks
        self.numtasks = len(tasks)
        self.keepgoing = True
        self.task_queue_index = 0
        self.done_queue_index = 0
        print(f" Setup for parallel mode: {numproc}")

    def read_prior_alignment(self):
        result_dict = {}
        with open(self.temp_align_file, "r") as f:
            line_list = f.readlines()
        for line in line_list:
            items = line.split(":")
            if len(items) == 3:
                result_dict[(int(items[0]), int(items[1]))] = float(items[2])
            elif len(items) > 0:
                logging.warning(f"Ignoring: {line}")
        return result_dict

    def process_tasks(self, results_matrix: numpy.array):
        """
        Start the execution of tasks by the processes.
        """
        self.keepgoing = True
        for (index0, index1), similarity in self._prior_results.items():
            results_matrix[index0, index1] = results_matrix[index1, index0] = similarity

        # Submit first set of tasks
        numprocstart = min(self.numprocesses, self.numtasks)
        for self.task_queue_index in range(numprocstart):
            self.task_queue.put(self.tasks[self.task_queue_index])

        self.done_queue_index = -1  # done queue index
        self.task_queue_index = numprocstart - 1  # task queue index
        while (self.done_queue_index < self.task_queue_index):
            # Get and print results
            self.done_queue_index += 1
            output = self.done_queue.get()
            _pid, index_set, ma_similarity = output
            results_matrix[index_set[0], index_set[1]] = results_matrix[index_set[1], index_set[0]] = ma_similarity
            print(f" -> Completed {self.done_queue_index:d} of {len(self.tasks)} pairs")
            with open(self.temp_align_file, "a") as f:
                f.write(f"{index_set[0]} : {index_set[1]} : {ma_similarity} \n")

            if ((self.keepgoing) and (self.task_queue_index + 1 < self.numtasks)):
                # Submit another task
                self.task_queue_index += 1
                self.task_queue.put(self.tasks[self.task_queue_index])

    def process_terminate(self):
        """
        Kill processes.
        """
        self.keepgoing = False
        for n in range(self.numprocesses):
            # Terminate any running processes
            self.processes[n].terminate()

        # Wait for all processes to terminate
        isalive = True
        time_count = 0
        while isalive and time_count < 40:
            time.sleep(0.25)
            isalive = numpy.any([x.is_alive() for x in self.processes])
            time_count += 1

    @classmethod
    def worker(cls, input_queue, output_queue, D, gap):
        """
        Calculate align results from the input_queue.
        """
        while True:
            args = input_queue.get()
            ma = align(args[0], args[1], D, gap)
            output_queue.put((multiprocessing.current_process().pid, args[2], ma.similarity))
