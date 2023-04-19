#!/usr/bin/env python3

import numpy as np
from typing import Tuple, List
from sys import stderr
from enum import Enum

""" Convert coordinates or intervals from genome space to transcript space.
    
    #__init__:
    #Initialize mapper object with the transcript informations:
        * exons_interval: list of tuples in the form (begin, end) representing exons
            * intervals are left inclusive and right exclusive.
            * order is not important
        * strand: "+" or "-"
    
    #convert_points_genome_to_transcript:
    #Convert a list of points from genome space to transcript space
    #Please note: if you have many coordinates to convert, calling the method once with
    # a list of points is way faster than calling it many times with single values.
        * positions: list of integers representing genome space positions.
        * outside_exons_points [optional]: allow to specify the behaviour for points outside of transcript:
            * OutOfExonsPoint.THROW_EXCEPTION [default behaviour]
            * OutOfExonsPoint.RETURN_MINUS_ONE
            * OutOfExonsPoint.RESTRICT_SMALLER: return first transcript position on the left of the specified one.
            * OutOfExonsPoint.RESTRICT_BIGGER: return first transcript position on the right of the specified one.

    #convert_interval_genome_to_transcript:
    #Convert an interval from genome space to transcript space.  
        * begin, end: integers.
        * restrict_interval_to_transcript: specify what to do if the interval is only partially
        *   included in the transcript:
            * True: interval is cut to fit inside the transcript.
            * False: if interval isn't completely included in the transcript, an execption is raised.
    #NOTE: all intervals are considered left side inclusive and right side exclusive
    #For example, for strand + : (1000, 2000) => from 1000 to 1999 included
    #             for strand - : (1999, 999) => from 1999 to 1000 included
                      strand - alternative: (1000, 2000) => again from 1999 to 1000
"""

#Enum class to express how to treat points outside of exons
class OutOfExonsPoint(Enum):
    THROW_EXCEPTION = 1
    RETURN_MINUS_ONE = 2
    RESTRICT_SMALLER = 3
    RESTRICT_BIGGER = 4

class GenomeToTranscriptMapper:

    #Initialize mapper object with the transcript exons and strand
    #Intervals are in the form [b .. e)
    def __init__(self, exons_interval: Tuple[int, int], strand: str):
        assert strand == '+' or strand == '-'

        self.strand = strand
        # sort the exons into "transcript" order
        exons_interval.sort()
        #Check that intervals aren't overlapped
        last_e = -1
        for exon in exons_interval:
            assert exon[0] > last_e
            last_e = exon[1]
        #Reverse intervals for negative strand
        if self.strand == "-":
            exons_interval.sort(reverse=True)
            exons_interval = [(y-1, x-1) for x, y in exons_interval]
        
        #Offset: coordinate of the start of first exon for strand +, end of last exon for strand -
        offset = 0
        intervals_with_offset = []

        #For each exon we store:
        # start, end -> In DNA coordinates
        # offset -> start of the exon in transcript coordinates
        for interval in exons_interval:
            intervals_with_offset.append(GenomeToTranscriptMapper.IntervalWithOffset(offset = offset, start = interval[0], end = interval[1]))
            offset = offset + abs(interval[1]-interval[0])

        self.intervals_with_offset = np.array(intervals_with_offset)

    #This function converts a list of points from genome space to transcript space
    #It raises ValueError if one of the points is outside of the transcript
    def convert_points_genome_to_transcript(self, positions: List[int], outside_exons_points: OutOfExonsPoint = OutOfExonsPoint.THROW_EXCEPTION):
        if (len(positions)==0):
            return []
        #Sort pos and keeps track of original position of elements
        if (self.strand=='-'):
            positions = np.array(positions)
            ordering = np.argsort(-positions)
            pos = -np.sort(-positions)
        else:
            ordering = np.argsort(positions)
            pos = np.sort(positions)
        
        results = np.zeros(len(pos))
        index = 0
        exon_index = 0
        while(exon_index < len(self.intervals_with_offset)):
            interval = self.intervals_with_offset[exon_index]
            if (index == len(pos)):
                break

            position = pos[index]
            if (self.strand=="-"):
                if (position <= interval.begin and position > interval.end):
                    results[index] = interval.offset + (interval.begin - position)
                    index += 1
                else:
                    if (outside_exons_points == OutOfExonsPoint.RESTRICT_SMALLER
                        and position > interval.begin):
                        results[index] = interval.offset
                        index += 1 
                    elif (outside_exons_points == OutOfExonsPoint.RESTRICT_BIGGER and (
                        position <= interval.end and (exon_index == len(self.intervals_with_offset)-1 or position > self.intervals_with_offset[exon_index+1].begin)
                    )):
                        results[index] = interval.offset + interval.begin - interval.end - 1
                        index += 1 
                    else:
                        exon_index += 1
            else:
                if (position >= interval.begin and position < interval.end):
                    results[index] = interval.offset + (position - interval.begin)
                    index += 1
                else:
                    if (outside_exons_points == OutOfExonsPoint.RESTRICT_BIGGER
                        and position < interval.begin):
                        results[index] = interval.offset
                        index += 1 
                    elif (outside_exons_points == OutOfExonsPoint.RESTRICT_SMALLER and (
                        position >= interval.end and (exon_index == len(self.intervals_with_offset)-1 or position < self.intervals_with_offset[exon_index+1].begin)
                    )):
                        results[index] = interval.offset + interval.end - interval.begin - 1
                        index += 1 
                    else:
                        exon_index += 1

        if (index == len(pos)):
            return results[ordering]

        if (outside_exons_points == OutOfExonsPoint.RETURN_MINUS_ONE):
            for index in range(index, len(pos)):
                pos[index] = -1 
            return pos
        else:
            raise ValueError(f"Transcript position {pos[index]} outside of transcript")


    def convert_interval_genome_to_transcript(self, begin: int, end: int, restrict_interval_to_transcript = True):
        #Makes interval with both sides inclusive - [x..y]
        if (self.strand=='-' and begin > end):
            interval = [end + 1, begin]
        else:
            interval = [begin, end-1]
        #Convert coordinates
        if (restrict_interval_to_transcript):
            interval[0] = self.convert_points_genome_to_transcript((interval[0],), OutOfExonsPoint.RESTRICT_BIGGER)[0]
            interval[1] = self.convert_points_genome_to_transcript((interval[1],), OutOfExonsPoint.RESTRICT_SMALLER)[0]
            if ((self.strand == "+" and interval[1] < interval[0]) or (self.strand == "-" and interval[1] > interval[0])):
                raise ValueError(f"Interval in {begin} - {end} is included inside an intron.")
        else:
            interval = self.convert_points_genome_to_transcript(interval, OutOfExonsPoint.THROW_EXCEPTION)
        #Makes interval again in the form [x..y)
        interval = sorted(interval)
        interval[1] += 1
        return interval

    #Inner class to represent intervals with offsets
    class IntervalWithOffset:
        def __init__(self, start: int, end: int, offset: int):
            #Offset is the distance between the start of the first interval and the start of this interval
            self.offset = offset;
            self.begin = start;
            self.end = end;
        
        def __str__(self):
            return f"{self.begin} - {self.end}"