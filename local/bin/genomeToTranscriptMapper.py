#!/usr/bin/env python3

import numpy as np
from typing import Tuple, List
from sys import stderr


class GenomeToTranscriptMapper:

    #Inner class to represent intervals with offsets
    class IntervalWithOffset:
        def __init__(self, start: int, end: int, offset: int):
            self.begin = start; self.end = end; self.offset = offset;
        
        def __str__(self):
            return f"{self.begin} - {self.end}"

    #Initialize mapper object with the transcript exons and strand
    #Intervals are in the form [b .. e)
    def __init__(self, exons_interval: Tuple[int, int], strand: str):
        assert strand == '+' or strand == '-'

        self.strand = strand
        # sort the exons into "transcript" order
        if self.strand == "-":
            exons_interval.sort(reverse=True)
            exons_interval = [(y-1, x-1) for x, y in exons_interval]
        else:
            exons_interval.sort(reverse=False)
        
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
    def convert_points_genome_to_transcript(self, pos: List[int]):
        if (len(pos)==0):
            return []

        ordering = np.argsort(pos)
        pos = np.sort(pos)
        
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
                    exon_index += 1
            else:
                if (position >= interval.begin and position < interval.end):
                    results[index] = interval.offset + (position - interval.begin)
                    index += 1
                else:
                    print(f"Pos {position} not in {interval}", file=stderr)
                    exon_index += 1

        if (index == len(pos)):
            return results[ordering]

        print(f"Position {pos[index]} - strand {self.strand} - min {self.intervals_with_offset[0]} - max {self.intervals_with_offset[len(self.intervals_with_offset)-1]}" ,file=stderr)
        raise ValueError(f"Transcript postion {pos[index]} outside of transcript")

    def convert_interval_genome_to_transcript(self, begin: int, end: int):
        #Makes interval with both sides inclusive - [x..y]
        if (self.strand=='-'):
            interval = (end - 1, begin)
        else:
            interval = (begin, end-1)

        transcript_list = self.convert_points_genome_to_transcript(interval)
        #Makes interval again in the form [x..y)
        transcript_list[1] += 1
        return sorted(transcript_list)