import genomeToTranscriptMapper as gttm

def test(exons, strand, begin, end, expected_begin, expected_end):
    converter = gttm.GenomeToTranscriptMapper(exons, strand)
    converted = converter.convert_interval_genome_to_transcript(begin,end)
    #print(f"{begin} - {end} - {expected_begin} - {expected_end} - {converted}")
    assert expected_begin == converted[0]
    assert expected_end == converted[1]


def run_test():
    exons = [(1000, 2000), (3000, 4000), (5000, 6000)]
    # + strand
    #Single exon
    test(exons, "+", 1000, 2000, 0, 1000)
    test(exons, "+", 1000, 1001, 0, 1)
    test(exons, "+", 3000, 4000, 1000, 2000)
    test(exons, "+", 3100, 3200, 1100, 1200)
    #Multiple exon
    test(exons, "+", 1100, 3100, 100, 1100)
    test(exons, "+", 1000, 6000, 0, 3000)
    #Partly outside of exon
    test(exons, "+", 500, 2000, 0, 1000)
    test(exons, "+", 500, 2100, 0, 1000)

    # - strand
    #Single exon
    test(exons, "-", 5999, 5998, 0, 1)
    test(exons, "-", 1999, 999, 2000, 3000)
    test(exons, "-", 5999, 5998, 0, 1)
    test(exons, "-", 3999, 3998, 1000, 1001)
    test(exons, "-", 3999, 3997, 1000, 1002)
    test(exons, "-", 3899, 3799, 1100, 1200)
    #Multiple exon
    test(exons, "-", 5899, 3899, 100, 1100)
    test(exons, "-", 5999, 999, 0, 3000)
    #Partly outside of exon
    test(exons, "-", 6500, 4999, 0, 1000)
    test(exons, "-", 6500, 4399, 0, 1000)



run_test()