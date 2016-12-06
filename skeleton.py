import sys
sys.path.append('./khmer')
import khmer

"""
bloom filter
"""

def build_bloom_filter(sequences, k):
    hbts = khmer.Nodegraph(k, 1e6, 1) # TODO: determine what those two other params are (they surely control false positives) 
    for dna in sequences:
        hbts.consume(dna)
    return hbts 

def test_bf():
    dna = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAG"
    hbts = build_bloom_filter([dna], 6)
    assert hbts.get("AGCTTT") == 1
    assert hbts.get("GATGAG") == 0

test_bf()


"""
Proposition 1
"""

def iterate_kmers(sequences, k):
    for seq in sequences:
        for i in xrange(len(seq)-k+1):
            yield seq[i:i+k]

def right_extensions(kmer):
    for c in "ACTG":
        yield kmer[1:]+c

def compute_branching_nodes(sequences, k):
    bf = build_bloom_filter(sequences, k)

    V_prime_plus = set()
    for kmer in iterate_kmers(sequences, k):
        kmer_right_extensions = [r_e for r_e in right_extensions(kmer)]
        kmer_right_extensions_in_bf = filter(lambda x: bf.get(x), kmer_right_extensions)
        if len(kmer_right_extensions_in_bf) >= 2:
            V_prime_plus.add(kmer)

    W_plus = set()
    for kmer in iterate_kmers(sequences,k):
        for r_e in right_extensions(kmer):
            if r_e in V_prime_plus:
                W_plus.add(r_e)

    # TODO do V_prime_minus 

    return W_plus

def test_compute_branching_nodes():
    dna = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAG"
    k=6
    print sorted(compute_branching_nodes([dna], k))

test_compute_branching_nodes()
