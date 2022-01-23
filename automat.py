from mml import *
import copy

PREFIX = "ggatgt"
A = "ctggct"
B = "cgcagc"
T = "tgtcgc"

CLEAVE_DISTANCE = (9, 13)
IS_ODD_SEQUENCE = "tgtcgc"

TRANSITIONS = [
    DNA("ggatgtag",   "cctacatgccga"),   # S0 -a-> S0
    DNA("ggatgacg",   "cctactgcgacc"),   # S1 -a-> S1
    DNA("ggatgacgac", "cctactgctggtcg"), # S0 -b-> S1
    DNA("ggatgg",     "cctaccgcgt"),     # S1 -b-> S0

]

CLEAVE_VALUE = ("ggatg", DNA.get_pair_string("ggatg"))


def automat_machine_even_b(input_sequence):
    # Generate the lab molucules
    print("[*] Generating molucule lab")
    transitions = copy.deepcopy(TRANSITIONS)
    for i in range(4):
        transitions = copy.deepcopy(transitions) + copy.deepcopy(transitions)
    lab = Lab(transitions)

    print("[*] Adding input dna sequence")
    dna = DNA(input_sequence, "", fill=True)
    lab.samples.append(dna)

    # Run while we have the cleaving sequence
    while True:
        # mimick Fok1 cleave
        trash, dna = dna.cleave(CLEAVE_VALUE, CLEAVE_DISTANCE)
        lab.samples.remove(trash)
        lab.samples.append(dna)
        lab.merge()

        # Failed merging our dna, that means the run ended
        if dna.removed == False:
            break

        # Get our dna to cleave again
        lab.length_sort()
        dna = lab.samples[0]

    # Using sequence check if an ending exists
    assert dna.sequence[0].strip() != "", "No ending state found"

    # Return if we are in the correct state at the end
    return not dna.has_sequence(IS_ODD_SEQUENCE)

if __name__ == "__main__":
    # Validate that the code works
    print("[*] Testing AT")
    assert automat_machine_even_b(PREFIX+A+T) == True
    print("[*] Testing ABT")
    assert automat_machine_even_b(PREFIX+A+B+T) == False
    print("[*] Testing ABBT")
    assert automat_machine_even_b(PREFIX+A+B+B+T) == True
    print("[*] Testing ABBBT")
    assert automat_machine_even_b(PREFIX+A+B+B+B+T) == False
    print("[*] Testing ABABABT")
    assert automat_machine_even_b(PREFIX+A+B+A+B+A+B+T) == False
    print("[*] Testing BABABABT")
    assert automat_machine_even_b(PREFIX+B+A+B+A+B+A+B+T) == False