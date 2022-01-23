from mml import *
import itertools
import random
import copy

SIZE = 6

variable_to_dna = {}
edges_to_node = {}

def parse_sat_formula(formula):
    """ Convert the 3Sat formula to a two dimensional array of of C's and their variables """
    formula = formula.replace("(", "").replace(")", "")
    sats = [[x for x in y.split("v")] for y in formula.split("^")]

    # Validate that the formula is really 3sat
    for C in sats:
        assert len(C) == 3, "Formula is not a 3Sat formula " + "v".join(C)

    return sats

def get_dna_value(var):
    """ Create a unique dna sequence for each variable """
    global variable_to_dna
    try:
        return copy.deepcopy(variable_to_dna[var])

    except KeyError:
        dna = DNA.generate_random(size=SIZE)
        dna.sequence[1] = ""
        variable_to_dna[var] = copy.deepcopy(dna)
        return dna

def get_all_dna_pairs(sats):
    """ Returns all variables with their negative counterpart """
    return [x for y in sats for x in [y, "~"+y]]

def get_nodes_from_edges(variables):
    global variable_to_dna
    nodes = []
    # Run on each variable
    for i, clean_var in enumerate(variables):
        # Each variable can also be negative
        for prefix in ["", "~"]:
            # Get the wanted variable
            var = prefix + clean_var
            half_len = len(variable_to_dna[var].sequence[0])//2

            # If we are on the first var, we need a short sequence to fill the start.
            # But it should also create a connection to the next dna
            if i == 0:
                base = variable_to_dna[var]
                sequence = base.sequence[0][:half_len]
                nodes.append(DNA("", DNA.get_pair_string(sequence)))

            # The last dna only needs half a sequence to fill the second half of the last dna
            if i == len(variables)-1:
                base = variable_to_dna[var]
                sequence = base.sequence[0][half_len:]
                nodes.append(DNA("", DNA.get_pair_string(sequence)))
                continue

            # All dna except the last should have a connection to the next variable (Positive or Negative)
            tail, head = variable_to_dna[var], variable_to_dna[variables[i+1]]
            sequence = tail.sequence[0][half_len:] + head.sequence[0][:half_len]
            d = DNA("", DNA.get_pair_string(sequence))
            nodes.append(d)

            # Add connection from current variable to the next one as negative
            tail, head = variable_to_dna[var], variable_to_dna["~"+variables[i+1]]
            sequence = tail.sequence[0][half_len:] + head.sequence[0][:half_len]
            nodes.append(DNA("", DNA.get_pair_string(sequence)))

    return nodes

def filter_by_sat(lab, sat):
    """
        Get a lab of samples, and check if it verifies the sat
        On each group of 3 of the sats, extract only the one that satisfies at least one of the variables
    """
    global variable_to_dna
    # Run on all Ci
    for x,y,z in sat:
        # Check if the samples satisfy the given 3 variables
        leftovers1 = Lab(lab.extract(variable_to_dna[x].sequence[0]))
        leftovers2 = Lab(leftovers1.extract(variable_to_dna[y].sequence[0]))
        leftovers2.extract(variable_to_dna[z].sequence[0])

        lab.samples += leftovers1.samples + leftovers2.samples

    return lab

def three_sat(formula):
    """ Calculate the given 3Sat formula using molecular programming """
    # Parse the formula to its variables
    print("[*] Parsing formula {}".format(formula))
    sat = parse_sat_formula(formula)
    variables = list(set([x.replace("~", "") for y in sat for x in y]))
    pairs = get_all_dna_pairs(variables)

    # Generate the needed dna's from the formula
    print("[*] Generating edges and nodes based on the formula")
    edges = [get_dna_value(x) for x in pairs]
    nodes = get_nodes_from_edges(variables)
    samples = edges + nodes

    # Create our lab from the given dna's for the formula
    print("[*] Creating the lab from the samples. Size: {}".format(len(samples)))
    lab = Lab(samples)

    # Amplify all samples using pcr. Split after running because we need to connect the edges and nodes
    print("[*] Amplifying samples using PCR")
    lab.amplify(['a', 'c', 't', 'g'], rounds=5)
    lab.split()

    # Merge all our dna samples
    print("[*] Merging our samples")
    lab.merge()

    # Search the sample for satisfying dna to the formula. Correct size and satisfiable
    print("[*] Filtering the unwanted data")
    lab.length_sort()
    expected_size = (SIZE*len(variables), SIZE*len(variables))
    lab.filter_by_length(expected_size)
    satisfiable = filter_by_sat(lab, sat).samples != []

    return satisfiable

formula = "(avbvc)^(~av~bvc)"
satisfiable = three_sat(formula)
print("{} is {} satisfiable".format(formula, "" if satisfiable else "not"))