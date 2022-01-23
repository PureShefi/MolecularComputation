from mml import *
import itertools
import random
import copy

SIZE = 8

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

def get_dna_value(var, create=False):
    """ Create a unique dna sequence for each variable """
    global variable_to_dna
    if not create:
        return copy.deepcopy(variable_to_dna[var])

    # If we don't care about creating run in try clause which genereates dna if it doesn't exist
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
    nodes = []
    # Run on each variable
    for i, clean_var in enumerate(variables):
        # Each variable can also be negative
        for prefix in ["", "~"]:
            # Get the wanted variable
            var = prefix + clean_var
            half_len = len(get_dna_value(var).sequence[0])//2

            # If we are on the first var, we need a short sequence to fill the start.
            # But it should also create a connection to the next dna
            if i == 0:
                base = get_dna_value(var)
                sequence = base.sequence[0][:half_len]
                nodes.append(DNA("", DNA.get_pair_string(sequence)))

            # The last dna only needs half a sequence to fill the second half of the last dna
            if i == len(variables)-1:
                base = get_dna_value(var)
                sequence = base.sequence[0][half_len:]
                nodes.append(DNA("", DNA.get_pair_string(sequence)))
                continue

            # All dna except the last should have a connection to the next variable (Positive or Negative)
            tail, head = get_dna_value(var), get_dna_value(variables[i+1])
            sequence = tail.sequence[0][half_len:] + head.sequence[0][:half_len]
            d = DNA("", DNA.get_pair_string(sequence))
            nodes.append(d)

            # Add connection from current variable to the next one as negative
            tail, head = get_dna_value(var), get_dna_value("~"+variables[i+1])
            sequence = tail.sequence[0][half_len:] + head.sequence[0][:half_len]
            nodes.append(DNA("", DNA.get_pair_string(sequence)))

    return nodes

def filter_by_sat(lab, sat):
    """
        Get a lab of samples, and check if it verifies the sat
        On each group of 3 of the sats, extract only the one that satisfies at least one of the variables
    """
    # Run on all Ci
    for x,y,z in sat:
        # Check if the samples satisfy the given 3 variables
        leftovers1 = Lab(lab.extract(get_dna_value(x).sequence[0]))
        leftovers2 = Lab(leftovers1.extract(get_dna_value(y).sequence[0]))
        leftovers2.extract(get_dna_value(z).sequence[0])

        lab.samples += leftovers1.samples + leftovers2.samples

    return lab

def sanitize_random_edges(edges):
    """ Make sure there are no two similar edges """
    for x in edges:
        for y in edges:
            v_x = x.sequence[0]
            v_y = y.sequence[0]
            half_len = len(v_x)//2
            if x != y:
                assert v_x[:half_len] != v_y[:half_len], "x:y - {}:{}".format(x, y)
                assert v_x[half_len:] != v_y[half_len:], "x:y - {}:{}".format(x, y)

def three_sat(formula):
    """ Calculate the given 3Sat formula using molecular programming """
    # Parse the formula to its variables
    print("[*] Parsing formula {}".format(formula))
    sat = parse_sat_formula(formula)
    variables = list(set([x.replace("~", "") for y in sat for x in y]))
    pairs = get_all_dna_pairs(variables)

    # Generate the needed dna's from the formula
    print("[*] Generating edges and nodes based on the formula")
    edges = [get_dna_value(x, True) for x in pairs]
    nodes = get_nodes_from_edges(variables)
    samples = copy.deepcopy(edges) + nodes

    # Make sure are generated dna's are unique enough
    sanitize_random_edges(edges)

    # Create our lab from the given dna's for the formula
    print("[*] Creating the lab from the samples. Size: {}".format(len(samples)))
    lab = Lab(samples)

    # Amplify all samples using pcr. Split after running because we need to connect the edges and nodes
    print("[*] Amplifying samples using PCR")
    lab.amplify(['a', 'c', 't', 'g'], rounds=6)
    lab.split()

    # Merge all our dna samples
    print("[*] Merging our samples")
    lab.merge(log=True)

    # Search the sample for satisfying dna to the formula. Correct size and satisfiable
    print("[*] Filtering the unwanted data")
    lab.length_sort()
    expected_size = (SIZE*len(variables), SIZE*len(variables))
    lab.filter_by_length(expected_size)
    tmp = copy.deepcopy(lab)
    satisfiable = filter_by_sat(lab, sat).samples != []


    # Log result and return its value
    print("[RESULT] {} is {}satisfiable".format(formula, "" if satisfiable else "not "))
    return satisfiable

if __name__ == "__main__":
    formula = "(avbvc)^(~av~bvc)"
    satisfiable = three_sat(formula)

    formula = "(avava)^(~av~av~a)"
    satisfiable = three_sat(formula)
