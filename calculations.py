from mml import *
import itertools
import random
import copy

SIZE = 6

variable_to_dna = {}
edges_to_node = {}

def parse_sat_formula(formula):
    formula = formula.replace("(", "").replace(")", "")
    sats = [[x for x in y.split("v")] for y in formula.split("^")]
    return sats

def get_dna_value(var):
    global variable_to_dna
    try:
        return copy.deepcopy(variable_to_dna[var])

    except KeyError:
        dna = DNA.generate_random(size=SIZE)
        dna.sequence[1] = ""
        variable_to_dna[var] = copy.deepcopy(dna)
        return dna

def get_connecting_node(var1, var2):
    # TODO: Think about setting the start and finish
    global edges_to_node

    try:
        return copy.deepcopy(edges_to_node[str(var1)+str(var2)])

    except KeyError:
        # Get corresponding edges
        edge1 = get_dna_value(var1)
        edge2 = get_dna_value(var2)

        # Calculate the connecting value of both edges
        half_len = len(edge1.sequence[0])//2
        sequence = edge1.sequence[0][half_len:] + edge2.sequence[0][:half_len]
        opposite = DNA.get_pair_string(sequence)

        # Save for future use
        dna = DNA("", opposite)
        edges_to_node[str(var1)+str(var2)] = copy.deepcopy(dna)
        return dna

def get_tail_nodes(var1, var2):
    dna = get_connecting_node(var1, var2)
    half_len = len(dna.sequence[1])//2
    return dna.cleave(("", dna.sequence[1][half_len:]))


def get_all_dna_pairs(sats):
    return [x for y in sats for x in [y, "~"+y]]

def get_nodes_from_edges(variables):
    global variable_to_dna
    nodes = []
    for i, clean_var in enumerate(variables):
        for prefix in ["", "~"]:
            var = prefix + clean_var
            half_len = len(variable_to_dna[var].sequence[0])//2

            if i == 0:
                base = variable_to_dna[var]
                sequence = base.sequence[0][:half_len]
                nodes.append(DNA("", DNA.get_pair_string(sequence)))

            if i == len(variables)-1:
                base = variable_to_dna[var]
                sequence = base.sequence[0][half_len:]
                nodes.append(DNA("", DNA.get_pair_string(sequence)))
                continue

            tail, head = variable_to_dna[var], variable_to_dna[variables[i+1]]
            sequence = tail.sequence[0][half_len:] + head.sequence[0][:half_len]
            d = DNA("", DNA.get_pair_string(sequence))
            nodes.append(d)
            print("TAIL|HEAD|SEQ", tail, head, d)

            tail, head = variable_to_dna[var], variable_to_dna["~"+variables[i+1]]
            sequence = tail.sequence[0][half_len:] + head.sequence[0][:half_len]
            nodes.append(DNA("", DNA.get_pair_string(sequence)))

    return nodes

def filter_by_sat(lab, sat):
    global variable_to_dna
    """ Get a lab of samples, and check if it verifies sat """
    for x,y,z in sat:
        leftovers1 = Lab(lab.extract(variable_to_dna[x].sequence[0]))
        leftovers2 = Lab(leftovers1.extract(variable_to_dna[y].sequence[0]))
        leftovers2.extract(variable_to_dna[z].sequence[0])

        lab.samples += leftovers1.samples + leftovers2.samples

    return lab

def three_sat(formula):
    sat = parse_sat_formula(formula)
    variables = list(set([x.replace("~", "") for y in sat for x in y]))
    pairs = get_all_dna_pairs(variables)

    edges = [get_dna_value(x) for x in pairs]
    nodes = get_nodes_from_edges(variables)
    samples = edges + nodes

    # Get a lot of each sample (2^n)
    for i in range(5):
        samples = samples + samples

    lab = Lab(samples)

    print(lab.samples)
    for i in range(2):
        lab.merge()

    expected_size = (SIZE*len(variables), SIZE*len(variables))
    lab.filter_by_length(expected_size)

    left = filter_by_sat(lab, sat)
    print("------------\n", left.samples)


def autamat_machine(automat):
    pass

three_sat("(avbvc)^(~av~bv~c)")