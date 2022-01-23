from itertools import chain
import random

PCR_FAIL_CHANCE = 0.0001

class Lab():
    def __init__(self, samples):
        self.samples = samples

    def polymerase(self):
        """ Fill the samples based on the longer existing sequence """
        for dna in self.samples:
            dna.fill()

        self.merge()

    def split(self):
        """ Split the dna into the two sequences """
        self.samples = list(chain.from_iterable(x.split() for x in self.samples))

    def merge(self):
        """ Try to merge two samples based on their sequence """
        for dna in self.samples:
            # Randomize second dna so each run will be unique and not connect the same dna's
            indexs = list(range(len(self.samples)))
            random.shuffle(indexs)
            for i in indexs:
                # If second dna still exists try to merge it
                connector = self.samples[i]
                if connector.removed == True:
                    continue

                dna.merge(connector)

        # Remove all merged dna's
        self.samples = [x for x in self.samples if x.removed == False]

    def extract(self, sequence):
        """ Remove from the lab all samples that don't have our sequence. Returns whats left """
        leftovers = [x for x in self.samples if not x.has_sequence(sequence)]
        self.samples = [x for x in self.samples if x.has_sequence(sequence)]

        return leftovers

    def length_sort(self):
        """ Sort the samples by len their len """
        self.samples.sort(key=lambda x: len(x.sequence[0]) + len(x.sequence[1]), reverse=True)

    def filter_by_length(self, expected):
        self.samples = [x for x in self.samples if len(x.sequence[0]) == expected[0] and len(x.sequence[1]) == expected[1]]

    def amplify(self, primers, rounds=10):
        """ Amplify a gene based on pcr """
        for i in range(rounds):
            # Split all dna's to allow connection of primers
            self.split()

            # For each dna, add the primers
            for dna in self.samples:
                for primer in primers:
                    if dna.add_primer(primer):
                        break

            # Extend the primers
            for dna in self.samples:
                dna.fill()

            # Merge remaining dna and get rid of what couldn't connect
            self.merge()
            self.cleanup()

    def cleanup(self):
        """ Remove all not full dna """
        self.samples = [x for x in self.samples if len(x.sequence[0].strip()) != 0 and len(x.sequence[1].strip()) != 0]

    def cleave(self, enzyme):
        """ Cut a dna sample based on an enzym """
        for dna in self.samples:
            # Try to cut the dna by the enzym, if it succeeded add the second half to the dna samples
            other = dna.cleave(enzyme)
            if other:
                self.samples.append(other)


class DNA():
    def __init__(self, sequence1, sequence2):
        self.sequence = [sequence1, sequence2]
        self.removed = False

    def has_sequence(self, sequence):
        """ Return if the sequence is inside the dna """
        return sequence in self.sequence[0] or sequence in self.sequence[1]

    def split(self):
        """ Split the dna strand into two dna's (Cuts the sequence) """
        return DNA(self.sequence[0], ""), DNA("", self.sequence[1])

    def merge(self, other):
        """ Check if we can merge with the other dna sequence """
        my_first_len = len(self.sequence[0].strip())
        my_second_len = len(self.sequence[1].strip())
        other_first_len = len(other.sequence[0].strip())
        other_second_len = len(other.sequence[1].strip())

        # Only fill objects with primers and not empty ones
        if my_first_len == my_second_len or other_second_len == other_first_len:
            return False

        # If we need to add to our second sequnce
        # aaaaaa
        # ttt       tttcccc
        my_larger_index = 0 if my_first_len > my_second_len else 1
        my_shorter_index = 1 if my_first_len > my_second_len else 0
        my_shorter_len = len(self.sequence[my_shorter_index])

        # Get the difference
        my_diff = self.sequence[my_larger_index][my_shorter_len:]

        # Make sure we have enough empty space
        if other.sequence[my_larger_index][:len(my_diff)].strip() != "":
            return False

        # Make sure that we match
        other_diff = other.sequence[my_shorter_index][:len(my_diff)]
        if not my_diff.startswith(self.get_pair_string(other_diff)):
            return False

        # Append the sequences
        self.sequence[my_shorter_index] += other.sequence[my_shorter_index].strip()
        self.sequence[my_larger_index]  += other.sequence[my_larger_index].strip()

        # Set the second dna as used so we wont use it again
        other.removed = True
        return True

    def fill(self):
        """ Fill the rest of the sequence based on the other sequence """
        first_len = len(self.sequence[0])
        second_len = len(self.sequence[1])

        # Only fill objects with primers and not empty ones
        if first_len == 0 or second_len == 0:
            return

        # If they are the same length we don't need to check anything
        if first_len == second_len:
            pass

        # Fill the second sequence
        elif first_len > second_len:
            self.sequence[1] += " " * (first_len-second_len)

        # Fill the first sequence
        else:
            self.sequence[0] += " " * (second_len-first_len)

        # Convert to list to allow item assignment
        self.sequence[0] = list(self.sequence[0])
        self.sequence[1] = list(self.sequence[1])

        # Fill with the pair of the base
        for i in range(len(self.sequence[0])):
            # Fill the pair with a small chance of failing
            if self.sequence[0][i] == " ":
                self.sequence[0][i] = self.get_pair_base(self.sequence[1][i])
                if random.random() < PCR_FAIL_CHANCE:
                    self.sequence[0][i], self.sequence[1][i] = self.sequence[1][i], self.sequence[0][i]

            elif self.sequence[1][i] == " ":
                self.sequence[1][i] = self.get_pair_base(self.sequence[0][i])
                if random.random() < PCR_FAIL_CHANCE:
                    self.sequence[0][i], self.sequence[1][i] = self.sequence[1][i], self.sequence[0][i]

        # Convert lists back to strings
        self.sequence[0] = "".join(self.sequence[0])
        self.sequence[1] = "".join(self.sequence[1])

    def add_primer(self, primer):
        """ Add the primer to the sequences. WILL ONLY WORK AFTER SPLIT """
        pair = self.get_pair_string(primer)
        # Try to see if the primer matches the first/second sequence
        if pair in self.sequence[0]:
            index = self.sequence[0].index(pair)
            self.sequence[1] = " " * index + primer

        elif pair in self.sequence[1]:
            index = self.sequence[1].index(pair)
            self.sequence[0] = " " * index + primer

        else:
            return False

        return True

    def cleave(self, enzyme):
        """ Cut the dna and sequence by the given enzym. If not found return none and do nothing """
        # Find the location of the enzyme sequence
        try:
            location1 = self.sequence[0].index(enzyme[0])
            location2 = self.sequence[1].index(enzyme[1])
        except ValueError:
            return None

        # If the cut location was found, cut it and return two dnas
        other = DNA(self.sequence[0][location1:], self.sequence[1][location2:])
        self.sequence[0] = self.sequence[0][:location1]
        self.sequence[1] = self.sequence[1][:location2]
        return self, other

    def __str__(self):
        return "(" + self.sequence[0] + "|" + self.sequence[1] + ")"

    def __repr__(self):
        return str(self)

    @staticmethod
    def generate_random(size, single_sequence=False):
        """ Return a randomly generated dna """
        first_sequence = "".join(random.choices(['a', 't', 'c', 'g'], k=size))
        second_sequence = ""
        if not single_sequence:
            second_sequnce = DNA.get_pair_string(first_sequence)

        return DNA(first_sequence, second_sequence)

    @staticmethod
    def get_pair_base(char):
        """ Return the pair for the base """
        opposite = {
            'a': 't',
            't': 'a',
            'c': 'g',
            'g': 'c'
        }
        return opposite[char]

    @staticmethod
    def get_pair_string(str):
        """ Return the pair for the base """
        return "".join([DNA.get_pair_base(x) for x in str])