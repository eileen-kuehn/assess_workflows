class MulticoreResult(dict):
    def __iadd__(self, other):
        to_check = [(self, other,)]
        while to_check:
            reference, other_reference = to_check.pop(0)
            try:
                for key, value in other_reference.items():
                    if key not in reference:
                        reference[key] = value
                    else:
                        to_check.append((reference.get(key), value,))
            except AttributeError:
                reference.extend(other_reference)
        return self


def multicore_factor(problem_size, min_size=100):
    if problem_size < min_size * 2:
        return 1
    result = 10
    while problem_size / result > min_size:
        result += 1
    return result
