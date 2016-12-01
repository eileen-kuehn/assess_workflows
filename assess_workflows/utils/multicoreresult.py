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
