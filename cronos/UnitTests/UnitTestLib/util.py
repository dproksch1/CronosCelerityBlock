# given an iterable of pairs return the key corresponding to the greatest value
def argmax(pairs):
    return max(pairs, key=lambda x: x[1])[0]

# given an iterable of values return the index of the greatest value
def argmax_index(values):
    return argmax(enumerate(values))

# given an iterable of keys and a function f, return the key with largest f(key)
def argmax_f(keys, f):
    return max(keys, key=f)

# given an iterable of pairs return the key corresponding to the smallest value
def argmin(pairs):
    return min(pairs, key=lambda x: x[1])[0]

# given an iterable of values return the index of the smallest value
def argmin_index(values):
    return argmin(enumerate(values))

# given an iterable of keys and a function f, return the key with smallest f(key)
def argmin_f(keys, f):
    return max(keys, key=f)