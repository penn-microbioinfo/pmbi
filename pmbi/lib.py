def dict_try_insert(d: dict, key, value):
    if not key in d:
        d[key] = value
    return d
