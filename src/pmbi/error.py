def try_except(func):
    def wrapper(*args, **kwargs):
        print("meh")
        func(*args, **kwargs)
        print("meh2")
    return wrapper
