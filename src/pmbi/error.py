def try_except(func):
    def wrapper(*args, **kwargs):
        try:
            res = func(*args, **kwargs)

            return res

        except Exception as e:
            raise e
    return wrapper

