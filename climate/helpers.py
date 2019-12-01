import ast


def evaluate(val):
    try:
        val = ast.literal_eval(val)
    except ValueError:
        pass
    return val