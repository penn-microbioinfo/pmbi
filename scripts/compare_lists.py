import argparse


def compare_lists(list1, list2):
    """
    Compare two lists and return elements missing in each.

    Args:
        list1 (list): First list to compare.
        list2 (list): Second list to compare.

    Returns:
        dict: A dictionary with keys 'missing_in_list1' and 'missing_in_list2',
              containing elements missing in list1 and list2 respectively.
    """
    missing_in_list1 = [item for item in list2 if item not in list1]
    missing_in_list2 = [item for item in list1 if item not in list2]
    return {
        'missing_in_list1': missing_in_list1,
        'missing_in_list2': missing_in_list2
    }

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("list", nargs=2, help = "Lists to compare")
    args = parser.parse_args()
    list1, list2 = map(lambda l: [x.strip() for x in open(l, 'r').readlines()], args.list)

    print(compare_lists(list1, list2))

