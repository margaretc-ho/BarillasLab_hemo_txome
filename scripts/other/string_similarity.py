import Levenshtein
import sys

def group_strings(strings, threshold, ignore_word=None):
    groups = []
    for string in strings:
        # Check if string contains the ignore word
        if ignore_word and ignore_word in string:
            index = string.index(ignore_word)
            # Omit the part of the string that includes the ignore word
            string = string[:index] + string[index+len(ignore_word):]

        # Check if string belongs to an existing group
        for group in groups:
            if any(Levenshtein.distance(string, s) <= threshold for s in group):
                group.append(string)
                break
        # If string does not belong to any existing group, create a new group
        else:
            groups.append([string])
    return groups

if __name__ == '__main__':
    # Get file path, threshold, and ignore word from command-line arguments
    file_path = sys.argv[1]
    threshold = int(sys.argv[2])
    ignore_word = sys.argv[3] if len(sys.argv) > 3 else None
    
    # Read strings from file
    with open(file_path, 'r') as f:
        strings = [line.strip() for line in f]

    # Group strings
    groups = group_strings(strings, threshold, ignore_word)

    # Print groups
    for i, group in enumerate(groups):
        print(f"Group {i+1}: {group}")