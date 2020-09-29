def cleanup_header(lines_text):
    """
    This function cleans up some lines of text contained in a list of strings.
    :param lines_text: lines of text
    :return: header of the file "cleaner"
    Example:
    >>> cleanup_header(['# hello',"\t,# ciao"])
    ['hello', 'ciao']

    """
    lines_text = [s.replace("#", '') for s in lines_text]
    lines_text = [s.replace("\n", '') for s in lines_text]
    lines_text = [s.replace("\t", ' ') for s in lines_text]
    lines_text = [s.replace(" ", '') for s in lines_text]
    header = list(lines_text)
    return header

