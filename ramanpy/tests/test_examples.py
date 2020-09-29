from ..generic_fit_class import GenericFit

from ..tools import cleanup_header


def test_cleanup_header():
    data = ['# hello', "\t,# ciao"]
    actual = cleanup_header(data)
    expected = ['hello', ',ciao']
    assert actual == expected


def test_bkg_model_selection():
    # check if bkg model selection returns a tuple, ugly test, but ok...
    data = 'linear'
    actual = GenericFit.choose_bkg_model(poly_type=data)
    expected = tuple

    assert isinstance(actual, expected)
