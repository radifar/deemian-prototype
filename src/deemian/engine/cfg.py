from lark import Lark


def parser(text):
    """Generate parser using Lark and grammar then parse the input text"""
    lark_parser = Lark.open("grammar.lark", rel_to=__file__, parser="lalr")
    return lark_parser.parse(text)
