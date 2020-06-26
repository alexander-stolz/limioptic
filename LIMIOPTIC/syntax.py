from PyQt5.QtCore import QRegExp

from PyQt5.QtGui import QColor, QTextCharFormat, QFont, QSyntaxHighlighter


def format(color, style=''):
    """ return QTextCharFormat """
    _color = QColor()
    _color.setNamedColor(color)

    _format = QTextCharFormat()
    _format.setForeground(_color)
    if 'bold' in style:     _format.setFontWeight(QFont.Bold)
    if 'italic' in style:   _format.setFontItalic(True)

    return _format


STYLES = {
    'keyword':      format('darkBlue', 'bold'),
    'function':     format('black', 'bold'),
    'brace':        format('black'),
    'defclass':     format('black', 'bold'),
    'string2':      format('darkGray'),
    'string':       format('darkMagenta'),
    'comment':      format('gray'),
    'numbers':      format('green'),
    'boldGreen':    format('green', 'bold')
}


class PythonHighlighter (QSyntaxHighlighter):
    # Python keywords
    keywords = [
        'and', 'assert', 'break', 'class', 'continue', 'def',
        'del', 'elif', 'else', 'except', 'exec', 'finally',
        'for', 'from', 'global', 'if', 'import', 'in',
        'is', 'lambda', 'not', 'or', 'pass', 'print',
        'raise', 'return', 'try', 'while', 'yield',
        'None', 'True', 'False', 'INPUT', "AddWaist()"]

    functions = [
        "AddThinLens", "AddMSA", "AddSlit", "AddESD",
        "AddVBFN", "AddAMSAcc", "AddBeam", "AddFNAccNeu",
        "ChangeBeamParameters", "AddGaussBeam", "AddBeamRandomGauss", "AddBeam3d",
        "AddBeamX", "AddQuadrupolRadFoc", "AddQuadrupolAxFoc", "AddAMSQPT", "AddBeamProfile", "AddFoil"]

    boldGreen = []

    braces = ['\{', '\}', '\(', '\)', '\[', '\]']

    def __init__(self, document):
        QSyntaxHighlighter.__init__(self, document)

        self.tri_single = (QRegExp("'''"), 1, STYLES['string2'])
        self.tri_double = (QRegExp('"""'), 2, STYLES['string2'])

        rules = []
        rules += [(r'\b%s\b' % w, 0, STYLES['keyword']) for w in PythonHighlighter.keywords]
        rules += [(r'\b%s' % o, 0, STYLES['function'])  for o in PythonHighlighter.functions]
        rules += [(r'\b%s' % o, 0, STYLES['boldGreen']) for o in PythonHighlighter.boldGreen]
        rules += [(r'%s' % b, 0, STYLES['brace'])       for b in PythonHighlighter.braces]

        rules += [
            (r'\b[+-]?[0-9]+[lL]?\b',            0, STYLES['numbers']),
            (r'\b[+-]?[0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?\b', 0, STYLES['numbers']),
            (r'"[^"\\]*(\\.[^"\\]*)*"',          0, STYLES['string']),
            (r"'[^'\\]*(\\.[^'\\]*)*'",          0, STYLES['string']),
            (r'#[^\n]*',                         0, STYLES['comment']),
        ]

        # QRegExp fuer jedes Muster
        self.rules = [(QRegExp(pat), index, fmt) for (pat, index, fmt) in rules]

    def highlightBlock(self, text):
        """ Apply syntax highlighting to the given block of text """

        for expression, nth, format in self.rules:
            index = expression.indexIn(text, 0)

            while index >= 0:
                index  = expression.pos(nth)
                length = len(expression.cap(nth))
                self.setFormat(index, length, format)
                index  = expression.indexIn(text, index + length)

        self.setCurrentBlockState(0)

        # multi-line strings
        in_multiline = self.match_multiline(text, *self.tri_single)
        if not in_multiline:
            in_multiline = self.match_multiline(text, *self.tri_double)

    def match_multiline(self, text, delimiter, in_state, style):
        """ Do highlighting of multi-line strings. ``delimiter`` should be a
        ``QRegExp`` for triple-single-quotes or triple-double-quotes, and
        ``in_state`` should be a unique integer to represent the corresponding
        state changes when inside those strings. Returns True if we're still
        inside a multi-line string when this function is finished. """
        # If inside triple-single quotes, start at 0
        if self.previousBlockState() == in_state:
            start = 0
            add = 0
        # Otherwise, look for the delimiter on this line
        else:
            start = delimiter.indexIn(text)
            # Move past this match
            add = delimiter.matchedLength()

        # As long as there's a delimiter match on this line...
        while start >= 0:
            # Look for the ending delimiter
            end = delimiter.indexIn(text, start + add)
            # Ending delimiter on this line?
            if end >= add:
                length = end - start + add + delimiter.matchedLength()
                self.setCurrentBlockState(0)
            # No; multi-line string
            else:
                self.setCurrentBlockState(in_state)
                length = len(text) - start + add
            # Apply formatting
            self.setFormat(start, length, style)
            # Look for the next match
            start = delimiter.indexIn(text, start + length)

        # Return True if still inside a multi-line string, False otherwise
        if self.currentBlockState() == in_state:
            return True
        else:
            return False
