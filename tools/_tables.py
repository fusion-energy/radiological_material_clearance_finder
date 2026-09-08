"""Shared helpers for parsing regulatory tables out of official sources.

Each source writes its numbers differently and none of them fail loudly when
misread, so the parsing of a value is kept in one place with an assertion that
the result looks like a regulatory limit.
"""
from __future__ import annotations

import html
import re

# legislation.gov.uk marks exponents with <Superior>, so 10^2 is
# "10<Superior>2</Superior>". Flattening tags first would turn that into "102".
_SUPERIOR = re.compile(r"<Superior[^>]*>(.*?)</Superior>", re.S)
_TAG = re.compile(r"<[^>]+>")
_TABULAR = re.compile(r"<Tabular[^>]*>(.*?)</Tabular>", re.S)
_NUMBER = re.compile(r"<Number>(.*?)</Number>", re.S)
_TITLE = re.compile(r"<Title>(.*?)</Title>", re.S)
_ROW = re.compile(r"<tr[^>]*>(.*?)</tr>", re.S)
_CELL = re.compile(r"<t[dh][^>]*>(.*?)</t[dh]>", re.S)

# 1 x 10^6, 1 × 10^6, 4 x 10^5, or a bare 10^-2
_TIMES_POWER = re.compile(r"^([\d.,]+)\s*[x×]\s*10\^(-?\d+)$")
_POWER = re.compile(r"^10\^(-?\d+)$")


def clean(fragment: str) -> str:
    """Flatten a markup fragment to text, preserving exponents as ``^n``."""
    text = _SUPERIOR.sub(lambda m: "^" + _TAG.sub("", m.group(1)).strip(), fragment)
    text = _TAG.sub("", text)
    text = html.unescape(text)
    # Non-breaking spaces are used as thousands separators, and an en dash is
    # used for a negative exponent in some sources.
    text = text.replace(" ", " ").replace("–", "-").replace("−", "-")
    return re.sub(r"\s+", " ", text).strip()


def tabulars(document: str) -> list[tuple[str, str, str]]:
    """Return ``(number, title, body)`` for each ``<Tabular>`` in the document."""
    out = []
    for match in _TABULAR.finditer(document):
        body = match.group(1)
        number = _NUMBER.search(body)
        title = _TITLE.search(body)
        out.append(
            (
                clean(number.group(1)) if number else "",
                clean(title.group(1)) if title else "",
                body,
            )
        )
    return out


def rows(table_body: str) -> list[list[str]]:
    """Return the cleaned cells of each row in a table."""
    return [[clean(cell) for cell in _CELL.findall(row)] for row in _ROW.findall(table_body)]


def parse_value(text: str, *, decimal_comma: bool = False) -> float | None:
    """Read a regulatory limit, returning ``None`` when the cell holds no number.

    Handles the notations these sources actually use: ``10^2``, ``4 x 10^5``,
    ``1 × 10^6``, plain decimals, comma decimal separators and space or
    non-breaking-space thousands separators.

    Args:
        text: The cleaned cell text.
        decimal_comma: Read ``0,1`` as 0.1, as EUR-Lex writes it.

    Raises:
        ValueError: If the cell holds something numeric-looking that does not
            parse, since silently skipping it would drop a real limit.
    """
    value = text.strip()
    if not value or value.lower() in {"no limit", "none", "-", "—", "not applicable"}:
        return None

    if decimal_comma:
        # "1 000,5" -> "1000.5". Thousands separators are spaces here, so they
        # can be removed before the comma becomes a point.
        value = value.replace(" ", "").replace(",", ".")
    else:
        value = re.sub(r"(?<=\d) (?=\d\d\d\b)", "", value)

    # The German source writes powers as "1 E-1" with a space before the E.
    value = re.sub(r"(?<=[\d.])\s*[Ee]\s*([+-]?\d+)$", r"e\1", value)

    match = _POWER.match(value)
    if match:
        return 10.0 ** int(match.group(1))
    match = _TIMES_POWER.match(value)
    if match:
        return float(match.group(1)) * 10.0 ** int(match.group(2))
    try:
        return float(value)
    except ValueError:
        if re.search(r"\d", value):
            raise ValueError(f"cannot read {text!r} as a limit") from None
        return None


def check_regulatory(values: dict[str, float], label: str) -> None:
    """Assert every limit is a single significant figure times a power of ten.

    These regulations tabulate values of that shape without exception, so a
    value like 102 (which is what 10^2 becomes when superscript markup is
    flattened) means the parse went wrong. Catching it here is the difference
    between a wrong table and a loud failure.
    """
    import math

    bad = []
    for name, value in values.items():
        if value <= 0:
            bad.append((name, value))
            continue
        exponent = math.floor(math.log10(value))
        mantissa = value / 10.0 ** exponent
        if abs(mantissa - round(mantissa)) > 1e-9 or not 1 <= round(mantissa) <= 9:
            bad.append((name, value))
    if bad:
        raise ValueError(
            f"{label}: {len(bad)} limit(s) are not one significant figure times a "
            f"power of ten, which usually means superscript markup was flattened: "
            f"{bad[:8]}"
        )
