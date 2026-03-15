# -*- coding: utf-8 -*-
"""
i18n.py — Lightweight internationalization loader for gstd TUI.

Usage:
    from genome_standardizer.tui.i18n import t, set_lang, get_i18n

    # Get a string
    label = t("input.gff_label")

    # Get a string with format placeholders
    msg = t("log.done", elapsed=3.14)

    # Switch language at runtime
    set_lang("zh")

Supported languages: "en" (default), "zh"
Locale files: genome_standardizer/locales/{lang}.json
"""

import json
import os
from pathlib import Path
from typing import Optional

# ── Constants ────────────────────────────────────────────────────────────────

LOCALES_DIR     = Path(__file__).parent.parent / "locales"
SUPPORTED_LANGS = ("en", "zh")
DEFAULT_LANG    = "en"


# ── Core class ───────────────────────────────────────────────────────────────

class I18n:
    """
    Internationalization helper.

    - Loads locale JSON files from genome_standardizer/locales/.
    - Always keeps English as a fallback: if a key is missing in the
      active language, the English value is returned instead.
    - Supports dot-separated nested key access: t("section.key").
    - Supports str.format() placeholders: t("log.done", elapsed=1.5).
    """

    def __init__(self, lang: str = DEFAULT_LANG) -> None:
        self._lang:    str  = DEFAULT_LANG
        self._en_data: dict = {}
        self._data:    dict = {}

        self._load_lang("en")           # always load English first (fallback)
        if lang != "en":
            self.set_lang(lang)

    # ── Public API ────────────────────────────────────────────────────────────

    @property
    def lang(self) -> str:
        """Currently active language code."""
        return self._lang

    def set_lang(self, lang: str) -> None:
        """
        Switch the active language at runtime.

        Args:
            lang: Language code, must be in SUPPORTED_LANGS.

        Raises:
            ValueError: If lang is not supported.
            FileNotFoundError: If the locale JSON file does not exist.
        """
        if lang not in SUPPORTED_LANGS:
            raise ValueError(
                f"Unsupported language '{lang}'. "
                f"Choose from: {', '.join(SUPPORTED_LANGS)}"
            )
        self._load_lang(lang)
        self._lang = lang

    def t(self, key: str, **kwargs) -> str:
        """
        Return the localized string for *key*.

        Lookup order:
          1. Active language data
          2. English fallback data
          3. Placeholder string "[MISSING:<key>]"

        Args:
            key:    Dot-separated key path, e.g. "input.gff_label".
            **kwargs: Optional format arguments, e.g. elapsed=3.14.

        Returns:
            Formatted localized string.
        """
        value = self._get_nested(self._data, key)
        if value is None:
            value = self._get_nested(self._en_data, key)
        if value is None:
            return f"[MISSING:{key}]"

        if kwargs:
            try:
                return value.format(**kwargs)
            except (KeyError, ValueError):
                return value   # return raw string if format fails

        return value

    # ── Private helpers ───────────────────────────────────────────────────────

    def _load_lang(self, lang: str) -> None:
        path = LOCALES_DIR / f"{lang}.json"
        if not path.exists():
            raise FileNotFoundError(
                f"Locale file not found: {path}\n"
                f"Expected at: {LOCALES_DIR / (lang + '.json')}"
            )
        try:
            with open(path, encoding="utf-8") as fh:
                data = json.load(fh)
        except json.JSONDecodeError as exc:
            raise ValueError(f"Invalid JSON in {path}: {exc}") from exc

        if lang == "en":
            self._en_data = data
        self._data = data

    @staticmethod
    def _get_nested(data: dict, key: str) -> Optional[str]:
        """Traverse nested dict using a dot-separated key path."""
        parts = key.split(".")
        node  = data
        for part in parts:
            if not isinstance(node, dict) or part not in node:
                return None
            node = node[part]
        return node if isinstance(node, str) else None


# ── Module-level singleton ────────────────────────────────────────────────────
# The TUI app holds one I18n instance and passes it around, but these
# module-level helpers are provided for convenience (e.g. in tests).

_instance: I18n = I18n(DEFAULT_LANG)


def get_i18n() -> I18n:
    """Return the shared module-level I18n singleton."""
    return _instance


def set_lang(lang: str) -> None:
    """Switch the module-level singleton to *lang*."""
    _instance.set_lang(lang)


def t(key: str, **kwargs) -> str:
    """Translate *key* using the module-level singleton."""
    return _instance.t(key, **kwargs)
