"""
Loading of the YAML configuration template.

The template supplies each parameter's default value only. When a block is
dropped onto the canvas it takes its own private copy of those defaults, so
editing one block's popup never affects another block of the same type, nor the
template itself.
"""

import os

import yaml


#The template lives at the repository root, next to config_HK.yaml
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DEFAULT_TEMPLATE = os.path.join(REPO_ROOT, "config.yaml")


class ConfigTemplate:
    """A loaded YAML template, queried per block for the parameters it uses."""

    def __init__(self, path=None):
        self.path = path or DEFAULT_TEMPLATE

        with open(self.path, "r") as stream:
            self.sections = yaml.safe_load(stream) or {}

    def parameter(self, section, key):
        """The full {Description, Type, Value, ...} record for one parameter."""
        return self.sections.get(section, {}).get(key)

    def describe(self, section, keys=None):
        """
        The parameter records for a block, in the order the block asked for
        them. Keys the template does not define are skipped rather than
        invented, so a block never shows a setting its script would ignore.
        """
        available = self.sections.get(section, {})
        keys = list(available) if keys is None else keys

        described = {}
        for key in keys:
            if key in available:
                described[key] = dict(available[key])

        return described

    def defaults(self, section, keys=None):
        """The default values only, ready to be copied into a block instance."""
        return {key: record.get("Value")
                for key, record in self.describe(section, keys).items()}


def coerce(value, declaredType):
    """
    Turn a popup's text back into the type the template declares.

    Only the value the user typed is converted; an unparseable entry is handed
    back untouched so the popup can report it rather than silently zeroing it.
    """
    if declaredType == "Boolean":
        if isinstance(value, bool):
            return value
        return str(value).strip().lower() in ("true", "1", "yes", "on")

    if declaredType == "Int":
        return int(str(value).strip())

    if declaredType == "Float":
        return float(str(value).strip())

    if declaredType == "BooleanList":
        if isinstance(value, list):
            return [bool(v) for v in value]

        text = str(value).strip().strip("[]")
        if not text:
            return []

        return [part.strip().lower() in ("true", "1", "yes", "on")
                for part in text.split(",")]

    return value
