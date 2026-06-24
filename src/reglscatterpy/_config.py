"""Package-level settings (currently just the default plot theme).

The theme controls the LIVE widget's card colours:

* ``"light"`` (default) — a white "figure card", portable and matching the
  static/exported HTML.
* ``"dark"`` — a dark card (light axes/legend), regardless of the host.
* ``"auto"`` (alias ``"system"``) — follow the notebook: a dark card only when
  the host theme (VS Code / JupyterLab) is dark, otherwise white.

Set it globally once::

    import reglscatterpy as rs
    rs.set_theme("auto")            # every plot now follows the notebook theme

or per call: ``rs.scatterplot(..., theme="dark")`` (overrides the global).
"""

_VALID = ("light", "dark", "auto", "system")
_THEME = "light"


def _normalize(value: str) -> str:
    if not isinstance(value, str) or value.lower() not in _VALID:
        raise ValueError(
            f"theme must be one of {_VALID!r}, got {value!r}. "
            "Use 'light' (white card), 'dark', or 'auto'/'system' (match the host)."
        )
    v = value.lower()
    return "auto" if v == "system" else v


def set_theme(value: str) -> str:
    """Set the default theme for all subsequent plots. Returns the value set."""
    global _THEME
    _THEME = _normalize(value)
    return _THEME


def get_theme() -> str:
    """The current default theme (``"light"``, ``"dark"`` or ``"auto"``)."""
    return _THEME


def resolve_theme(value=None) -> str:
    """Per-call resolution: an explicit ``value`` wins, else the global default."""
    return _normalize(value) if value is not None else _THEME
