"""``python -m combobatch`` — the same entry point as the ``combobatch`` script.

The dispatcher launches its workers this way rather than by name, because the console
script is not guaranteed to be on ``PATH`` inside a pod or a virtualenv invoked by
absolute path, while ``sys.executable -m`` always is.
"""

from combobatch.cli.main import main

if __name__ == "__main__":
    raise SystemExit(main())
